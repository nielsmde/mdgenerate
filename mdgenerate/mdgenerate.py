"""
Python module to generate GROMACS simulation input
and submit them to SLURM.
"""

import re
import os
import shutil
import time
import yaml

import jinja2

env = jinja2.Environment(loader=jinja2.PackageLoader('mdgenerate', 'templates'))


def time_to_ps(tstr):
    """
    Convert a time with unit to a float in pico seconds.
    Supported units: fs, ps, ns, us, ms
    """
    prefactors = ['', 'm', 'u', 'n', 'p', 'f']
    m = re.match('([\d.]+)([{}]?)s'.format(''.join(prefactors)), tstr)
    if m is None:
        raise ValueError('Could not parse time: {}'.format(tstr))
    val, prefactor = m.groups()
    print(prefactor)
    decade = -3 * prefactors.index(prefactor) + 12
    return float(val) * 10**decade


def locate_template(template):
    """Locate the template file of a given name."""
    *base, ext = template.split('.')
    if ext != 'yaml':
        template += '.yaml'
    return os.path.abspath(os.path.join(os.environ['HOME'], '.mdgenerate/templates', template))


class MetaDict(dict):
    """
    Dictionary of simulation metadata with special update method that merges mdp-options
    and raises meanigful KeyErrors.
    """

    time_keys = ['time', 'timestep']

    def _timestamp(self):
        self['timestamp'] = time.strftime('%D.%m.%Y %H:%M')

    @property
    def as_mdp(self):
        self._timestamp()
        return env.get_template('template.mdp').render(**self)

    @property
    def as_slurmstr(self):
        self._timestamp()
        return env.get_template('template.slurm').render(**self)

    def load_yaml(yaml_file):
        """Load a yaml file and process its template."""
        with open(yaml_file) as f:
            yaml_dict = MetaDict(yaml.load(f.read()), yaml_file)

        if 'template' in yaml_dict:
            template_dict = MetaDict(locate_template(yaml_dict['template']))
            template_dict.update(yaml_dict)
            yaml_dict = template_dict

        return yaml_dict

    def update(self, update_dict):
        """Update the dictionary, by updating 'mdp-options' recursively."""
        dict_items = []
        for key in list(update_dict.keys()):
            if isinstance(update_dict[key], dict):
                dict_items.append((key, update_dict.pop(key)))
        super().update(update_dict)

        for key, val in dict_items:
            if key in self:
                self[key].update(val)
            else:
                self[key] = val

    def __init__(self, yaml_file):
        self.directory, self.filename = os.path.split(yaml_file)
        with open(yaml_file) as f:
            yaml_dict = yaml.load(f.read())

        if 'template' in yaml_dict:
            if isinstance(yaml_dict['template'], list):
                templates = yaml_dict['template']
                template_dict = MetaDict(locate_template(templates[0]))
                for template in templates:
                    with open(locate_template(template)) as f:
                        template_dict.update(yaml.load(f.read()))
            else:
                template_dict = MetaDict(locate_template(yaml_dict['template']))
            template_dict.update(yaml_dict)
            yaml_dict = template_dict

        super().update(yaml_dict)
        for key in self.time_keys:
            if key in self:
                val = self[key]
                if isinstance(val , str):
                    self[key] = time_to_ps(val)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            raise KeyError('Mandatory key {} not found in file: {}'.format(key, self.filename))


def generate_directory(meta, override=False):
    """
    Generate the simulation directory based on a MetaDict.

    Args:
        meta (MetaDict): The metadata of the simulation.
        override (bool, opt.): Override an exsiting directory.
    """
    _, dirs, files = next(os.walk(meta.directory))
    if dirs or files.remove(meta.filename):
        if not override:
            raise FileExistsError('Simulation directory is not empty: {}'.format(meta.directory))
        else:
            print('Simulation directory is not empty, but override=True')

    indir = os.path.join(meta.directory, meta.indir)
    outdir = os.path.join(meta.directory, meta.outdir)
    os.mkdir(indir)
    os.mkdir(outdir)
    shutil.copyfile(meta.top, indir)
    shutil.copyfile(meta.gro, indir)
