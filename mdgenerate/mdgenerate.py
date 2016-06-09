"""
Python module to generate GROMACS simulation input
and submit them to SLURM.
"""

import re
import os
import shutil
import time
import subprocess

# from subprocess import run, PIPE, STDOUT

import yaml
import jinja2

GMX_VERSION = '5.1'

DEFAULT_KEYS = {
    'mdp-file': 'mdpin.mdp',
    'slurm-file': 'slurm.sh',
    'tpr-file': 'topol.tpr',

}

TEMPLATE_DIR = os.path.join(os.environ['HOME'], '.mdgenerate/templates')

env = jinja2.Environment(loader=jinja2.PackageLoader('mdgenerate', 'templates'))


class GMXError(OSError):
    pass


def run_gmx(command):
    process = subprocess.Popen(
        ['/bin/bash', '-l'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    stdin = """
        module purge gromacs
        module load gromacs/{gmx}
        {cmd}
    """.format(gmx=GMX_VERSION, cmd=command)
    stdin, stderr = process.communicate(stdin)
    process.terminate()
    if process.returncode is not 0:
        raise GMXError('Error in gmx command: {}'.format(command), stderr)
    return stdin, stderr


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
    return os.path.abspath(os.path.join(TEMPLATE_DIR, template))


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
    def as_slurm(self):
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

    def _init_with_yaml(self, yaml_file):
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
                if isinstance(val, str):
                    self[key] = time_to_ps(val)

        for key, val in DEFAULT_KEYS.items():
            self.setdefault(key, val)

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], str):
            self._init_with_yaml(args[0])
        else:
            super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            raise KeyError('Mandatory key {} not found in file: {}'.format(key, self.filename))


def generate_files(meta, override=False):
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

    indir = os.path.join(meta.directory, meta.get('indir', ''))
    outdir = os.path.join(meta.directory, meta.get('outdir', ''))
    os.mkdir(indir)
    os.mkdir(outdir)

    if meta.get('copy-topology', False):
        for f in meta['topology-files']:
            shutil.copyfile(f, indir)

    with open(meta['mdp-file'], 'w') as f:
        f.write(meta.as_mdp)


def generate_slurm(meta):
    slurm_file = os.path.join(meta.directory, meta.get('slurm-file', 'slurm.sh'))
    with open(slurm_file, 'w') as f:
        f.write(meta.as_slurm)


def grompp(meta, generate=True):
    if generate:
        generate_files(meta)

    indir = os.path.join(meta['directory'], meta['indir'])
    args = [
        'gmx', 'grompp',
        '-f', os.path.join(indir, meta['mdp-file']),
    ]
    for fname in meta['topology-files']:
        ext = fname.split('.')[-1]
        if ext in ['gro', 'g96', 'pdb', 'brk', 'ent', 'esp', 'tpr']:
            args += ['-c', fname]
        elif ext == 'ndx':
            args += ['-n', fname]
        elif ext == 'top':
            args += ['-p', fname]

    assert '-c' in args, 'No structure file specified.'
    assert '-p' in args, 'No topology file specified.'

    args.append('-o')
    args.append(os.path.join(indir, meta['tpr-file']))

    try:
        run_gmx(' '.join(args))
    except GMXError as error:
        raise error
        # _, stderr = error.args
        # if 'Error' in error.stderr:
        #    raise error



def submit(meta):
    pass
