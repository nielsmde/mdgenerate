"""
Python module to generate GROMACS simulation input
and submit them to SLURM.
"""
from .utils import save_open as open, write_gro

import re
import os
import shutil
import time
import subprocess

import numpy as np

import mdevaluate as md
# from subprocess import run, PIPE, STDOUT

import yaml
import jinja2

GMX_VERSION = '5.1.3'

DEFAULT_KEYS = {
    'mdp_file': 'mdpin.mdp',
    'slurm_file': 'slurm.sh',
    'tpr_file': 'topol.tpr',
    'indir': '',
    'outdir': '',
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
        self['timestamp'] = time.strftime('%d.%m.%Y %H:%M')

    def _todict(self):
        d = dict()
        for key in self:
            d[key] = self[key]
        return d

    @property
    def as_mdp(self):
        self._timestamp()
        return env.get_template('template.mdp').render(**self._todict())

    @property
    def as_slurm(self):
        self._timestamp()
        return env.get_template('template.slurm').render(**self._todict())

    def load_yaml(yaml_file):
        """Load a yaml file and process its template."""
        with open(yaml_file) as f:
            yaml_dict = MetaDict(yaml.load(f.read()), yaml_file)

        if 'extends' in yaml_dict:
            template_dict = MetaDict(locate_template(yaml_dict['extends']))
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
        with open(yaml_file, 'r') as f:
            yaml_dict = yaml.load(f.read())

        if 'extends' in yaml_dict:
            if isinstance(yaml_dict['extends'], list):
                templates = yaml_dict['extends']
                template_dict = MetaDict(locate_template(templates[0]))
                for template in templates:
                    template_dict.update(MetaDict(locate_template(template)))
            else:
                template_dict = MetaDict(locate_template(yaml_dict['extends']))
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
        self['root'] = self.directory
        self['filename'] = self.filename

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], str):
            self._init_with_yaml(args[0])
        else:
            super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        if key == 'directory':
            return self.directory
        try:
            val = super().__getitem__(key)
            if 'file' in key or key in ['indir', 'outdir']:
                if isinstance(val, str):
                    val = os.path.join(self.directory, val)
                elif isinstance(val, list):
                    val = [os.path.join(self.directory, v) for v in val]
            return val
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
#    if dirs or files.remove(meta.filename):
#        if not override:
#            raise FileExistsError('Simulation directory is not empty: {}'.format(meta.directory))
#        else:
#            print('Simulation directory is not empty, but override=True')

    # indir = os.path.join(meta.directory, meta.get('indir', ''))
    # outdir = os.path.join(meta.directory, meta.get('outdir', ''))
    os.makedirs(meta['indir'], exist_ok=True)
    os.makedirs(meta['outdir'], exist_ok=True)

    if meta.get('copy-topology', False):
        for f in meta['topology-files']:
            shutil.copy(f, meta['indir'])

    with open(meta['mdp_file'], 'w') as f:
        f.write(meta.as_mdp)


def generate_slurm(meta):
    slurm_file = os.path.join(meta.directory, meta.get('slurm_file', 'slurm.sh'))
    with open(slurm_file, 'w') as f:
        f.write(meta.as_slurm)


def grompp(meta, generate=True):
    if generate:
        generate_files(meta)

    indir = os.path.join(meta['directory'], meta['indir'])
    args = [
        'gmx', 'grompp',
        '-f', meta['mdp_file'],
        '-po', os.path.join(indir, 'mdout.mdp')
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
    args.append(os.path.join(indir, meta['tpr_file']))

    for arg, val in meta.get('grompp', {}).items():
        args.append('-{}'.format(arg))
        args.append(str(val))

    try:
        run_gmx(' '.join(args))
    except GMXError as error:
        raise error
        # _, stderr = error.args
        # if 'Error' in error.stderr:
        #    raise error


def submit(meta):
    pass


def process(yaml_file, submit=False, run_grompp=True):
    """
    Process a YAML file and generate all files. If submit=True submit it to slurm.
    """
    meta = MetaDict(yaml_file)
    global GMX_VERSION
    if 'gmx_version' in meta:
        GMX_VERSION = meta['gmx_version']
    else:
        meta['gmx_version'] = GMX_VERSION
    if run_grompp:
        grompp(meta)
    generate_slurm(meta)
    if submit:
        submit(meta)


def nvt_from_npt(trajectory, topology, outfile, edrfile=None, window=0.3, accuracy=1e-5, use_best=False):
    """
    Find the mean boxsize of a NpT trajectory and write out the latest frame
    with this boxsize as gro file.

    Args:
        trajectory: Trajectory file of the NpT simulation
        topology: Topology of the NpT simulation
        outfile: Filename of the output gro file
        window (opt.): Fraction of the trajectory the box size is averaged over
        accuracy (opt.): Maximum deviation of the boxsize from the mean value
    """
    npt_traj = md.open('', trajectory=trajectory, topology=topology, cached=True, verbose=False)
    if edrfile is not None:
        energy = md.open_energy(edrfile)
        boxsize = np.eye(3) * [energy[dim][-int(window * len(energy.time)):].mean() for dim in ('Box-X', 'Box-Y', 'Box-Z')]
    else:
        boxsize = np.mean([f.box for f in npt_traj[-int(window * len(npt_traj)):]], axis=0)
    fbox = 0
    i = 0
    best = (np.inf, 0)
    difference = np.inf
    while difference > accuracy:
        i += 1
        fbox = npt_traj[-i].box
        difference = np.absolute(fbox - boxsize).max()
        if difference < best[0]:
            best = (difference, i)
        if i > window * len(npt_traj):
            if use_best:
                i = best[1]
                break
            else:
                raise ValueError("Adequate Box size was not found. Best: {}".format(best))

    time = npt_traj[-i].time
    command = """gmx trjconv -o {out} -s {top} -f {xtc} -pbc whole -b {t} -e {t} <<EOF
    0
    EOF
    """
    run_gmx(command.format(out=outfile, top=topology, xtc=trajectory, t=time))


def make_short_sims(directory, N, dt=0.01, timeshift=1000.0, queue='nodes'):
    """
    Create some short simulations for a given simulation.
    """
    traj = md.open(directory, trajectory='out/*.xtc')
    time_max = traj[-1].time
    timestep = round(traj[1].time - traj[0].time, 3)
    time = timestep * 11
    if time < 2000 and queue == 'nodes':
        queue = 'short'
    if timestep <= dt:
        print("Timestep of trajectory is {:.3f}, dt={:.3f}. Skipping...".format(timestep, dt))
        return
    timeshift = max(timeshift, 10 * timestep)

    if N * timeshift > time_max / 2:
        timeshift = time_max / (2 * N)

    print('timestep=', timestep, 'time=', time, 'timeshift=', timeshift)

    with open(os.path.join(directory, 'meta.yaml'), mode='r') as f:
        meta = yaml.load(f.read())

    meta['name'] += '_sh'
    meta['time'] = '{:.3f}ps'.format(time)
    meta['mdp']['nstxout-compressed'] = int(dt // time_to_ps(meta['timestep']))
    topfiles = meta['topology-files']
    meta['topology-files'] = []
    for top in topfiles:
        if '.top' in top:
            top = '../../' + top
        meta['topology-files'].append(top)
    meta['extends'] = [queue, 'nvt']

    for i in range(N):
        shift = - int(((i * timeshift) // timestep + 1))
        step = len(traj) + shift
        print('shifting:', step, shift, traj[shift].time)
        short_dir = os.path.join(directory, 'short', str(step))
        short_yaml = os.path.join(short_dir, 'meta.yaml')
        if os.path.exists(os.path.join(short_dir, 'topol.tpr')):
            continue

        os.makedirs(os.path.join(short_dir), exist_ok=True)
        with open(short_yaml, mode='w') as f:
            # f.write(yaml.dump(meta))
            yaml.dump(meta, stream=f, default_flow_style=False)

        frame = traj[step]
        atoms = [{'res': res, 'resnr': resnr, 'atm': atm, 'x': x} for (res, resnr, atm, x) in zip(
            frame.residue_names, frame.residue_ids, frame.atom_names, frame
        )]
        write_gro(
            os.path.join(short_dir, 'biwater.gro'),
            atoms, 'Configuration at step={}'.format(step), traj[step].box.diagonal()
        )
        process(short_yaml)
