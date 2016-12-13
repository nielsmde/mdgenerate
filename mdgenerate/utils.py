import os
import shutil

import numpy as np
from mdevaluate import pbc


def save_open(filename, mode='w'):
    """
    Open a file savely by backing up an existing file with this name.
    """

    if os.path.exists(filename) and mode in ['w', 'wb']:
        i = 0
        backup_format = '{}/#{}.{}'.format(*os.path.split(filename), '{i}')
        while os.path.exists(backup_format.format(i=i)):
            i += 1
        shutil.copy(filename, backup_format.format(i=i))
        # print('Backing up file: {} -> {}'.format(filename, backup_format.format(i=i)))
    return open(filename, mode=mode)

GRO_HEAD = '{name}\n{nmol}\n'
GRO_BODY = '{resnr:>5}{res:<5}{atm:>5}{atmnr:>5}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{v[0]}{v[1]}{v[2]}\n'
GRO_BOX = '{box[0]:10.05} {box[1]:10.05} {box[2]:10.05}\n'


def write_gro(gro_file, atoms, name, box):
    """
    Write atoms to a gro file.

    Args:
        gro_file: Gro file to write
        atoms: List of dicts with the keys:
            - 'res': residue name
            - 'resnr': residue index
            - 'atm': atom name
            - 'x': coordinates
            - 'v': velocities (optional)
        name: Name of the configuration
        box: Box size vector of the configuration

    """
    with save_open(gro_file) as gro:
        gro.write(GRO_HEAD.format(name=name, nmol=len(atoms)))
        for i, atom in enumerate(atoms):
            atom.setdefault('v', ('', '', ''))
            gro.write(GRO_BODY.format(atmnr=i+1, **atom))
        gro.write(GRO_BOX.format(box=box))


def write_index(groups, ndx_file='index.ndx'):
    """
    Append index groups to an index file.

    Args:
        groups: Dict of index groups.
        ndx_file (opt.): Index file that will be written or appended if it exists.
    """

    with open(ndx_file, 'a') as ndx:
        for grp_name, grp_indices in groups.items():
            ndx.write()
