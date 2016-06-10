"""
This module provides functionality to generate simulations of confined systems.
"""
import os
import numpy as np
import mdevaluate as md

from .utils import save_open, write_gro
from .mdgenerate import env

POSRE_LINE = '{ind:>6d}  {func}  {params[0]:e} {params[1]:e} {params[2]:e}\n'


def write_posre(indices, itp_file='psore.itp', func_type=1, params=(1e6, 1e6, 1e6)):
    """
    Write indices of constrained atoms to an itp file.

    Args:
        indices: List of constraint atom indices, in gromacs format (starting at 1).
        itp_file (opt.): Filename of the itp file.
        func_type (opt.): Function type of the constraint.
        params (opt.): Parameters of the position restraint potential.
    """

    with save_open(itp_file) as itp:
        itp.write('[ position_restraints ]\n')
        fdict = {'func': func_type, 'params': params}
        for ind in indices:
            fdict['ind'] = ind
            itp.write(POSRE_LINE.format(**fdict))


def make_spherical_conf(trajectory, constraint_subset, step, outfile, radius,
                        method='residue', **kwargs):
    """
    Generate an initial configuration of spherically pinned molecules.

    Args:
        trajectory (mdevaluate.Coordinates):
            Bulk simulation from which the configuration is taken
        constraint_subset (dict):
            Definition of a subset of the atoms defining the constraints
        step: Timestep at which the configuration is taken from the source trajectory
        outfile: Output file of the new configuration
        radius: Radius of the sphercial confinement
        method: Method by which molecules are constraint, possible values are:
            - residue: Change the residue of constraint molecules to the name given by
                       the keyword 'constrained_residue'
            - posres: Use position_restraints to constrain molecules
            - freeze: Use freeze groups to constrain molecules
    Returns:
        Tuple of number of unconstrained and constrained molecules.
    """
    subset = trajectory.subset(**constraint_subset)
    coords = subset[step]
    coords %= coords.box.diagonal()
    coords -= coords.box.diagonal()/2
    r, _, _ = md.coordinates.spherical_coordinates(coords[:, 0], coords[:, 1], coords[:, 2])

    # Check if it is a water system with only one residue
    residue_ids = trajectory.atoms.residue_ids
    atom_names = trajectory.atoms.atom_names
    if len(set(residue_ids)) == 1:
        len_mol = len(set(atom_names))
        nr_mol = len(atom_names) // len_mol
        residue_ids[0::len_mol] = residue_ids[1::len_mol] = residue_ids[2::len_mol] = range(1, nr_mol + 1)

    constr_res = np.where(r >= radius, subset.atom_subset.residue_ids, -np.ones(r.shape, dtype=int))
    unconstrained = sum(constr_res == -1)
    constrained = len(constr_res) - unconstrained

    if method == 'residue':
        atoms = []
        atom_crds = trajectory[step] % trajectory[step].box.diagonal()
        for atm, res, resnr, x in zip(atom_names,
                                      trajectory.atoms.residue_names,
                                      residue_ids,
                                      atom_crds):
            atm_dict = {'atm': atm, 'resnr': resnr, 'x': x}
            if resnr in constr_res:
                atm_dict['res'] = kwargs['constrained_residue']
            else:
                atm_dict['res'] = res
            atoms.append(atm_dict)
        atoms = sorted(atoms, key=lambda x: x['res'])
        cur_resnr = -1
        i = 0
        for atm in atoms:
            if atm['resnr'] != cur_resnr:
                cur_resnr = atm['resnr']
                i += 1
            atm['resnr'] = i
        write_gro(outfile, atoms, kwargs['name'], trajectory[step].box.diagonal())
    else:
        raise NotImplementedError('method={} not implemented at the moment.'.format(method))

    return unconstrained, constrained


def write_water_top(topout, **kwargs):
    """
    Write a water topology with constrained molecules.

    Args:
        topout: Output file of the topology
        name: Name of the system
        nr_sol: Number of unconstrained molecules (SOL)
        nr_wal: Number of constrained molecules (WAL)
    """
    with save_open(topout) as top:
        top.write(
            env.get_template('water.top').render(**kwargs)
        )


def generate_spherical_water(outdir, trajectory, step, radius):
    """
    Generate gromacs topology for water in spherical neutral confinement.

    Args:
        outdir: Output directory of the new topology
        trajectory: Trajectory from which the starting configuration is taken
        step: Step at which the configuration is taken
        radius: Radius (in nm) of the spherical confinement
    """
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass
    name = 'Water in spherical confinement, r={}nm'.format(radius)
    nr_sol, nr_wal = make_spherical_conf(
        trajectory,
        {'atom_name': 'OW'},
        step,
        os.path.join(outdir, 'water.gro'),
        radius,
        constrained_residue='WAL',
        name=name
    )
    write_water_top(
        os.path.join(outdir, 'water.top'),
        name=name,
        nr_sol=nr_sol, nr_wal=nr_wal
    )
