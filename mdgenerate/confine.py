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


class ResizedBoxFrame(md.coordinates.CoordinateFrame):
    _box = None
    _residue_ids = None

    @property
    def box(self):
        if self._box is None:
            return super().box
        else:
            return self._box

    @box.setter
    def box(self, val):
        if val.shape == (3,):
            self._box = np.eye(3) * val
        else:
            self._box = val

    @property
    def residue_ids(self):
        if self._residue_ids is None:
            return super().residue_ids
        else:
            return self._residue_ids

    @residue_ids.setter
    def residue_ids(self, val):
        self._residue_ids = val


def resize_to_box(frame, box, len_mol):
    """
    Cut a smaller box out of the frame.
    """
    if box.shape == (3, 3):
        box = box.diagonal()
    nr_res = len(frame) // len_mol

    residues = md.pbc.whole(frame % frame.box.diagonal(), len_res=len_mol).reshape(nr_res, len_mol, 3)
    masses = frame.masses.reshape(nr_res, len_mol, 1)
    masses /= masses.sum() / nr_res
    com = (residues * masses).sum(axis=1)
    res_in_box = np.where((com <= box).all(axis=1))
    new_frame = residues[res_in_box].view(ResizedBoxFrame).reshape(-1, 3)
    new_frame.box = box
    new_frame.residue_ids = np.zeros(len(new_frame), dtype=int)
    for i in range(len_mol):
        new_frame.residue_ids[i::len_mol] = res_in_box[0] + 1
    return new_frame


def make_spherical_conf(trajectory, constrained_subset, step, outfile, radius,
                        resize_box=False, method='residue', **kwargs):
    """
    Generate an initial configuration of spherically pinned molecules.

    Args:
        trajectory (mdevaluate.Coordinates):
            Bulk simulation from which the configuration is taken
        constrained_subset (dict):
            Definition of a subset of the atoms defining the constraints
        step: Timestep at which the configuration is taken from the source trajectory
        outfile: Output file of the new configuration
        radius: Radius of the sphercial confinement
        resize_box:
            If the siulation box should be resized accoriding to the size of confinement.
            When this is True, the box size is set to L = 2*radius + 2.0.
        method: Method by which molecules are constraint, possible values are:
            - residue: Change the residue of constraint molecules to the name given by
                       the keyword 'constrained_residue'
            - posres: Use position_restraints to constrain molecules
            - freeze: Use freeze groups to constrain molecules
    Returns:
        Tuple of number of unconstrained and constrained molecules.
    """
    subset = trajectory.subset(**constrained_subset)
    coords = subset[step]

    # Check if it is a water system with only one residue
    residue_ids = trajectory.atoms.residue_ids
    atom_names = trajectory.atoms.atom_names
    len_mol = len(set(atom_names))
    nr_mol = len(atom_names) // len_mol
    if len(set(residue_ids)) == 1:
        residue_ids[0::len_mol] = residue_ids[1::len_mol] = residue_ids[2::len_mol] = range(1, nr_mol + 1)
        coords = coords.view(ResizedBoxFrame)
        coords.residue_ids = residue_ids[0::len_mol]

    if resize_box:
        L = min(2*radius + 3.0, coords.box.max())
        coords = resize_to_box(coords, np.array([L, L, L]), len(coords)//nr_mol)

    coords %= coords.box.diagonal()
    coords -= coords.box.diagonal()/2
    r = md.coordinates.spherical_radius(coords, origin=0)

    constr_res = np.where(r >= radius, coords.residue_ids, -np.ones(r.shape, dtype=int))
    unconstrained = sum(constr_res == -1)
    constrained = len(constr_res) - unconstrained

    if method == 'residue':
        atoms = []
        if resize_box:
            frame = trajectory[step]
            atom_crds = md.pbc.whole(frame % frame.box.diagonal(), len_res=len_mol)
        else:
            atom_crds = trajectory[step] % trajectory[step].box.diagonal()

        for atm, res, resnr, x in zip(atom_names,
                                      trajectory.atoms.residue_names,
                                      residue_ids,
                                      atom_crds):
            atm_dict = {'atm': atm, 'resnr': resnr, 'x': x}
            if resnr not in coords.residue_ids:
                continue
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
        #TODO: write correct box in gro file!!! (do energy minimazation w/o posres first)
        write_gro(outfile, atoms, kwargs['name'], coords.box.diagonal()+kwargs.get('box_buffer', 0))
    else:
        raise NotImplementedError('method={} not implemented at the moment.'.format(method))

    return unconstrained, constrained


def make_slit_conf(trajectory, constrained_subset, step, outfile, thickness, **kwargs):
    """
    Generate an initial configuration of a slit pore of pinned molecules.

    Args:
        trajectory (mdevaluate.Coordinates):
            Bulk simulation from which the configuration is taken
        constrained_subset (dict):
            Definition of a subset of the atoms defining the constraints
        step: Timestep at which the configuration is taken from the source trajectory
        outfile: Output file of the new configuration
        radius: Radius of the sphercial confinement
        resize_box:
            If the siulation box should be resized accoriding to the size of confinement.
            When this is True, the box size is set to L = 2*radius + 2.0.
        method: Method by which molecules are constraint, possible values are:
            - residue: Change the residue of constraint molecules to the name given by
                       the keyword 'constrained_residue'
            - posres: Use position_restraints to constrain molecules
            - freeze: Use freeze groups to constrain molecules
    Returns:
        Tuple of number of unconstrained and constrained molecules.
    """
    subset = trajectory.subset(**constrained_subset)
    coords = subset[step]

    # Check if it is a water system with only one residue
    residue_ids = trajectory.atoms.residue_ids
    atom_names = trajectory.atoms.atom_names
    len_mol = len(set(atom_names))
    nr_mol = len(atom_names) // len_mol
    if len(set(residue_ids)) == 1:
        residue_ids[0::len_mol] = residue_ids[1::len_mol] = residue_ids[2::len_mol] = range(1, nr_mol + 1)
        coords = coords.view(ResizedBoxFrame)
        coords.residue_ids = residue_ids[0::len_mol]

    coords %= coords.box.diagonal()
    z = coords[:, 2]

    constr_res = np.where(z < thickness, coords.residue_ids, -np.ones(z.shape, dtype=int))
    unconstrained = sum(constr_res == -1)
    constrained = len(constr_res) - unconstrained

    atoms = []
    atom_crds = (trajectory[step] % trajectory[step].box.diagonal()).whole

    for atm, res, resnr, x in zip(atom_names,
                                  trajectory.atoms.residue_names,
                                  residue_ids,
                                  atom_crds):
        atm_dict = {'atm': atm, 'resnr': resnr, 'x': x}
        if resnr not in coords.residue_ids:
            continue
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
    write_gro(outfile, atoms, kwargs['name'], coords.box.diagonal())

    return unconstrained, constrained


def make_cylindrical_conf(trajectory, constrained_subset, step, outfile, radius, **kwargs):
    """
    Generate an initial configuration of a cylindrical pore of pinned molecules.

    Args:
        trajectory (mdevaluate.Coordinates):
            Bulk simulation from which the configuration is taken
        constrained_subset (dict):
            Definition of a subset of the atoms defining the constraints
        step: Timestep at which the configuration is taken from the source trajectory
        outfile: Output file of the new configuration
        radius: Radius of the sphercial confinement
        resize_box:
            If the siulation box should be resized accoriding to the size of confinement.
            When this is True, the box size is set to L = 2*radius + 2.0.
        method: Method by which molecules are constraint, possible values are:
            - residue: Change the residue of constraint molecules to the name given by
                       the keyword 'constrained_residue'
            - posres: Use position_restraints to constrain molecules
            - freeze: Use freeze groups to constrain molecules
    Returns:
        Tuple of number of unconstrained and constrained molecules.
    """
    subset = trajectory.subset(**constrained_subset)
    coords = subset[step]

    # Check if it is a water system with only one residue
    residue_ids = trajectory.atoms.residue_ids
    atom_names = trajectory.atoms.atom_names
    len_mol = len(set(atom_names))
    nr_mol = len(atom_names) // len_mol
    if len(set(residue_ids)) == 1:
        residue_ids[0::len_mol] = residue_ids[1::len_mol] = residue_ids[2::len_mol] = range(1, nr_mol + 1)
        coords = coords.view(ResizedBoxFrame)
        coords.residue_ids = residue_ids[0::len_mol]

    coords %= coords.box.diagonal()
    coords -= coords.box.diagonal()/2
    r, _ = md.coordinates.polar_coordinates(coords[:, 0], coords[:, 1])

    constr_res = np.where(r > radius, coords.residue_ids, -np.ones(r.shape, dtype=int))
    unconstrained = sum(constr_res == -1)
    constrained = len(constr_res) - unconstrained

    atoms = []
    atom_crds = (trajectory[step] % trajectory[step].box.diagonal()).whole

    for atm, res, resnr, x in zip(atom_names,
                                  trajectory.atoms.residue_names,
                                  residue_ids,
                                  atom_crds):
        atm_dict = {'atm': atm, 'resnr': resnr, 'x': x}
        if resnr not in coords.residue_ids:
            continue
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
    write_gro(outfile, atoms, kwargs['name'], coords.box.diagonal())

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


def generate_spherical_water(outdir, trajectory, step, radius, resize_box=False, **kwargs):
    """
    Generate gromacs topology for water in spherical neutral confinement.

    Args:
        outdir: Output directory of the new topology
        trajectory: Trajectory from which the starting configuration is taken
        step: Step at which the configuration is taken
        radius: Radius (in nm) of the spherical confinement
        resize_box (opt.): Resize the simulation bix to reduce the number of wall atoms
        **kwargs:aditional keyword arguments are passed to `make_spherical_conf`
    """
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass
    name = 'Water in spherical confinement, r={:.3}nm'.format(radius)
    nr_sol, nr_wal = make_spherical_conf(
        trajectory,
        {'atom_name': 'OW'},
        step,
        os.path.join(outdir, 'water.gro'),
        radius,
        constrained_residue='WAL',
        name=name,
        resize_box=resize_box,
        **kwargs
    )
    write_water_top(
        os.path.join(outdir, 'water.top'),
        name=name,
        nr_sol=nr_sol, nr_wal=nr_wal
    )


def generate_slit_water(outdir, trajectory, step, thickness, **kwargs):
    """
    Generate a neutral slit of water.
    """
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass
    name = 'Water in slit confinement, L={:.3}nm'.format(thickness)

    nr_sol, nr_wal = make_slit_conf(
        trajectory,
        {'atom_name': 'OW'},
        step,
        os.path.join(outdir, 'water.gro'),
        thickness,
        constrained_residue='WAL',
        name=name,
        **kwargs
    )
    write_water_top(
        os.path.join(outdir, 'water.top'),
        name=name,
        nr_sol=nr_sol, nr_wal=nr_wal
    )


def generate_cylindrical_water(outdir, trajectory, step, radius, **kwargs):
    """
    Generate a neutral slit of water.
    """
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass
    name = 'Water in cylindrical confinement, R={:.3}nm'.format(radius)

    nr_sol, nr_wal = make_cylindrical_conf(
        trajectory,
        {'atom_name': 'OW'},
        step,
        os.path.join(outdir, 'water.gro'),
        radius,
        constrained_residue='WAL',
        name=name,
        **kwargs
    )
    write_water_top(
        os.path.join(outdir, 'water.top'),
        name=name,
        nr_sol=nr_sol, nr_wal=nr_wal
    )
