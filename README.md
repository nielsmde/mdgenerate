# MDGenerate: GROMACS Simulation generation

This Python package provides functionality to generate input for Gromacs simulations
and submit them to a SLURM queue.
It is based on a YAML-templating system which defines the simulation parameters.
The templating system supports multiple inheritance which allows the user to define
several parameters of the simulation by using predefined templates.

## YAML-files

The parameters of the simulation are defined in YAML files.
For this package only a subset of the [YAML syntax](http://pyyaml.org/wiki/PyYAMLDocumentation#YAMLsyntax) is used:

1. Key value pairs are defined with a colon: `key: value`

2. Lists are defined with a dash: `- item`

A simulation should be defined through a YAML file, e.g. `meta.yaml` in an empty directory, which will be the root directory of the new simulation.
All simulation properties are defined through key values pairs,
where the value may be a literal (like string or float) or a nested dictionary of
key value pairs or a list.

### Templates

Each YAML file can have a key `extends`, which defines one or more templates.
Templates can be given by full path, or a name (or relative path) which will be
looked for in the standard template location `~/.mdgenerate/templates`.
If a list of templates is given, the inheritance is carried out from top to
bottom of the list.
All keys of the templates will be defined in the resulting simulation but can be overwritten in the YAML file of the simulation or any subsequent template.
Keys that define a dictionary (e.g. mdp, so below) will be merged when templates are applied.

### Example

This is an example how a YAML file can look like:
```
extends:   
  - base
  - nodes

name:         MySim
time:         3ns
timestep:     1fs
temperature:  300

topology-files:
  - /path/to/topology.top
  - /path/to/conf.gro
  - /path/to/index.ndx

indir:        in
outdir:       out

mdp:
  - tcouple:    Nose-Hoover
  - pcouple:    Berendsen

```

This YAML file is based on two templates, which could for example define options for mdrun or sbatch.

## Mandatory keys

All keys mentioned in this section have to be defined in the YAML file
(or a template file) before the simulation can be generated.

- **name**: The name identifying the simulation.
- **time**: The whole simulation time.
- **timestep**: The time step of the simulation.
- **topology-files**: A list of topology files, e.g. top, gro, ndx.

## Optional keys

The following keys are optional in the sense, that `mdgenerate` will work if they are not defined, but they may still be necessary for the simulation to run properly.
Default values are given in brackets, if they are not None.

- **temperature**
- **pressure**
- **indir**: Sub-Directory  where input files will be located.
- **outdir**: Sub-Directory where output files will be located.
- **mdp-file** [mdpin.mdp]: Name of the generated mdp file.
- **slurm-file** [slurm.sh]: Name of the generated slurm script.

## Parameter dictionaries

There are three sets of parameters which can be defined for the simulations:
**mdp**, **sbatch** and  **mdrun**.
As shown in the example, these dictionary are defined in the manner of nested key value pairs.


## Time values

For time values `mdgenerate` supports time units (fs, ps, ns, us, ms, s) which will be converted to the float value in pico seconds in the generation process.
