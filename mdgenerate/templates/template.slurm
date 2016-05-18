#!/bin/bash -l

### Generated with mdgenerate.py, {{ timestamp }} ###

#SBATCH --workdir={{ root }}
#SBATCH --job-name={{ shortname }}
{%  for key, val in sbatch.items() %}
#SBATCH --{{ key }}={{ val }}
{% endfor %}

OUTDIR="{{ root }}/{{ outdir }}"
NAME="{{ shortname }}"
module purge
module load gromacs/5.1

gmx mdrun \
  -s {{ tpr }} \
  -o $OUTDIR/$NAME.trr \
  -c $OUTDIR/$NAME.gro \
  -x $OUTDIR/traj_full_$NAME.xtc \
  -e $OUTDIR/energy_$NAME.edr \
  -g $OUTDIR/log_$NAME.log \
  -cpo $OUTDIR/$NAME.cpt \
  -cpi $OUTDIR/$NAME.cpt \
  {% for key, val in mdrun.items() %} -{{ key }} {{ val }}{% endfor %}
