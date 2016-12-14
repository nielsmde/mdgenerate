#!/bin/bash -l

### Generated with mdgenerate.py, {{ timestamp }} ###

#SBATCH --workdir={{ root }}
#SBATCH --job-name={{ name }}
{% if sbatch is defined %}
{%  for key, val in sbatch.items() %}#SBATCH --{{ key }}={{ val }}
{% endfor %}
{% endif %}

OUTDIR="{{ outdir }}"
NAME="{{ name }}"
TPRFILE="{{ tpr_file }}"
XTCFILE="$OUTDIR/traj_full_$NAME.xtc"
module purge
module load gromacs/{{ gmx_version }}

gmx mdrun \
  -s $TPRFILE \
  -o $OUTDIR/$NAME.trr \
  -c $OUTDIR/$NAME.gro \
  -x $XTCFILE \
  -e $OUTDIR/energy_$NAME.edr \
  -g $OUTDIR/log_$NAME.log \
  -cpo $OUTDIR/$NAME.cpt \
  -cpi $OUTDIR/$NAME.cpt \
 {% for key, val in mdrun.items() %} -{{ key }} {{ val }}{% endfor %}

touch $OUTDIR/finished

{% if do_index is defined and do_index %}
module load mdevaluate/dev
index-xtc --nojump --tpr=$TPRFILE $XTCFILE
{% endif %}
