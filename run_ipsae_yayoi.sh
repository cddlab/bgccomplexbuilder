#!/bin/bash
#PBS -q default
#PBS -l nodes=1:ppn=1:cpu
#PBS -l walltime=12:00:00
#PBS -m ae
#PBS -M moriwaki@bilab.sakura.ne.jp

test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
# run the environment module
. /home/apps/Modules/init/profile.sh

. .venv/bin/activate
test $START || { echo "START is not set." >&2 ; exit 1 ; }
test $END || { echo "END is not set." >&2 ; exit 1 ; }

.venv/bin/python3.12 complexbuilder/analysis/run_ipsae_yayoi.py -s $START -e $END
