#!/bin/bash
#PBS -N fitpot
#PBS -o out
#PBS -q batch
#PBS -j oe
#PBS -l nodes=12:ppn=6
#-----------------------------------------------------------------------
# Usage:
#   $ qsub run_vasp.sh
#-----------------------------------------------------------------------

export LANG=en_US
cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
cat $PBS_NODEFILE
echo 'NPROCS=' $NPROCS

for dir in learning_set*
do
  cp $PBS_NODEFILE $dir/nodelist.txt
done

fitpot=$HOME/src/nap/fitpot/fitpot.py
python=/usr/local/bin/python
jobid=$PBS_JOBID

$python $fitpot > out.fitpot 2>&1

#touch ${jobid}.done
