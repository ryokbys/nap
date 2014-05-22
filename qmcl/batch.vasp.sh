#!/bin/sh
#PBS -N vasp-Al-cluster
#PBS -o out
#PBS -q batch
#PBS -j oe
#PBS -l nodes=16:ppn=2
#-----------------------------------------------------------------------
# Usage:
#   $ qsub run_vasp.sh
#-----------------------------------------------------------------------

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
MPIRUN=/usr/local/openmpi-1.2.8-intel64-v11.0.081/bin/mpirun
jobid=$PBS_JOBID

rm OUTCAR
#rm CHG* CONTCAR WAVECAR
$MPIRUN -np $NPROCS vasp

touch ${jobid}.done
