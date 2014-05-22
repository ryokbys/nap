#!/bin/bash
#-----------------------------------------------------------------------
# run_QM.sh
# This script runs vasp as a QM calculator.
#-----------------------------------------------------------------------

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
MPIRUN=/usr/local/openmpi-1.2.8-intel64-v11.0.081/bin/mpirun
jobid=$PBS_JOBID

cwd=$(cd $(dirname $0); pwd)
vaspdir=${cwd}/vasp
if [ ! -d $vaspdir ]; then
    echo " There is no vasp directory !!!"
    exit 1
fi

cp QM_input $vaspdir/POSCAR
cd $vaspdir
vaspid=`qsub run_vasp.sh`
cd $cwd

#.....Wait until both programs end
until [ -f ${vaspdir}/${vaspid}.done ]
do
  sleep 1
done
cd $vaspdir
$qmcldir/getforce.vasp.sh
$qmcldir/geterg.vasp.sh
cd $cwd

#.....Copy force files to cwd
cp $vaspdir/frc.vasp ./frc.QM
cp $vaspdir/erg.vasp ./erg.QM
cp $vaspdir/OUTCAR ./OUTCAR.`printf "%03d" $istp`
