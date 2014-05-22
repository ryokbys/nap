#!/bin/bash
#-----------------------------------------------------------------------
# run_CL+QM.sh
# This script runs pmd and vasp in the directories of the same names
# using CL_input and QM_input, respectively.
#-----------------------------------------------------------------------
# USAGE:
#   $ ./run_pmd+vasp.sh 10
#     * $1 = step-id
#

cwd=$(cd $(dirname $0); pwd)
pmddir=${cwd}/pmd
vaspdir=${cwd}/vasp
qmcldir=${cwd}/src

istp= $1

#.....Error check
if [ ! -e $pmddir/pmd ]; then
    echo " Could not find pmd executable file !!!"
    exit 1
fi

if [ ! -d $vaspdir ]; then
    echo " There is no vasp directory !!!"
    exit 1
fi

#.....Run vasp first, because it takes a lot of time
cp QM_input $vaspdir/POSCAR
cd $vaspdir
vaspid=`qsub run_vasp.sh`
cd $cwd

#.....Run pmd next
mkdir -p $pmddir/0000
cp CL_input $pmddir/0000/
cd $pmddir
./pmd > out.pmd
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
cp $pmddir/frc0000 ./frc.CL
cp $pmddir/erg0000 ./erg.CL
cp $vaspdir/frc.vasp ./frc.QM
cp $vaspdir/erg.vasp ./erg.QM
cp $vaspdir/OUTCAR ./OUTCAR.`printf "%03d" $istp`
