#!/bin/bash
#-----------------------------------------------------------------------
# run_CL.sh
# This script runs pmd as a CL calculator.
#-----------------------------------------------------------------------
# USAGE:
#   $ ./run_pmd.sh 10
#     * $1 = step-id
#

cwd=$(cd $(dirname $0); pwd)
pmddir=$cwd/pmd
istp=$1

#.....Error check
if [ ! -e $pmddir/pmd ]; then
    echo " Could not find pmd executable file !!!"
    exit 1
fi

#.....Run pmd next
cp CL_input $pmddir/0000/pmd00000
cd $pmddir
./pmd > out.pmd
#./pmd.local > out.pmd
wait $!

cd $cwd

#.....Copy force files to cwd
cp $pmddir/frc0000 ./frc.CL
cp $pmddir/erg0000 ./erg.CL
