#!/bin/bash

script=$HOME/src/nap/nappy/vasp/make_random_deform_POSCARs.py

echo "making random deform POSCARs..."
python $script -n 1000 --deform-range=0.01,0.2 --displace-range=0.1 --num-displace=2 POSCAR

echo "making directories and moving POSCARs..."
for f in POSCAR_*
do
    dname=smpl_$(basename `pwd`)_random_`echo $f | sed 's/POSCAR_//'`
    mkdir -p $dname
    mv $f $dname/POSCAR
done

