#!/bin/bash

script=$HOME/src/nap/nappy/vasp/make_deformed_POSCARs.py

echo "isotropic..."
python $script isotropic --fmax=1.5 --fmin=0.8 -n 50 --offset 0 POSCAR
echo "uniaxial..."
python $script uniaxial --fmax=1.2 --fmin=0.9 -n 50 --offset 100 POSCAR
echo "orthorhombic..."
python $script orthorhombic --fmax=1.2 --fmin=0.9 -n 50 --offset 200 POSCAR
echo "shear..."
python $script shear --fmax=1.2 --fmin=0.9 -n 50 --offset 300 POSCAR

echo "making directory and moving POSCARs..."
for f in POSCAR-*
do
    dname=smpl_$(basename `pwd`)_deform_`echo $f | sed 's/POSCAR-//'`
    mkdir -p $dname
    mv $f $dname/POSCAR
done

