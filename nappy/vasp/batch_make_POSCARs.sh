#!/bin/bash

script=$HOME/src/nap/nappy/vasp/make_deformed_POSCARs.py

python $script isotropic --fmax=1.5 --fmin=0.8 -n 20 --offset 0 POSCAR
python $script uniaxial --fmax=1.1 --fmin=0.9 -n 15 --offset 100 POSCAR
python $script orthorhombic --fmax=1.1 --fmin=0.9 -n 15 --offset 200 POSCAR
python $script shear --fmax=1.1 --fmin=0.9 -n 15 --offset 300 POSCAR

for f in POSCAR-*
do
    dname=smpl_$(basename `pwd`)_deformed_`echo $f | sed 's/POSCAR-//'`
    mkdir -p $dname
    mv $f $dname/POSCAR
done
