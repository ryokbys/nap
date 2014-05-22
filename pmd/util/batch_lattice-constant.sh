#!/bin/sh

PMDDIR=/home/kobayashi/2012/06-25_shear-mod-LJ2D/src/
for r in `seq 0.99 0.001 1.01`
do
  for file in orig*
  do
    $PMDDIR/boxsize $file $r `echo $file | sed 's/orig/pmd/'`
  done > tmp
  echo "`printf \"%8.3f\" $r`  `mpirun -n 1 $PMDDIR/pmd | sed -n 's/Potential energy=//p'`"
  #echo "`printf \"%8.3f\" $r`  `mpirun -n 1 $PMDDIR/pmd | sed -n 's/arho,phi =//p'`"
done

