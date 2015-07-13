#!/bin/bash
#PJM -N smd
#PJM -L "rscgrp=fx-small"
#PJM -L "node=8"
#PJM --mpi "proc=8"
#PJM --mpi "rank-map-bynode"
#PJM -L "elapse=10:00:00"
#PJM -X
#------------------------------------------------------
# Usage:
#   $ pjsub run.cx.sh
#------------------------------------------------------
#

pmdsrc="/home/usr5/z48425a/src/fx/nap"
smd="${pmdsrc}/pmd/smd"

NNODE=8
NPROC=8

fparam=in.params.NN

#.....for helios, uncomment following line
# module load intel intelmpi

for i in `seq 0 $NPROC`
do
    echo "(`expr $i % $NNODE`)" > rank$i
#    echo "($i)" > rank$i
done

echo "job started at `date`"

n=0
for dir in 0????
do
  # #...1st argument is not directory
  # if [ $dir = $1 ]; then
  #     continue
  # fi
  rest=`expr $n % $NPROC`
  mkdir -p $dir/smd
  sed "s/cutoff_radius.*/cutoff_radius   $(head -n1 $fparam | awk '{print $2}')/" in.smd > $dir/smd/in.smd
  cp $fparam $dir/smd/
  cd $dir/smd
  cp ../pos ./smd0000
  #.....run pmd
  echo -n "."
  mpiexec --vcoordfile ../../rank$rest -n 1 --std out.smd $smd &
  cd ../
  cp smd/erg.smd ./
  cp smd/frc.smd ./
  cd ../
  n=`expr $n + 1`
  if test `expr $n % $NPROC` -eq 0; then
      wait
      n=0
  fi
done


echo "job finished at `date`"
