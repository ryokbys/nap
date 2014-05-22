#!/bin/sh
#-----------------------------------------------------------------------
# This script works only at helios supercomputer.
# To use sbatch system, we need to specify which node a qmcl uses.
#-----------------------------------------------------------------------

#.....Num of tasks for each vasp calculation
ntask=32
mfile="hostfile"

cwd=$(cd $(dirname $0); pwd)
qmcldir=$cwd
echo "run_force.sh: cwd= $cwd"
echo "run_force.sh: qmcldir= $qmcldir"

if [ ! -e $qmcldir/qmcl.local ]; then
    echo " Could not find qmcl.local executable file !!!"
    exit 1
fi

if [ $# -ne 2 ]; then
    name=`basename $0`
    echo " USAGE: $name <# ini> <# fin>" 1>&2
    exit 1
fi

nini=$1
nfin=$2

for islc in `seq $nini $nfin`
do
  cnum=`printf "%03d" $islc`

  if [ ! -e ./$cnum ]; then
      cp -r qmcldir $cnum
      #cd $cnum
      #$qmcldir/selcl
      #cd $cwd
  fi
  #.....make hostfile for each slice which will be used by vasp
  off0=`expr \( $islc - 1 \) \* $ntask`
  off1=`expr $islc \* $ntask`
  head -n $off1 $mfile | tail -n `expr $off1 - $off0` > $cnum/$mfile
  cp neb$cnum $cnum/qmcl000
  #cp in.qmcl$cnum $cnum/in.qmcl
  cd $cnum
  $qmcldir/qmcl.local > out.qmcl &
  pidarr[islc]=$!
  cd $cwd
done

for pid in ${pidarr[@]}
do
  wait $pid
done

for islc in `seq $nini $nfin`
do
  cnum=`printf "%03d" $islc`
  cp $cnum/frc000 ./frc$cnum
  cp $cnum/erg000 ./erg$cnum
done

