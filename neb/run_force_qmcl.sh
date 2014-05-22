#!/bin/sh

cwd=$(cd $(dirname $0); pwd)
qmcldir=/home/kobayashi/src/qmcl


if [ ! -e $qmcldir/qmcl ]; then
    echo " Could not find qmcl executable file !!!"
    exit 1
fi

if [ $# -ne 2 ]; then
    name=`basename $0`
    echo " USAGE: $name <# start> <# end>" 1>&2
    exit 1
fi

nstart=$1
nend=$2

for islc in `seq $nstart $nend`
do
  cnum=`printf "%03d" $islc`

  if [ ! -e ./$cnum ]; then
      cp -r qmcldir $cnum
      cd $cnum
      $qmcldir/selcl
      cd $cwd
  fi
  cp neb$cnum $cnum/qmcl000
  cp in.qmcl$cnum $cnum/in.qmcl
  cd $cnum
  $qmcldir/qmcl > out.qmcl &
  pidarr[islc]=$!
  cd $cwd
done

for pid in ${pidarr[@]}
do
  wait $pid
done

for islc in `seq $nstart $nend`
do
  cnum=`printf "%03d" $islc`
  cp $cnum/frc000 ./frc$cnum
  cp $cnum/erg000 ./erg$cnum
done

