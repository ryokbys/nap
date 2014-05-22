#!/bin/sh

cwd=$(cd $(dirname $0); pwd)
pmddir=${cwd}/../pmd/src


if [ ! -e $pmddir/pmd ]; then
    echo " Could not find pmd executable file !!!"
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
      mkdir $cnum
  fi
  cp neb$cnum $cnum/pmd00000-0000
  cp in.pmd $cnum/
  cd $cnum
  $pmddir/pmd > out.pmd &
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
  cp $cnum/frc0000 ./frc$cnum
  cp $cnum/erg0000 ./erg$cnum
done

