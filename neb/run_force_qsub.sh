#!/bin/sh
#-----------------------------------------------------------------------
# Please use this script by renaming to run_force.sh
#-----------------------------------------------------------------------

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
  cp neb$cnum $cnum/pmd000-000
  cp in.pmd 30pmdrun.king $cnum/
  cd $cnum
  # $pmddir/pmd > out.pmd &
  pmdid=`qsub 30pmdrun.king`
  pidarr[islc]=$pmdid
  cd $cwd
done

#.....Wait until all the programs end
# for pid in ${pidarr[@]}
# do
#   wait $pid
# done

for islc in `seq $nstart $nend`
do
  cnum=`printf "%03d" $islc`
  #.....Wait until the program ends
  until [ -f $cnum/${pidarr[islc]}.done ]
  do
    sleep 0.1
  done
  cp $cnum/frc000 ./frc$cnum
  cat $cnum/erg000 >> ./erg$cnum
  rm $cnum/${pidarr[islc]}.done
done

