#!/bin/sh
#-----------------------------------------------------------------------
# Please use this script by renaming to run_force.sh
#-----------------------------------------------------------------------

cwd=$(cd $(dirname $0); pwd)

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
  pmd2vasp < neb$cnum > $cnum/POSCAR
  cp INCAR POTCAR KPOINTS run_vasp.sh $cnum/
  cd $cnum
  pmdid=`qsub run_vasp.sh`
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
    sleep 1
  done
  cd $cnum
  $cwd/getforce.vasp.sh
  $cwd/geterg.vasp.sh
  cp frc.vasp $cwd/frc$cnum
  cp erg.vasp $cwd/erg$cnum
  cd $cwd
  rm $cnum/${pidarr[islc]}.done
done

