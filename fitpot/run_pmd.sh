#!/bin/bash
#

pmdsrc="${HOME}/src/pmd"
expand_pmd="${pmdsrc}/util/expand_pmd.rb"
pmd="${pmdsrc}/src/pmd"
reduce="${pmdsrc}/fitpot/reduce_erg_frc.py"

python='/usr/local/bin/python'
mpirun='/opt/intel/impi/4.0.0.028/intel64/bin/mpirun'

usage=" Usage: $ ./run_pmd.sh 00??? [in.params.????]"

if [ $# -ne 2 ]; then
    echo "[Error] Number of arguments was wrong."
    echo $usage
    exit 1
fi

dir=$1
fparam=$2

mkdir -p $dir/pmd/0000
sed "s/cutoff_radius.*/cutoff_radius   $(head -n1 $fparam | awk '{print $2}')/" in.pmd > $dir/pmd/in.pmd
cp $fparam $dir/pmd/
cd $dir/pmd
$expand_pmd ../pos $(head -n1 $fparam | awk '{print $2}') > 0000/pmd00000
#.....run pmd
echo -n "."
#echo " running pmd on $dir/pmd ..."
$pmd > out.pmd
cd ../
$python $reduce
cd ../
