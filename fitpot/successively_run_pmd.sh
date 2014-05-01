#!/bin/bash
#

pmdsrc="${HOME}/src/pmd"
expand_pmd="${pmdsrc}/util/expand_pmd.rb"
pmd="${pmdsrc}/src/pmd"
reduce="${pmdsrc}/fitpot/reduce_erg_frc.py"

for dir in 0????
do
    mkdir -p $dir/pmd/0000
    sed "s/cutoff_radius.*/cutoff_radius   $(head -n1 in.params.SW_Si | awk '{print $2}')/" in.pmd > $dir/pmd/in.pmd
    cp in.params.SW_Si $dir/pmd/
    cd $dir/pmd
    $expand_pmd ../pos $(head -n1 in.params.SW_Si | awk '{print $2}') > 0000/pmd00000
    #.....run pmd
    echo -n "."
    #echo " running pmd on $dir/pmd ..."
    $pmd > out.pmd
    cd ../
    python $reduce
    cd ../
done
echo " running pmd done."
