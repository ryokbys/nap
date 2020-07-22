#!/bin/bash

for d in dataset_Si/smpl_*
do
    echo $d
    mkdir -p $d/pmd
    cp $d/pos $d/pmd/pmdini
    cp in.params.desc $d/pmd/
done
