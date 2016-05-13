#!/bin/bash

smpls=$*

for dir in $smpls
do
   echo $dir
   mkdir -p $dir/smd
   cp in.const.NN $dir/smd/
done

