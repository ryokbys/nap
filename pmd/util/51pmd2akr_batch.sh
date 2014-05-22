#!/bin/bash
#
# USAGE:
#   ./51pmd2akr_batch.sh #start #end #skip
#

execdir=`dirname $0`
#for arg in `seq $1 $3 $2` # for GNU
for arg in `jot 1000 $1 $2 $3` # for BSD
do
  ${execdir}/50pmd2akr.py ${arg} 7 8 9 10 11
done
