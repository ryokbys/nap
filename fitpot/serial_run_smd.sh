#!/bin/bash
#

usage=' Usage: ./serial_run_smd.sh in.params.SW_Si'

if [ $# -ne 1 ]; then
    echo "[Error] Number of arguments was wrong."
    echo $usage
    exit 1
fi

for dir in 0????
do
  ./run_smd.sh $1 $dir 
done
echo " running smd done."
