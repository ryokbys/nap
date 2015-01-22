#!/bin/bash
#PJM -N fitpot-NN1        # job name
#PJM -L "rscgrp=cx-debug"
#PJM -L "vnode=40, vnode-core=1"
#PJM -P "vn-policy=abs-unpack=20"
# #PJM -P "vn-policy=abs-unpack=12"
#PJM -L "elapse=0:10:00"    # execution time
#PJM -j
#PJM -X
#PJM -S
#-----------------------------------------------------------------------
# Usage:
#   $ pjsub batch.fitpot.sh
#-----------------------------------------------------------------------
# #PJM -P "vn-policy=abs-unpack=10"

# source /center/local/apl/cx/intel/bin/compilervars.sh intel64
# source /center/local/apl/cx/intel/mkl/bin/mklvars.sh intel64
# source /center/local/apl/cx/intel/impi/4.1.1.036/intel64/bin/mpivars.sh

export LANG=en_US
cp $PJM_O_NODEINF ./pjm_o_nodeinf
NPROCS=`wc -l < ./pjm_o_nodeinf`
NNODES=`uniq ./pjm_o_nodeinf | wc -l`
echo 'NNODES=' $NNODES
echo 'NPROCS=' $NPROCS

for dir in learning_set*
do
    for i in `seq 0 $(wc -l < ./pjm_o_nodeinf)`
    do
	printf '(%d)' $i > $dir/rankfile$i
    done
done

# mpdboot -n $NNODES -f ./pjm_o_nodeinf -r /bin/pjrsh

fitpot=$HOME/src/nap/fitpot/fitpot.py
python=/usr/bin/python
jobid=$PBS_JOBID

$python $fitpot > out.fitpot 2>&1

#touch ${jobid}.done
# mpdallexit

#......remove all the temporally files 
# rm learning_set*/rankfile*
