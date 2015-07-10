#!/bin/bash
#PJM -N fitpot-LLZO
#PJM -L "rscgrp=cx-debug"
#PJM -L "vnode=4"
#PJM -L "vnode-core=12"
#PJM -L "elapse=0:30:00"
#PJM -P "vn-policy=pack"
#PJM -X
#------------------------------------------------------
# Usage:
#   $ pjsub run_fitpot.fx.sh
#------------------------------------------------------
#

source /center/local/apl/cx/intel/composer_xe_2013_sp1/bin/compilervars.sh intel64
source /center/local/apl/cx/intel/impi/4.1.1.036/intel64/bin/mpivars.sh intel64
source /center/local/apl/cx/intel/mkl/bin/mklvars.sh intel64

pmdsrc="/home/usr5/z48425a/src/cx/nap"
fitpot="${pmdsrc}/fitpot/fitpot"

echo "job started at `date`"
mpiexec.hydra -n 48 -print-rank-map $fitpot > out.fitpot 2>&1
echo "job finished at `date`"
