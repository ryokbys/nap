#!/bin/bash
#PJM -N fitpot-LLZO
#PJM -L "rscgrp=fx-small"
#PJM -L "node=8"
#PJM --mpi "proc=40"
#PJM --mpi "rank-map-bychip"
#PJM -L "elapse=1:00:00"
#PJM -X
#------------------------------------------------------
# Usage:
#   $ pjsub run_fitpot.fx.sh
#------------------------------------------------------
#

pmdsrc="/home/usr5/z48425a/src/fx/nap"
fitpot="${pmdsrc}/fitpot/fitpot"

echo "job started at `date`"
mpiexec --std out.fitpot $fitpot
echo "job finished at `date`"
