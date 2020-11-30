#!/bin/bash

mpirun -np 1 ../../fitpot/fitpot 2>&1 > out.fitpot
grep -e ENERGY -e FORCE -e STRESS out.fitpot.REF | tail -n3 | awk '{printf "%s %12.6f %12.6f %12.6f \n",$1,$4,$6,$8}' > /tmp/out.REF
grep -e ENERGY -e FORCE -e STRESS out.fitpot | tail -n3 | awk '{printf "%s %12.6f %12.6f %12.6f \n",$1,$4,$6,$8}' > /tmp/out
diff -q /tmp/out /tmp/out.REF
