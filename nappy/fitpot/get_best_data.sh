#!/bin/bash

if [ -d best_data ]; then
    rm -rf best_data
fi
num=$(grep 'step,time' out.fp | tail -n1 | awk '{print $4}')
#num=$(ls in.vars.fitpot.[0-9]* | sed 's/in.vars.fitpot.//' | sort -n | tail -n 1)
echo "copying iid_${num} ==> ./best_data"
cp -r $(find . -name "iid_${num}") ./best_data
