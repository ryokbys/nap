#!/bin/bash

rm -rf best_data
num=$(ls in.vars.fitpot.[0-9]* | sed 's/in.vars.fitpot.//' | sort -n | tail -n 1)
cp -r $(find . -name "iid_${num}") ./best_data
