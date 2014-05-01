#!/bin/bash
#
# USAGE: $ ./geterg.vasp.sh 
#

outfile='erg.ref'

grep 'energy without' OUTCAR  | awk '{print $5}' | tail -n1 > $outfile
