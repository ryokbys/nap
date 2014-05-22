#!/bin/bash
#
# USAGE: $ ./geterg.vasp.sh 
#

ev2hrt=0.0367495737

#tail -n1 OSZICAR | sed "s/.*F=\(\s.*\)E0.*/\1/" > erg.vasp
tail -n1 OSZICAR | sed "s/.*F=\s\(.*\)\s.*/\1/" | awk '{printf "%22.14e\n", $1}' > erg.vasp

#awk "
#BEGIN {num=0}
#/F=/ { if (num==0) {print \$3*$ev2hrt}; num+=1}
#" < OSZICAR > erg.vasp
