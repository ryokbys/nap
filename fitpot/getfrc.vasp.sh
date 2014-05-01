#!/bin/bash
#
# USAGE: $ ./getforce.vasp.sh 
#

tmpfile='tmp.frc'
outfile='frc.ref'

nfrc=`grep 'TOTAL-FORCE' OUTCAR | wc -l`
ev2hrt=0.0367495737
aa2bohr=1.88972616356
convf=0.019447

awk "
BEGIN { num=0 }
/total drift:/ { start=0 }
!/----------/  { if ( start==1 && num==$nfrc ) { printf \"%14.8f  %14.8f  %14.8f\n\",\$4,\$5,\$6  }}
/TOTAL-FORCE/  { start=1; num+=1 }
" < OUTCAR > $tmpfile

sed "1 i\\
$(wc -l < $tmpfile)
" $tmpfile > $outfile
rm $tmpfile
