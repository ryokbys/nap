#!/bin/sh

rm out.size-energy
aa2bohr=1.88972616356
#.....set range
a0=10.0
a1=20.0

#for i in 0.5 0.6 0.7 0.8 0.9 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7
for i in `seq 0 99`
do
r=`echo "scale=4; 1.0/100*($a1-$a0)*$i" | bc` 
a=`echo "scale=4; ($a0+$r)*$aa2bohr" | bc`
echo "a= $a"

cat >pmd00000-0000 <<EOF
2
${a}E+00  0.00E+00  0.00E+00
0.00E+00  38.0E+00  0.00E+00
0.00E+00  0.00E+00  25.0E+00
0.00E+00  0.00E+00  0.00E+00
0.00E+00  0.00E+00  0.00E+00
0.00E+00  0.00E+00  0.00E+00
EOF

cat pos >> pmd00000-0000
./pmd > out.pmd
erg=`grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'`
echo $a $erg >> out.size-energy
done

