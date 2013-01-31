#!/bin/sh

rm out.size-energy
a0=29.9

for i in 0.9 0.92 0.94 0.96 0.98 1.0 1.02 1.04 1.06 1.08 1.1
do
a=`echo "scale=4; $a0*$i" | bc`
echo "a= $a"

cat >pmd000-000 <<EOF
250
${a}E+00  0.00E+00  0.00E+00
0.00E+00  ${a}E+00  0.00E+00
0.00E+00  0.00E+00  ${a}E+00
0.00E+00  0.00E+00  0.00E+00
0.00E+00  0.00E+00  0.00E+00
0.00E+00  0.00E+00  0.00E+00
EOF

cat pos >> pmd000-000
./pmd > out.pmd
erg=`grep 'Potential energy' out.pmd | awk '{print $3}'`
echo $a $erg >> out.size-energy
done

