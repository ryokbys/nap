#!/opt/local/bin/python

u"""
Make in.const.linreg to be used in linreg potential.
"""

import sys,os

constfname='in.const.linreg'

max_nexp= 2
#.....Gaussian
max_exp_gauss= 5.0
num_exp_gauss= 10
max_len_gauss= 10.0
num_len_gauss= 10
#.....cosine
max_q_cos= 10.0
num_q_cos= 100
#.....polynomial
a0_max= 1.0
a0_min=-1.0
a2_max= 1.0
a2_min=-1.0
a6_max= 1.0
a6_min=-1.0
a10_max= 1.0
a10_min=-1.0
num_a0= 3
num_a2= 5
num_a6= 5
num_a10= 5

f=open(constfname,'w')
nelem= num_exp_gauss *num_len_gauss +num_q_cos \
    +num_a0 *num_a2 *num_a6 *num_a10
#.....1st line
f.write(' {0:10d} {1:4d}\n'.format(nelem,max_nexp))
#.....Gaussian
dexp= max_exp_gauss/num_exp_gauss
dlen= max_len_gauss/num_len_gauss
for i1 in range(1,num_exp_gauss+1):
    for i2 in range(1,num_len_gauss+1):
        f.write(' {0:3d} {1:10.4f} {2:10.4f}\n'.format(1,i1*dexp,i2*dlen))
#.....cosine
dq= max_q_cos/num_q_cos
for i1 in range(1,num_q_cos+1):
    f.write(' {0:3d} {1:10.4f}\n'.format(2,i1*dq))
#.....polynomial
d0= (a0_max-a0_min)/(num_a0-1)
d2= (a2_max-a2_min)/(num_a2-1)
d6= (a6_max-a6_min)/(num_a6-1)
d10= (a10_max-a10_min)/(num_a10-1)
for i0 in range(num_a0):
    for i2 in range(num_a2):
        for i6 in range(num_a6):
            for i10 in range(num_a10):
                f.write(' {0:3d} {1:10.4f}'.format(3,i0*d0+a0_min))
                f.write(' {0:10.4f}'.format(i2*d2+a2_min))
                f.write(' {0:10.4f}'.format(i6*d6+a6_min))
                f.write(' {0:10.4f}'.format(i10*d10+a10_min))
                f.write('\n')
f.close()
