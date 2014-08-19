#!/opt/local/bin/python

u"""
Make in.const.linreg to be used in linreg potential.
"""

import sys,os

constfname='in.const.linreg'

nexp= 2 # which includes [0.5, 1]
exps=[1.0, 2.0, 0.5]
#.....Gaussian
range_exp_gauss=[0.0,5.0]
num_exp_gauss= 5
range_len_gauss=[0.0,10.0]
num_len_gauss= 5
#.....cosine
range_q_cos=[0.0,10.0]
range_exp_cos=[-2.0,0.0]
num_q_cos= 10
num_exp_cos= 3
#.....polynomial
range_a0_poly= [-1.0,1.0]
range_a2_poly= [-1.0,1.0]
range_a6_poly= [-1.0,1.0]
range_a10_poly= [-1.0,1.0]
num_a0= 3
num_a2= 3
num_a6= 3
num_a10= 3
#.....angular
range_angle=[-1.0,1.0]
num_angle= 10

f=open(constfname,'w')
nelem= num_exp_gauss *num_len_gauss +num_q_cos*num_exp_cos \
    +num_a0 *num_a2 *num_a6 *num_a10 +num_angle
#.....1st line
f.write(' {0:10d} {1:4d}\n'.format(nelem,nexp))
#.....repeat max_nexp times
for iexp in range(nexp):
    aexp= exps[iexp]
    #.....Gaussian
    dexp= (range_exp_gauss[1]-range_exp_gauss[0])/num_exp_gauss
    dlen= (range_len_gauss[1]-range_len_gauss[0])/num_len_gauss
    for i1 in range(1,num_exp_gauss+1):
        for i2 in range(1,num_len_gauss+1):
            egauss= range_exp_gauss[0] +i1*dexp 
            lgauss= range_len_gauss[0] +i1*dlen
            f.write(' {0:3d} {1:5.1f} {2:10.4f} {3:10.4f}\n'.format(1,
                                                                    aexp,
                                                                    egauss,
                                                                    lgauss))
    #.....cosine
    dq= (range_q_cos[1]-range_q_cos[0])/num_q_cos
    decs=(range_exp_cos[1]-range_exp_cos[0])/(num_exp_cos-1)
    for i0 in range(num_exp_cos):
        ecs= range_exp_cos[0] +i0*decs
        for i1 in range(1,num_q_cos+1):
            q= range_q_cos[0] +i1*dq
            f.write(' {0:3d} {1:5.1f} {2:10.4f} {3:10.4f}\n'.format(2,
                                                                  aexp,
                                                                  q,
                                                                  ecs))
    #.....polynomial
    d0= (range_a0_poly[1] -range_a0_poly[0])/(num_a0-1)
    d2= (range_a2_poly[1] -range_a2_poly[0])/(num_a2-1)
    d6= (range_a6_poly[1] -range_a6_poly[0])/(num_a6-1)
    d10=(range_a10_poly[1]-range_a10_poly[0])/(num_a10-1)
    for i0 in range(num_a0):
        a0= range_a0_poly[0] +i0*d0
        for i2 in range(num_a2):
            a2= range_a2_poly[0] +i2*d2
            for i6 in range(num_a6):
                a6= range_a6_poly[0] +i6*d6
                for i10 in range(num_a10):
                    a10= range_a10_poly[0] +i10*d10
                    f.write(' {0:3d} {1:5.1f} {2:10.4f}'.format(3,aexp,a0))
                    f.write(' {0:10.4f}'.format(a2))
                    f.write(' {0:10.4f}'.format(a6))
                    f.write(' {0:10.4f}'.format(a10))
                    f.write('\n')
    #.....angular
    dangle= (range_angle[1]-range_angle[0])/num_angle
    for i in range(1,num_angle+1):
        ang= range_angle[0] +i*dangle
        f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(4,
                                                    aexp,
                                                    ang))
f.close()
