#!/opt/local/bin/python

u"""
Make in.const.linreg and in.params.linreg to be used in linreg potential.
"""

import sys,os,random

constfname='in.const.linreg'
paramfname='in.params.linreg'

#.....cutoff radius in Angstrom
rcut= 3.0
#.....min,max of parameters
pmin= -0.00001
pmax=  0.00001
#.....exponent of the basis func
rexp=[0.5, 1.0]
nexp=len(rexp)
#.....Gaussian
rexp_gauss=[1.0, 2.0]
nexp_gauss= len(rexp_gauss)
rlen_gauss=[0.5,0.75,1.0,1.25,1.5,1.75,2.0]
nlen_gauss= len(rlen_gauss)
rsft_gauss=[0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4]
nsft_gauss= len(rsft_gauss)
#.....cosine
rq_cos=[2.0,4.0,6.0,8.0,10.0]
nq_cos= len(rq_cos)
rexp_cos=[0.0, 1.0]
nexp_cos= len(rexp_cos)
#.....polynomial
ra0_poly= [1.0]
ra1_poly= [1.0]
ra2_poly= [1.0]
ra4_poly= [1.0]
ra6_poly= [1.0]
ra8_poly= [1.0]
ra10_poly= [1.0]
ra12_poly= [1.0]
na1= len(ra1_poly)
na2= len(ra2_poly)
na4= len(ra4_poly)
na6= len(ra6_poly)
na8= len(ra8_poly)
na10= len(ra10_poly)
na12= len(ra12_poly)
#.....angular
#rangle=[0.0, 0.25, 1.0/3, 0.5, 0.75, 1.0]
rangle=[0.25, 1.0/3, 0.5, 0.75]
nangle= len(rangle)

f=open(constfname,'w')
nelem= nexp_gauss*nlen_gauss*nsft_gauss \
    +nq_cos*nexp_cos \
     +nangle \
     +na1 +na2 +na4 +na6 +na8 +na10 +na12
#.....1st line
f.write(' {0:10d} {1:5d}\n'.format(nelem*nexp,nexp))
#.....repeat max_nexp times
for iexp in range(nexp):
    aexp= rexp[iexp]
    #.....Gaussian
    for egauss in rexp_gauss:
        for lgauss in rlen_gauss:
            for sgauss in rsft_gauss:
                f.write(' {0:3d} {1:5.1f} '.format(1,aexp))
                f.write('{0:10.4} {1:10.4f} {2:10.4f}\n'.format(egauss,
                                                                lgauss,
                                                                sgauss))
    #.....cosine
    for ecs in rexp_cos:
        for qcs in rq_cos:
            f.write(' {0:3d} {1:5.1f} {2:10.4f} {3:10.4f}\n'.format(2,
                                                                    aexp,
                                                                    qcs,
                                                                    ecs))
    #.....angular
    for ang in rangle:
        f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(3,
                                                    aexp,
                                                    ang))
    #.....polynomial
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(4,aexp,ra1_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(5,aexp,ra2_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(6,aexp,ra4_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(7,aexp,ra6_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(8,aexp,ra8_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(9,aexp,ra10_poly[0]))
    f.write(' {0:3d} {1:5.1f} {2:10.4f}\n'.format(10,aexp,ra12_poly[0]))

f.close()

#.....output parameter file
g=open(paramfname,'w')
g.write(' {0:6d} {1:10.4f}\n'.format(nelem*nexp,rcut))
for ie in range(nelem*nexp):
    g.write(' {0:15.6e}'.format(random.uniform(pmin,pmax)))
    g.write(' {0:13.4e} {1:13.4e}\n'.format(pmin,pmax))
g.close()
