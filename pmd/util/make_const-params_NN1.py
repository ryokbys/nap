#!/opt/local/bin/python

u"""
Make in.const.NN1 and in.params.NN1 to be used in NN1 potential.
"""

import sys,os,random

constfname='in.const.NN1'
paramfname='in.params.NN1'

#.....cutoff radius in Angstrom
rcut= 7.0
#.....min,max of parameters
pmin= -1.0
pmax=  1.0
#.....num of eta in 2-body symmetry function
neta= 5
reta=[0.01, 0.05, 0.1, 0.2, 0.4, 1.0]
#.....num of Rs in 2-body symmetry function
nrs= 3
rrs=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
#.....num of 3-body angular symmetry functions
nsf3= 5
rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2, 1.0]
#.....num of nodes in a layer
nhl1= 10

#.....num of 2-body Gaussian-type symmetry functions
nsf2= neta*nrs

f= open(constfname,'w')
nsf= nsf2+nsf3
f.write(' {0:5d} {1:5d}\n'.format(nsf,nhl1))
#.....2-body Gaussian-type
for ieta in range(neta):
    for irs in range(nrs):
        f.write(' 1 {0:10.4f} {1:10.4f}\n'.format(reta[ieta],rrs[irs]))
#.....3-body
for isf3 in range(nsf3):
    f.write(' 2 {0:10.4f}\n'.format(rsf3[isf3]))
f.close()

g= open(paramfname,'w')
nc= (nsf+1)*nhl1 +(nhl1+1)
g.write(' {0:6d} {1:10.4f}\n'.format(nc,rcut))
for ic in range(nc):
    g.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
    g.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
g.close()
