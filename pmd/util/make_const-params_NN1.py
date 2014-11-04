#!/opt/local/bin/python

u"""
Make in.const.NN1 and in.params.NN1 to be used in NN1 potential.
"""

import sys,os,random

constfname='in.const.NN1'
paramfname='in.params.NN1'

#.....cutoff radius in Angstrom
rcut= 5.0
#.....min,max of parameters
pmin= -1.0
pmax=  1.0
#.....num of eta in 2-body symmetry function
#reta=[0.01, 0.05, 0.1, 0.2, 0.4, 1.0]
reta=[1.0]
#reta=[0.05]
#.....num of Rs in 2-body symmetry function
#rrs=[1.0, 1.5, 2.0, 2.5, 3.0]
rrs=[1.0, 2.0, 3.0]
#.....num of 3-body angular symmetry functions
#rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2]
rsf3=[1.0/5, 1.0/3]
#.....num of nodes in a layer
nhl1= 3
#nhl1= 3

#.....num of 2-body Gaussian-type symmetry functions
nsf2= len(reta)*len(rrs)
nsf3= len(rsf3)

f= open(constfname,'w')
nsf= nsf2+nsf3
f.write(' {0:5d} {1:5d}\n'.format(nsf,nhl1))
#.....2-body Gaussian-type
for eta in reta:
    for rs in rrs:
        f.write(' 1 {0:10.4f} {1:10.4f}\n'.format(eta,rs))
#.....3-body
for sf3 in rsf3:
    f.write(' 2 {0:10.4f}\n'.format(sf3))
f.close()

g= open(paramfname,'w')
nc= (nsf+1)*nhl1 +(nhl1+1)
g.write(' {0:6d} {1:10.4f}\n'.format(nc,rcut))
for ic in range(nc):
    g.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
    g.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
g.close()
