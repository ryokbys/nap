#!/opt/local/bin/python

u"""
Make in.const.NN1 and in.params.NN1 to be used in NN1 potential.
And the file storing combination information, 
in.comb.NN1, is also written.
"""

import sys,os,random,math

#=========================================================== Functions
def comb(n,m):
    '''
    Calculate nCm.
    '''
    return math.factorial(n)/math.factorial(m)

def make_combination(nsp,fname='in.comb.NN1'):
    '''
    Make combinations and write them into file, in.comb.NN1.
    '''
    f= open(fname,'w')
    #.....2-body
    n=0
    pairs=[]
    for i in range(1,nsp+1):
        for j in range(i,nsp+1):
            pairs.append([i,j])
            n += 1
            f.write(' {0:4d} {1:3d} {2:3d}\n'.format(n,i,j))
    
    #.....3-body
    n= 0
    for i in range(1,nsp+1):
        for pair in pairs:
            n += 1
            f.write(' {0:4d} {1:3d} {2:3d} {3:3d}\n'.format(n,i, \
                                                            pair[0], \
                                                            pair[1] ))
    f.close()

#=========================================================== Constants
constfname='in.const.NN1'
paramfname='in.params.NN1'

#.....cutoff radius in Angstrom
rcut= 5.0
#.....number of species
nsp= 2
#.....min,max of parameters
pmin= -0.01
pmax=  0.01
#.....num of eta in 2-body symmetry function
reta=[0.5, 1.0]
#.....num of Rs in 2-body symmetry function
rrs=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
#.....num of 3-body angular symmetry functions (cosine value)
rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2, 2.0/3, 3.0/5]
#.....num of nodes in a layer
nhl1= 2

#=========================================================== Routines

#..........................................compute num of combinations
ncmb2= nsp+ comb(nsp,2)
ncmb3= ncmb2*nsp

make_combination(nsp)

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
#nc= (nsf+1)*nhl1 +(nhl1+1)
nsf= nsf2*ncmb2 +nsf3*ncmb3
nc= nsf*nhl1 +nhl1
g.write(' {0:6d} {1:10.4f}\n'.format(nc,rcut))
for ic in range(nc):
    g.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
    g.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
g.close()
