#!/opt/local/bin/python

u"""
Make in.const.NN and in.params.NN to be used in NN potential.
And the file storing combination information, 
in.comb.NN, is also written.
"""

import sys,os,random,math
import optparse

#=========================================================== Constants
#.....cutoff radius in Angstrom
rcut= 3.5
#.....number of species
nsp= 2
#.....number of hidden layers
nl = 1
#.....num of nodes in a layer
nhl= [0,2]
#.....min,max of parameters
pmin= -1.0
pmax=  1.0
#.....num of eta in Gaussian symmetry function, f(r)=exp(-eta*(dij-rs)**2)
type_gauss= 1
reta=[0.5, 1.0, 1.5, 2.0, 2.5]
#.....num of Rs in 2-body symmetry function
#rrs=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
#rrs=[0.5, 1.0, 1.5, 2.0, 2.5]
rrs=[0.0]
#.....num of k in cosine func, f(r)= (1d0+cos(k*r))
type_cos= 2
r2pi= 2.0*math.pi/rcut
rk= [r2pi, r2pi*1.5]
#.....num of a in polynomial func, f(r)= 1.0/r**a
type_poly= 3
rpoly= [1.0, 2.0, 3.0, 6.0, 9.0, 12.0]
#.....num of 3-body angular symmetry functions (cosine value)
type_angle= 101
#rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2, 2.0/3, 3.0/5]
rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2]

constfname='in.const.NN'
paramfname='in.params.NN'
combfname='in.comb.NN'


#=========================================================== Functions
def ncomb(n,m):
    '''
    Calculate nCm.
    '''
    return math.factorial(n)/math.factorial(m)

def make_combination(nsp,fname):
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
            f.write(' {0:3d} {1:3d} {2:4d}\n'.format(i,j,n))
    
    #.....3-body
    n= 0
    for i in range(1,nsp+1):
        for pair in pairs:
            n += 1
            f.write(' {0:3d} {1:3d} {2:3d} {3:4d}\n'.format(i, \
                                                            pair[0], \
                                                            pair[1],n ))
    f.close()

def get_comb(nsp):
    """
    Make combinations.
    """
    pairs= []
    for i in range(1,nsp+1):
        for j in range(i,nsp+1):
            pairs.append([i,j])
    triplets= []
    for i in range(1,nsp+1):
        tmp= [i]
        for pair in pairs:
            triplets.append(tmp+pair)
    return pairs,triplets
    

#========================================================= main routine
if __name__ == "__main__":

    if not nl in (1,2):
        print " [Error] nl is not 1 nor 2, nl=",nl
        exit()

    #....compute num of combinations
    ncmb2= nsp+ ncomb(nsp,2)
    ncmb3= ncmb2*nsp

    print ' ncmb2, ncmb3 = ',ncmb2,ncmb3

    pairs,triplets= get_comb(nsp)
    print ' pairs=',pairs
    print ' triplets=',triplets
    if len(pairs) != ncmb2:
        print '[Error] len(pairs) != ncmb2'
        print 'len(pairs),ncmb2=',len(pairs),ncmb2
        exit()
    if len(triplets) != ncmb3:
        print '[Error] len(triplets) != ncmb3'
        print 'len(triplets),ncmb3=',len(triplets),ncmb3
        exit()
    
    #.....num of 2-body Gaussian-type symmetry functions
    nsf2 = len(reta)*len(rrs)
    nsf2+= len(rk)
    nsf2+= len(rpoly)
    nsf3= len(rsf3)

    print ' nsf2, nsf3 = ',nsf2,nsf3
    
    f= open(constfname,'w')
    nsf= nsf2+nsf3
    nhl[0]= nsf2*ncmb2 +nsf3*ncmb3
    #nhl[0]= nsf
    f.write(' {0:5d} {1:5d}'.format(nl,nsp))
    for il in range(nl+1):
        f.write(' {0:5d}'.format(nhl[il]))
    f.write('\n')

    for pair in pairs:
        ia= pair[0]
        ja= pair[1]
        for eta in reta: # Gaussian
            for rs in rrs:
                f.write(' {0:3d}'.format(type_gauss) \
                        +' {0:3d} {1:3d}'.format(ia,ja) \
                        +' {0:10.4f} {1:10.4f}\n'.format(eta,rs))
        for k in rk:  # cosine
            f.write(' {0:3d}'.format(type_cos) \
                    +' {0:3d} {1:3d}'.format(ia,ja) \
                    +' {0:10.4f}\n'.format(k))
        for p in rpoly:  # polynomial
            f.write(' {0:3d}'.format(type_poly) \
                    +' {0:3d} {1:3d}'.format(ia,ja) \
                    +' {0:10.4f}\n'.format(p))
    for triple in triplets:
        ia= triple[0]
        ja= triple[1]
        ka= triple[2]
        for sf3 in rsf3: # 3-body
            f.write(' {0:3d}'.format(type_angle) \
                    +' {0:3d} {1:3d} {2:3d}'.format(ia,ja,ka) \
                    +' {0:10.4f}\n'.format(sf3))
    f.close()
    
    g= open(paramfname,'w')
    #nc= (nsf+1)*nhl1 +(nhl1+1)
    #nsf= nsf2*ncmb2 +nsf3*ncmb3
    #nhl[0]= nsf2*ncmb2 +nsf3*ncmb3
    if nl == 1:
        nc= nhl[0]*nhl[1] +nhl[1]
    elif nl == 2:
        nc= nhl[0]*nhl[1] +nhl[1]*nhl[2] +nhl[2]
    g.write(' {0:6d} {1:10.4f}\n'.format(nc,rcut))
    for ic in range(nc):
        g.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
        g.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
    g.close()
