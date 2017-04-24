#!/usr/bin/env python
"""
Make in.const.NN and in.params.NN to be used in NN potential.
And the file storing combination information, 
in.comb.NN, is also written.
"""

import random,math

#=========================================================== Constants
#.....cutoff radius in Angstrom
rcut= 5.8
rc3 = 3.5
#.....number of species
nsp= 4
#.....number of hidden layers
nl = 1
#.....num of nodes in a layer
nhl= [0,30]
#.....min,max of parameters
pmin= -0.1
pmax=  0.1
#.....num of eta in Gaussian symmetry function, f(r)=exp(-eta*(dij-rs)**2)
type_gauss= 1
#.....num of Rs in 2-body symmetry function
#rrs=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
#rrs=[0.5, 1.0, 1.5, 2.0, 2.5]
nrs= 20
rmin= 1.5
dr= ((rcut-0.1)-rmin)/(nrs-1)
rrs=[ rmin+dr*i for i in range(nrs) ]
#rrs=[0.0]
nr= 1
#dr= (10.0-0.01)/(nr-1)
#reta=[ 0.1*10.0**i for i in range(nr) ]
#reta=[ 0.01+dr*i for i in range(nr) ]
# reta= [0.01, 0.1, 1.0, 10.0]
# reta=[0.01, 0.1, 1.0]
#reta= [ 1.0/2/(2*dr)**2 ]
reta= [10.0]

# #.....num of k in cosine func, f(r)= (1d0+cos(k*r))
# type_cos= 2
# nk= 50
# drk= (rcut-1.0)*math.pi/rcut/(nk-1)
# rk= [ math.pi/rcut +drk*i for i in range(nk) ]
# r2pi= 2.0*math.pi/rcut
#rk= [0.5, 1.0, 2.0, 3.0, 4.0]
rk= []
# #.....num of Morse potential
# type_morse= 4
# nm1= 4
# dnm1= (10.0-0.1)/(nm1-1)
# rm1= [ 0.1 +dnm1*i for i in range(nm1) ]
# nm2= 20
# dnm2= (rcut-1.0-0.5)/(nm2-1)
# rm2= [ 0.5 +dnm2*i for i in range(nm2) ]
#.....num of 3-body angular symmetry functions (cosine value)
type_angle= 101
nang= 3
dang= 1.0/(nang-1)
rsf3= [ dang*i for i in range(nang) ]
#rsf3= []
#rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2, 2.0/3, 3.0/5]
# rsf3=[0.0, 1.0/5, 1.0/3, 1.0/2, 1.0]
#choose_triplet = []
choose_triplet = [(3, 1, 1), (3, 1, 2), (3, 1, 3), (3, 2, 2), (3, 2, 3), (3, 3, 3),]

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
            if choose_triplet and not (tmp[0],pair[0],pair[1]) in choose_triplet:
                continue
            triplets.append(tmp+pair)
    return pairs,triplets
    

#========================================================= main routine
if __name__ == "__main__":

    if not nl in (1,2):
        print " [Error] nl is not 1 nor 2, nl=",nl
        exit()

    print 'rrs = ',rrs
    print 'reta = ',reta
    print 'rk = ',rk
    print 'rsf3 = ',rsf3

    #....compute num of combinations
    ncmb2= nsp+ ncomb(nsp,2)/2
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
        print 'Since len(triplets) != ncmb3, set ncmb3 = len(triplets)'
        ncmb3 = len(triplets)
    
    #.....num of 2-body Gaussian-type symmetry functions
    nsf2 = len(reta)*len(rrs)
    nsf2+= len(rk)
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
    g.write(' {0:6d} {1:10.4f} {2:6.2f}\n'.format(nc,rcut,rc3))
    for ic in range(nc):
        g.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
        g.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
    g.close()
