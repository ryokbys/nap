#!/usr/bin/env python
"""
Make a polycrystal system.

Usage:
  make_polycrystal.py [options]

Options:
  -h, --help  Show this help message and exit.
  --ng NG     Number of grains in the system. [default: 3]
  -s, --size=SIZE
              Number of copies of the unit cell, int the format, nx,ny,nz. [default: 1,1,1]
  --struct STRUCTURE
              Crystal structure name. [default: fcc]
  -l, --lattice-constant=LATCONST
              Lattice constant of an axis. [default: 5.427]
  -o OUTFILE  Output file name. Format is detected automatically. [default: POSCAR]
  --two-dim   Whether or not 2D system. In the case of 2D, the system is a thin-sliced supercell with z-axis as the thinest direction.
"""
from __future__ import print_function

from docopt import docopt
import numpy as np
from random import random
import time

from nappy.napsys import NAPSystem
from nappy.atom import Atom
import cell_maker as cm


#...constants
RCUT = 3.5   # Angstrom
DMIN_RATE = 0.7

class Grain(object):
    
    def __init__(self,point,angle):
        self.point= point
        self.angle= angle
        self.rmat= self.get_rotation_matrix(self.angle)
        
    def get_rotation_matrix(self,angle):
        a= angle[0]
        b= angle[1]
        c= angle[2]
        rmat= np.zeros((3,3))
        #...Eular angle?
        # rmat[0,0]= cos(a)*cos(b)*cos(c) -sin(a)*sin(c)
        # rmat[0,1]=-cos(a)*cos(b)*sin(c) -sin(a)*cos(c)
        # rmat[0,2]= cos(a)*sin(b)
        # rmat[1,0]= sin(a)*cos(b)*cos(c) +cos(a)*sin(c)
        # rmat[1,1]=-sin(a)*cos(b)*sin(c) +cos(a)*sin(c)
        # rmat[1,2]= sin(a)*sin(b)
        # rmat[2,0]= -sin(b)*cos(c)
        # rmat[2,1]= sin(b)*sin(c)
        # rmat[2,2]= cos(b)
        #...yaw-pitch-rolling rotation
        rmx= np.zeros((3,3))
        rmy= np.zeros((3,3))
        rmz= np.zeros((3,3))
        rmx[0,:]= [ 1.0, 0.0,     0.0 ]
        rmx[1,:]= [ 0.0, np.cos(a), -np.sin(a) ]
        rmx[2,:]= [ 0.0, np.sin(a),  np.cos(a) ]
        rmy[0,:]= [ np.cos(b), 0.0, np.sin(b) ]
        rmy[1,:]= [ 0.0,    1.0,  0.0 ]
        rmy[2,:]= [-np.sin(b), 0.0, np.cos(b) ]
        rmz[0,:]= [ np.cos(c), -np.sin(c), 0.0 ]
        rmz[1,:]= [ np.sin(c),  np.cos(c), 0.0 ]
        rmz[2,:]= [ 0.0,     0.0,    1.0 ]
        rmat= np.dot(rmx,rmy)
        rmat= np.dot(rmat,rmz)
        return rmat

def anint(x):
    if x >= 0.5:
        x = 1.0
    elif x < -0.5:
        x =-1.0
    else:
        x = 0.0
    return x

def distance(p1,p2,two_dim=False):
    a = np.zeros(3)
    a[0]= p1[0]-p2[0] -anint(p1[0]-p2[0])
    a[1]= p1[1]-p2[1] -anint(p1[1]-p2[1])
    a[2]= p1[2]-p2[2] -anint(p1[2]-p2[2])
    if two_dim:
        a[2] = 0.0
    return np.linalg.norm(a)

def pbc(x):
    if x >= 1.0:
        x -= 1.0
    elif x < 0.0:
        x += 1.0
    return x

def shift_vector(two_dim=False):
    if not two_dim:
        nsv = 27
        sv= np.zeros((nsv,3))
        n=0
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    sv[n,0]= float(i)
                    sv[n,1]= float(j)
                    sv[n,2]= float(k)
                    n += 1
    else:
        nsv = 9
        sv= np.zeros((nsv,3))
        n= 0
        k= 0
        for i in range(-1,2):
            for j in range(-1,2):
                sv[n,0]= float(i)
                sv[n,1]= float(j)
                sv[n,2]= float(k)
                n += 1
    return sv,nsv

def uniq(lst):
    uniqarr = []
    for l in lst:
        if l not in uniqarr:
            uniqarr.append(l)
    return uniqarr

def make_polycrystal(grns,uc,n1,n2,n3,two_dim=False):
    """
    THIS ROUTINE IS NOT THAT UNIVERSAL.
    Each grain has to have neighboring grains within a supercell,
    otherwise there will be some unexpecting grain boundries.
    In order to do so, the system should be large enough and
    the number of grains should be large enough.
    """
    #...Calc the minimum bond distance in unit cell and use it as penetration depth
    dmin = 1.0e+30
    for i in range(uc.num_atoms()-1):
        for j in range(i+1,uc.num_atoms()):
            dij = uc.get_distance(i,j)
            dmin = min(dij,dmin)
    print(' Minimum bond distance in the unitcell: ',dmin)
    dmin = dmin *DMIN_RATE
    penetration_depth = dmin*2
    print(' Minimum bond distance allowed in the new system: ',dmin)
            
    sv,nsv= shift_vector(two_dim)
    # print(' nsv =',nsv)
    # for i in range(nsv):
    #     print(' i,sv[i]=',i,sv[i])
    system= NAPSystem(specorder=uc.specorder)
    system.set_lattice(uc.alc,uc.a1*n1,uc.a2*n2,uc.a3*n3)
    hmat= np.zeros((3,3))
    hmat[0]= system.a1 *system.alc
    hmat[1]= system.a2 *system.alc
    hmat[2]= system.a3 *system.alc
    hmati= np.linalg.inv(hmat)
    ix0 = -n1/2-1
    ix1 =  n1/2+2
    iy0 = -n2/2-1
    iy1 =  n2/2+2
    iz0 = -n3/2-1
    iz1 =  n3/2+2
    if two_dim:
        if n3 != 1:
            raise ValueError('n3 should be 1 in case two_dim is ON.')
        iz0 = 0
        iz1 = 1
    print(' x range = ',ix0,ix1)
    print(' y range = ',iy0,iy1)
    print(' z range = ',iz0,iz1)
    for ig in range(len(grns)):
        grain= grns[ig]
        rmat= grain.rmat  # Rotation matrix of the grain
        pi= grain.point   # Grain center in reduced coordinate
        api= np.dot(hmat,pi)  # Grain center in Cartessian coordinate
        print(' grain-ID = ',ig+1)
        for ix in range(ix0,ix1):
            # print('ix=',ix)
            for iy in range(iy0,iy1):
                for iz in range(iz0,iz1):
                    for m in range(len(uc.atoms)):
                        rt= np.zeros((3,))
                        rt[0]= (uc.atoms[m].pos[0]+ix)/n1
                        rt[1]= (uc.atoms[m].pos[1]+iy)/n2
                        rt[2]= (uc.atoms[m].pos[2]+iz)/n3
                        #...rt to absolute position
                        art= np.dot(hmat,rt)
                        #...Rotate
                        ari= np.dot(rmat,art)
                        #...Shift origin to the grain center
                        ari[0]= ari[0]+api[0]
                        ari[1]= ari[1]+api[1]
                        ari[2]= ari[2]+api[2]
                        #...check distance from all the grain points
                        di= distance(ari,api,two_dim)
                        isOutside= False
                        for jg in range(len(grns)):
                            gj= grns[jg]
                            for isv in range(nsv):
                                pj= gj.point
                                if jg == ig:
                                    if not two_dim and isv == 13:
                                        continue
                                    elif two_dim and isv == 4:
                                        continue
                                svi= sv[isv]
                                pj= pj +svi
                                apj = np.dot(hmat,pj)
                                dj= distance(ari,apj,two_dim)
                                if dj +penetration_depth < di:  # Allow some penetration here
                                    isOutside= True
                                    break
                            if isOutside:
                                break
                        if isOutside:
                            break
                        #...here ri is inside this grain, register it
                        atom= Atom()
                        #...Cartessian coord to reduced coord
                        ri = np.dot(hmati,ari)
                        ri[0]= pbc(ri[0])
                        ri[1]= pbc(ri[1])
                        ri[2]= pbc(ri[2])
                        atom.set_pos(ri[0],ri[1],ri[2])
                        atom.set_symbol(uc.atoms[m].symbol)
                        system.add_atom(atom)

    #...remove too-close atoms at the grain boundaries
    print(' Making pair list in order to remove close atoms...')
    print(' Number of atoms: ',system.num_atoms())
    system.make_pair_list(RCUT)
    system.write('POSCAR_tmp')
    short_pairs = []
    # dmin2= dmin**2
    # xij= np.zeros((3,))
    print(' Making the list of SHORT pairs...')
    for ia in range(system.num_atoms()):
        # ai= system.atoms[ia]
        # pi= ai.pos
        nlst= system.nlspr[ia]
        lst= system.lspr[ia]
        for j in range(nlst):
            ja= lst[j]
            dij = system.get_distance(ia,ja)
            if dij < dmin:
                short_pairs.append((ia,ja,dij))
            # aj= system.atoms[ja]
            # pj= aj.pos
            # xij[0]= pj[0]-pi[0] -anint(pj[0]-pi[0])
            # xij[1]= pj[1]-pi[1] -anint(pj[1]-pi[1])
            # xij[2]= pj[2]-pi[2] -anint(pj[2]-pi[2])
            # xij= np.dot(hmat,xij)
            # d2= xij[0]**2 +xij[1]**2 +xij[2]**2
            # if d2 < dmin2:
            #     if not ia in ls_remove:
            #         ls_remove.append(ia)
            #     elif not ja in ls_remove:
            #         ls_remove.append(ja)
        # print('ia,len(ls_remove)=',ia,len(ls_remove))

    print(' Number of short pairs: ',len(short_pairs))

    #...Remove only relevant atoms, not all the atoms in the short_pairs.
    ls_remove = []
    ls_not_remove = []
    for pair in short_pairs:
        ia = pair[0]
        ja = pair[1]
        if ia not in ls_not_remove and ja not in ls_not_remove:
            ls_remove.append(ia)
            ls_not_remove.append(ja)
        elif ia not in ls_not_remove:
            ls_remove.append(ia)
        elif ja not in ls_not_remove:
            ls_remove.append(ja)
        else:  # Both atoms are already in not_remove list, which should be avoided.
            ls_not_remove.remove(ia)
            ls_remove.append(ia)
            ls_not_remove.append(ja)
    #...Remove double registered IDs
    ls_remove = uniq(ls_remove)
    print(' Number of to be removed atoms: ',len(ls_remove))
            
    # print(' Number of to be removed atoms: ',len(ls_remove))
    #...one of two will survive
    # print(' One of two too-close atoms will survive...')
    # count= [ ls_remove.count(ls_remove[i]) for i in range(len(ls_remove))]
    # for i in range(0,len(ls_remove),2):
    #     if count[i] > count[i+1]:
    #         ls_remove[i+1]= -1
    #     elif count[i] < count[i+1]:
    #         ls_remove[i]= -1
    #     else:
    #         n= int(random()*2.0) # 0 or 1
    #         ls_remove[i+n]= -1
    ls_remove.sort()
    # for ia in range(len(ls_remove)-1,-1,-1):
    #     n= ls_remove[ia]
    #     if ia != len(ls_remove)-1:
    #         if n == nprev: continue
    #     system.atoms.pop(n)
    #     nprev= n
    for i in reversed(range(len(ls_remove))):
        ia = ls_remove[i]
        system.atoms.pop(ia)
    return system


if __name__ == '__main__':

    args = docopt(__doc__)
    ngrain = int(args['--ng'])
    nx,ny,nz = [ int(x) for x in args['--size'].split(',') ]
    ofname = args['-o']
    latconst= float(args['--lattice-constant'])
    struct = args['--struct']
    two_dim = args['--two-dim']

    print(' {0:=^72}'.format(' make_polycrystal_py '))
    t0= time.time()

    print(" number of grains= {0:3d}".format(ngrain))
    print(" number of unit cells= {0:3d} {1:3d} {2:3d}".format(nx,ny,nz))

    makestruct = None
    if struct == 'bcc':
        makestruct = cm.make_bcc
    elif struct == 'bcc110':
        makestruct = cm.make_bcc110
    elif struct == 'bcc111':
        makestruct = cm.make_bcc111
    elif struct == 'fcc':
        makestruct = cm.make_fcc
    elif struct == 'fcc110':
        makestruct = cm.make_fcc110
    else:
        raise ValueError('No such structure available: ',struct)

    grains= []
    for i in range(ngrain):
        pi= np.zeros((3,))
        ai= np.zeros((3,))
        pi[0]= random()
        pi[1]= random()
        pi[2]= random()
        if two_dim:
            ai[0]= 0.0
            ai[1]= 0.0
            ai[2]= random()*np.pi*2 -np.pi
        else:
            ai[0]= random()*np.pi*2 -np.pi
            ai[1]= random()*np.pi/2 -np.pi/2
            ai[2]= random()*np.pi*2 -np.pi
        print(' point=',pi)
        print(' angle=',ai)
        gi= Grain(pi,ai)
        grains.append(gi)
    unitcell= makestruct(latconst)
    unitcell.write('POSCAR_uc')
    system= make_polycrystal(grains,unitcell,nx,ny,nz,two_dim)
    system.write(ofname)

    print(' Elapsed time = {0:12.2f}'.format(time.time()-t0))
    print(' Wrote a file: {0:s}'.format(ofname))
