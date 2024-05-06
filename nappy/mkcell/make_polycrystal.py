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
  --min-gdist=MIN_GDIST
              Minimum grain-center distance in reduced unit. Note that the x,y-lengths are assumed to be comparable, 
              otherwise the distance in the reduced unit looses its sense. [default: 0.1]
  --angle-range=ANGLE_RANGE
              Range of angle for grains to be distributed within. [default: 180.0]
"""
from docopt import docopt
import numpy as np
from random import random
import time

from nappy.napsys import NAPSystem
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
    nsys= NAPSystem(specorder=uc.specorder)
    nsys.set_lattice(uc.alc,uc.a1*n1,uc.a2*n2,uc.a3*n3)
    hmat = nsys.get_hmat()
    hmati = nsys.get_hmat_inv()
    nmax = n1*n2*n3 *uc.num_atoms()
    sidsl = np.zeros(nmax,dtype=int)
    symsl = []
    possl = np.zeros((nmax,3))
    velsl = np.zeros((nmax,3))
    frcsl = np.zeros((nmax,3))
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
    inc = 0
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
                    for m in range(uc.num_atoms()):
                        sidt = uc.get_atom_attr(m,'sid')
                        rt= np.zeros((3,))
                        pm = uc.get_atom_attr(m,'pos')
                        rt[0]= (pm[0]+ix)/n1
                        rt[1]= (pm[1]+iy)/n2
                        rt[2]= (pm[2]+iz)/n3
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
                        #...Cartessian coord to reduced coord
                        ri = np.dot(hmati,ari)
                        ri[0]= pbc(ri[0])
                        ri[1]= pbc(ri[1])
                        ri[2]= pbc(ri[2])
                        sidsl[inc] = sidt
                        possl[inc] = ri
                        velsl[inc,:] = 0.0
                        frcsl[inc,:] = 0.0
                        symsl.append(nsys.specorder[sidt-1])
                        inc += 1
                        if inc > nmax:
                            raise ValueError('inc > nmax')
    #...Create filled arrays from non-filled ones
    poss = np.array(possl[:inc])
    vels = np.array(velsl[:inc])
    frcs = np.array(frcsl[:inc])
    nsys.add_atoms(symsl,poss,vels,frcs)

    #...remove too-close atoms at the grain boundaries
    print(' Making pair list in order to remove close atoms...')
    print(' Number of atoms: ',nsys.num_atoms())
    nsys.make_pair_list(RCUT)
    nsys.write('POSCAR_orig')
    short_pairs = []
    # dmin2= dmin**2
    # xij= np.zeros((3,))
    print(' Making the list of SHORT pairs...')
    for ia in range(nsys.num_atoms()):
        lst= nsys.get_atom_attr(ia,'lspr')
        for j in range(len(lst)):
            ja= lst[j]
            if ja > ia:
                continue
            dij = nsys.get_distance(ia,ja)
            if dij < dmin:
                short_pairs.append((ia,ja,dij))

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

    nsys.remove_atoms(*ls_remove)
    return nsys

def get_grains(ng,gdmin,angrange0,angrange1,two_dim):
    """
    Get specified number of grains with conditions of minimum distance and angle range.
    """
    dang = (angrange1-angrange0) /180.0 *np.pi /ng
    grains= []
    ig = 0
    dmin = 1e+30
    while True:
        if ig >= ng: break
        pi= np.zeros((3,))
        ai= np.zeros((3,))
        pi[0]= random()
        pi[1]= random()
        if two_dim:
            too_close = False
            dminl = 1e+30
            for i,gi in enumerate(grains):
                dlt = pi - gi.point
                dlt = dlt - np.round(dlt)
                d = np.sqrt( dlt[0]**2+dlt[1]**2 )
                dminl = min(d,dminl)
                if d < gdmin:
                    too_close = True
                    break
            if too_close:
                continue
            dmin = min(dminl,dmin)
            pi[2]= 0.0
            ai[0]= 0.0
            ai[1]= 0.0
            ai[2]= angrange0 +dang*ig +random()*dang
                    
        else:
            pi[2]= random()
            too_close = False
            for gi in grains:
                dlt = pi - gi.point
                dlt = dlt - np.round(dlt)
                d = np.sqrt( dlt[0]**2+dlt[1]**2+dlt[2]**2 )
                if d < gdmin:
                    too_close
                    break
            if too_close:
                continue
            ai[0]= random()*np.pi*2 -np.pi
            ai[1]= random()*np.pi/2 -np.pi/2
            ai[2]= random()*np.pi*2 -np.pi
        print(' point,angle =',pi,ai)
        gi= Grain(pi,ai)
        grains.append(gi)
        ig += 1
    print(' Minimum distance between grains and limit = ',dmin,gdmin)
    return grains

if __name__ == '__main__':

    args = docopt(__doc__)
    ngrain = int(args['--ng'])
    nx,ny,nz = [ int(x) for x in args['--size'].split(',') ]
    ofname = args['-o']
    latconst= float(args['--lattice-constant'])
    struct = args['--struct']
    two_dim = args['--two-dim']
    gdistmin = float(args['--min-gdist'])
    anglerange = float(args['--angle-range'])

    #...Set angle range within [-180.0:180.0]
    angrange1 = anglerange % 360.0
    angrange0 = max(0.0, angrange1-180.0)
    angrange1 = angrange1 -angrange0

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


    grains= get_grains(ngrain,gdistmin,angrange0,angrange1,two_dim)
    uc= makestruct(latconst)
    uc.write('POSCAR_uc')
    nsys= make_polycrystal(grains,uc,nx,ny,nz,two_dim)
    nsys.write(ofname)

    print(' Elapsed time = {0:12.2f}'.format(time.time()-t0))
    print(' Wrote a file: {0:s}'.format(ofname))
