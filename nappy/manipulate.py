#!/usr/bin/env python
"""
Functions for manipulating napsys objects.
This script itself does nothing by calling from commandline.
Import it and use some functions from other scripts or notebooks.

Usage:
  manipulate.py [options]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

from docopt import docopt
import numpy as np
import copy
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "200424"

def substitute(nsys0,spc1,spc2,num=1):
    """
    Substitute an atom of SPC1 in the system NSYS with SPC2.

    Parameters
    ----------
    nsys0 : NAPSystem object
           System in which atoms are to be replaced. (not inplace)
    spc1 : integer or string
           An atom ID or atom species to be replaced.
           If it is integer, the atom specified by SPC1 is replaced.
           If it is string, NUM atoms of species SPC1 are to be replaced with atoms of SPC2.
    spc2 : integer or string
           An atom ID to be replaced with spcs1, which is available only if spcs1 is integer.
           Atom species that replace atoms of spc1.
    num : integer, optional
           Number of atoms to be replaced.
    
    Returns
    -------
    nsys_out : NAPSystem object
           Result of the system in which atoms are substituted.
    """

    if type(spc1) not in (int,str):
        raise TypeError('SPC1 must be either int or str.')
    if type(spc2) not in (int,str):
        raise TypeError('SPC2 must be either int or str.')
    if type(spc2) is int and type(spc1) is not int:
        raise ValueError('SPC2 must be string if SPCS1 is string.')
    if not isinstance(nsys0,nappy.napsys.NAPSystem):
        raise ValueError('NSYS0 is not an instance of NAPSystem.')
    
    import copy
    nsys = copy.deepcopy(nsys0)

    if type(spc1) is int:
        if spc1 > nsys.num_atoms():
            raise ValueError('An integer SPC1 is greater than num of atoms in NSYS0')
        sid1 = nsys.atoms.sid[spc1]
        if type(spc2) is int:
            sid2 = nsys.atoms.sid[spc2]
            if sid1 == sid2:
                # do nothing
                return nsys
            else:
                nsys.atoms.at[spc1,'sid'] = sid2
                nsys.atoms.at[spc2,'sid'] = sid1
                return nsys
        else:
            if spc2 not in nsys.specorder:
                nsys.specorder.append(spc2)
            sid2 = nsys.specorder.index(spc2) +1
            if sid1 == sid2:
                # do nothing
                return nsys
            else:
                nsys.atoms.at[spc1,'sid'] = sid2
                return nsys
    else:  # type(spc1) is str
        import random
        if spc1 not in nsys.specorder:
            raise ValueError('The species {0:s} is not in the system.'.format(spc1))
        sid1 = nsys.specorder.index(spc1) +1
        if spc2 not in nsys.specorder:
            nsys.specorder.append(spc2)
        sid2 = nsys.specorder.index(spc2) +1
        indices = nsys.atoms[nsys.atoms.sid == sid1].index.tolist()
        for n in range(num):
            if len(indices) == 0:
                break
            idx = random.choice(indices)
            nsys.atoms.at[idx,'sid'] = sid2
            indices.remove(idx)
        return nsys
    
    raise RuntimeError('Something is wrong.')

def insert(nsys0,num=1,spc=None,rule='maxmindist',num_trial=100):
    """
    Add atoms into the sytem of reasonable positions.
    Look for the reasonable position by Monte-Carlo simulation with Levy flight motion.
    
    Parameters
    ----------
    nsys0 : NAPSystem object
          The system into which atoms are inserted (not inplace).
    num : integer
          Number of atoms to be inserted.
    rule : str
          The rule to be applied when placing atoms into the system.
          - maxmindist: maximize min distance
    num_trial : integer
          Number of trials of insertion positions.

    Returns
    -------
    nsys : NAPSystem object
         The system in which atoms are inserted.
    """
    if not isinstance(nsys0,nappy.napsys.NAPSystem):
        raise ValueError('NSYS0 must be an instance of NAPSystem.')
    if type(spc) is not str:
        raise ValueError('SPC must be string.')
    if num < 1:
        raise ValueError('NUM must be greater than 0.')
    if num_trial < 1:
        raise ValueError('NUM_TRIAL must be greater than 0.')

    import copy
    #...Levy flight, see RK's memo about Cuckoo search
    from scipy.special import gamma
    beta = 1.5
    betai = 1.0/beta
    vsgm = 1.0
    usgm = (gamma(1+beta)*np.sin(np.pi*beta/2)/ \
            gamma((1+beta)/2)*beta*2.0**((beta-1)/2))**betai

    
    nsys = copy.deepcopy(nsys0)
    if spc not in nsys.specorder:
        nsys.specorder.append(spc)
    sid = nsys.specorder.index(spc) +1
    hmat = nsys.get_hmat()

    for n in range(num):
        pi = np.random.rand(3)
        pi_prev = copy.copy(pi)
        d2min_prev = 0.0
        spos = nsys.get_scaled_positions()
        for it in range(num_trial):
            #...Move the position using Levy flight
            for ixyz in range(3):
                p = pi_prev[ixyz]
                u = np.random.normal() *usgm
                v = max( abs(np.random.normal()*vsgm), 1.0e-8 )
                w = u /v**betai
                zeta = 0.01 *w
                p += zeta*np.random.normal()
                if p < 0.0:
                    p += 1.0
                elif p >= 1.0:
                    p -= 1.0
                pi[ixyz] = p

            #...Calc minimum distance
            d2min = 1.0e+30
            for j in range(len(nsys.atoms)):
                pj = spos[j]
                xij = pj -pi
                xij = xij -np.round(xij)
                rij = np.dot(hmat,xij)
                dij2 = rij[0]**2 +rij[1]**2 +rij[2]**2
                d2min = min(d2min,dij2)

            #...Decide whether or not to employ the new position
            # print('pi=',pi,', d2min,d2min_prev=',d2min,d2min_prev)
            if d2min > d2min_prev:
                pi_prev = copy.copy(pi)
                d2min_prev = d2min


        #...Insert the atom at the position
        symbols = [spc]
        poss = [pi_prev,]
        vels = [[0., 0., 0.],]
        frcs = [[0., 0., 0.],]
        # print('Atom added to {0:6.3f} {1:6.3f} {2:6.3f}, dist={3:.3f}'.format(*pi_prev,np.sqrt(d2min)))
        nsys.add_atoms(symbols,poss,vels,frcs)

    return nsys

def replicate(nsys0,n1o,n2o,n3o,n1m=0,n2m=0,n3m=0):
    """
    Return the system multiplied by n1o,n2o,n3o.
    """
    import pandas as pd
    nsys = copy.deepcopy(nsys0)
    #...Convert to int
    n1 = int(n1o)
    n2 = int(n2o)
    n3 = int(n3o)
    if n1 == 0: n1 = 1
    if n2 == 0: n2 = 1
    if n3 == 0: n3 = 1
    if n1 == n2 == n3 == 1:
        return None
    #...unit vectors to be repeated
    m1 = n1-n1m
    m2 = n2-n2m
    m3 = n3-n3m
    nsys.a1= nsys.a1*m1
    nsys.a2= nsys.a2*m2
    nsys.a3= nsys.a3*m3
    #n123= m1*m2*m3
    maxsid = nsys.atoms.sid.max()
    #natm0= nsys.num_atoms()
    # atoms0= copy.copy(nsys.atoms)
    newnatm = len(nsys.atoms) *m1*m2*m3
    newsids = [ 0 for i in range(newnatm) ]
    newposs = np.zeros((newnatm,3))
    newvels = np.zeros((newnatm,3))
    newfrcs = np.zeros((newnatm,3))
    colnames = list(nsys.atoms.columns)
    #...Labels except (sid,pos,vel,frc) are all auxiliary data
    auxnames = colnames.copy()
    auxnames.remove('sid')
    auxnames.remove('x')
    auxnames.remove('y')
    auxnames.remove('z')
    auxnames.remove('vx')
    auxnames.remove('vy')
    auxnames.remove('vz')
    auxnames.remove('fx')
    auxnames.remove('fy')
    auxnames.remove('fz')
    newauxs = {}
    for auxname in auxnames:
        newauxs[auxname] = []
    inc = 0
    poss = nsys.get_scaled_positions()
    vels = nsys.get_scaled_velocities()
    frcs = nsys.get_scaled_forces()
    for i1 in range(n1m,n1):
        for i2 in range(n2m,n2):
            for i3 in range(n3m,n3):
                for i0 in range(len(nsys.atoms)):
                    pi0 = poss[i0]
                    x= pi0[0]/m1 +1.0/m1*i1
                    y= pi0[1]/m2 +1.0/m2*i2
                    z= pi0[2]/m3 +1.0/m3*i3
                    newsids[inc] = nsys.atoms.sid[i0]
                    newposs[inc,:] = [x,y,z]
                    newvels[inc,:] = vels[i0]
                    newfrcs[inc,:] = frcs[i0]
                    for auxname in auxnames:
                        newauxs[auxname].append(nsys.atoms.loc[i0,auxname])
                    inc += 1
    #...Use DataFrame nsys.atoms
    nsys.atoms = pd.DataFrame(columns=colnames)
    nsys.atoms[['x','y','z']] = newposs
    nsys.atoms[['vx','vy','vz']] = newvels
    nsys.atoms[['fx','fy','fz']] = newfrcs
    nsys.atoms['sid'] = newsids
    for auxname in auxnames:
        nsys.atoms[auxname] = newauxs[auxname]
    return nsys

def change_cell(nsys0,X0):
    """
    Change the cell to the new one whose lattice vectors are given by,
      (anew,bnew,cnew) = (a,b,c)*X0
    where a,b,c are column vectors that consist the original cell.
    And the atoms are reduced to those within the new cell.
    X is 3x3 matrix whose components are integer.
    """
    if X0.dtype != int:
        raise TypeError('X0.dtype is wrong.')
    if X0.shape != (3,3):
        raise TypeError('X0 has wrong shape.')
    X = np.array(X0,dtype=float)
    ncp = np.zeros(3,dtype=int)
    ncp[0] = X0[0,:].max()
    ncp[1] = X0[1,:].max()
    ncp[2] = X0[2,:].max()

    nsys = replicate(nsys0,ncp[0],ncp[1],ncp[2])
    hmat0 = nsys0.get_hmat()
    hmat = np.dot(hmat0,X0)
    print(nsys)
    sposs = nsys.get_scaled_positions()
    spnews = np.array(sposs)
    nsys.set_hmat(hmat)
    #...Since hmat0 is that of extended system,
    #...X should correspond to it.
    X[0,:] /= ncp[0]
    X[1,:] /= ncp[1]
    X[2,:] /= ncp[2]
    Xi = np.linalg.inv(X)
    for i,p in enumerate(sposs):
        pnew = np.dot(Xi,p)
        for l in range(3):
            pnew[l] = nappy.util.pbc(pnew[l])
        spnews[i,:] = pnew[:]
    nsys.set_scaled_positions(spnews)
    return nsys

def join(nsys1,nsys2,axis):
    """Join two systems into one along the given axis (0,1,or 2).
The size of the interface area is determined from the cell size of system 1 perpendicular to the axis.
It is assumed that the lattice vectors of both systems perpendicular to the given axis are very similar.

Input
-----
nsys1,nsys2 : NAPSystem
    Systems to be joined.
axis : int (0, 1 or 2)
    The axis of two systems which two systems are joined along.

Output
------
newsys1,newsys2 : NAPSystem
    The systems of nsys1 and nsys2 with vacuum region in given axis.
newsys : NAPSystem
    The system composed of nsys1 and nsys2 joined along the given axis.
    """
    if type(axis) is not int:
        raise ValueError('AXIS must be an integer.')
    if axis < 0 or axis > 2:
        raise ValueError('AXIS must be either 0, 1, or 2.')
    len1 = [0.0, 0.0, 0.0]
    len2 = [0.0, 0.0, 0.0]
    len1[0],len1[1],len1[2] = nsys1.get_lattice_lengths()
    len2[0],len2[1],len2[2] = nsys2.get_lattice_lengths()
    #...Add vacuum of length of the joining system.
    newsys1 = copy.copy(nsys1)
    newsys2 = copy.copy(nsys2)
    newsys1.add_vacuum(axis,len2[axis])
    newsys2.add_vacuum(axis,len1[axis])
    #...specorder is an union of the two systems
    specorder1 = newsys1.specorder
    specorder2 = newsys2.specorder
    specorder = copy.copy(specorder1)
    for s in specorder2:
        if s not in specorder:
            specorder.append(s)
    #...Shift positions of system2 on top of system1
    newsys = copy.copy(newsys1)
    spos1 = newsys1.get_scaled_positions()
    slen1 = spos1[:,axis].max() -spos1[:,axis].min()
    spos2 = newsys2.get_scaled_positions()
    slen2 = spos2[:,axis].max() -spos2[:,axis].min()
    syms2 = newsys2.get_symbols()
    shift2 = slen1 +(1.0 -slen1 -slen2)/2
    for i in range(len(spos2)):
        spos2[i,axis] += shift2
    #...Finaly, add them to the NEWSYS
    newsys.add_atoms(syms2,spos2)
    
    return newsys1,newsys2,newsys
    

if __name__ == "__main__":

    print(__doc__)
