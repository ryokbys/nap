#!/usr/bin/env python
"""
Manipulation functions.

Usage:
  manipulate.py [options]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

from docopt import docopt
import numpy as np

from nappy.napsys import NAPSystem

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
    if not isinstance(nsys0,NAPSystem):
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
    if not isinstance(nsys0,NAPSystem):
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
                pj = nsys.atoms.pos[j]
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
        poss = [pi_prev]
        vels = [[0., 0., 0.]]
        frcs = [[0., 0., 0.]]
        # print('Atom added to {0:6.3f} {1:6.3f} {2:6.3f}, dist={3:.3f}'.format(*pi_prev,np.sqrt(d2min)))
        nsys.add_atoms(symbols,poss,vels,frcs)

    return nsys


if __name__ == "__main__":

    args = docopt(__doc__)

    print('manipulate.py does nothing...')