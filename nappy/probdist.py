#!/usr/bin/env python
"""
Compute probability distribution of specified species.
The size of the system cell is determined as that of the system read 1st.

Usage:
  probdist.py [options] FILES [FILES...]

Options:
  -h, --help  Show this message and exit.
  -w WIDTH    Mesh width in Ang. [default: 0.2]
  --sigma SIGMA
              Sigma value of the gaussian smearing. If it is positive,
              use this value in Ang. If it is negative, adopt |SIGMA|*WIDTH. [default: -2.0]
  --sid SID   ID of species whose distributions are computed.
              If it is 0, all the atoms are to be taken into account. [default: 0]
  --specorder SPECORDER
              Species order separated by comma. [default: None]
"""
from __future__ import print_function

import os
from docopt import docopt
import numpy as np

from nappy.common import get_key
from nappy.napsys import NAPSystem

__author__ = "RYO KOBAYASHI"
__version__ = "181225"

def get_num_division(nsys,width):
    ndivs = np.zeros(3,dtype=int)
    la,lb,lc = nsys.get_lattice_lengths()
    ndivs[0] = la/width +1
    ndivs[1] = lb/width +1
    ndivs[2] = lc/width +1
    return ndivs

def get_prob_dist(ndivs,nsys,sid,sgm):
    tsgm2 = 2.0 *sgm *sgm
    tsgm2i = 1.0 /tsgm2
    pdist = np.zeros(ndivs,dtype=float)
    la,lb,lc = nsys.get_lattice_lengths()
    ra = 1.0/ndivs[0]
    rb = 1.0/ndivs[1]
    rc = 1.0/ndivs[2]
    dv = la*ra*lb*rb*lc*rc
    prefactor = dv / (np.pi*tsgm2)**1.5
    for ia in range(nsys.num_atoms()):
        if nsys.get_atom_attr(ia,'sid') != sid:
            continue
        pi = nsys.get_atom_attr(ia,'pos')
        for i in range(ndivs[0]):
            rai = ra*i -pi[0]
            rai = rai if abs(rai) < 0.5 else rai -1.0*np.sign(rai)
            da = la*rai
            if abs(da) > 2.5*sgm:
                continue
            for j in range(ndivs[1]):
                rbi = rb*j -pi[1]
                rbi = rbi if abs(rbi) < 0.5 else rbi -1.0*np.sign(rbi)
                db = lb*rbi
                if abs(db) > 2.5*sgm:
                    continue
                for k in range(ndivs[2]):
                    rci = rc*k -pi[2]
                    rci = rci if abs(rci) < 0.5 else rci -1.0*np.sign(rci)
                    dc = lc*rci
                    if abs(dc) > 2.5*sgm:
                        continue
                    dr2 = da*da +db*db +dc*dc
                    pdist[i,j,k] = pdist[i,j,k] \
                                   +prefactor*np.exp(-dr2 *tsgm2i)
    return pdist

def write_CHGCAR(nsys,pdist,fname='CHGCAR'):
    poscar = 'POSCAR_tmp'
    nsys.write_POSCAR(fname=poscar)
    with open(poscar,'r') as f:
        lines = f.readlines()
    with open(fname,'w') as f:
        for line in lines:
            f.write(line)
        f.write('\n')
        shape = np.shape(pdist)
        f.write(' {0:d} {1:d} {2:d}\n'.format(*shape))
        n = 0
        for k in range(shape[2]):
            for j in range(shape[1]):
                for i in range(shape[0]):
                    p = pdist[i,j,k]
                    f.write('  {0:12.4e}'.format(p))
                    if n % 5 == 4:
                        f.write('\n')
                    n += 1
    os.system('rm '+poscar)
    return None

def normalize_pdist(pdist,nsys,sid):
    print(' Normalize pdist...')
    s = 0.0
    s = np.sum(pdist)
    na = nsys.num_atoms(sid=sid)
    pdist *= float(na)/s
    #print('   Max in pdist = ',np.max(pdist))
    #print('   Sum of pdist = ',np.sum(pdist))
    #print('   Num of atoms considered = ',na)
    return pdist
    
if __name__ == "__main__":

    args = docopt(__doc__)
    files = args['FILES']
    sid = int(args['--sid'])
    width = float(args['-w'])
    sgm = float(args['--sigma'])
    if sgm < 0.0:
        sgm = abs(sgm)*width
    specorder = args['--specorder'].split(',')
    if 'None' in specorder and not files[0].find('POSCAR') > -1:
        raise ValueError('ERROR: specorder must be specified, unless files are POSCAR format.')

    if sid == 0:
        raise ValueError('ERROR: sid should be specified and should not be 0.')

    files.sort(key=get_key,reverse=True)

    nsys0 = NAPSystem(fname=files[0],specorder=specorder)
    a,b,c = nsys0.get_lattice_angles()
    if abs(a-np.pi/2) > 180.0/np.pi or abs(b-np.pi/2) > 180.0/np.pi \
       or abs(c-np.pi/2) > 180.0/np.pi:
        raise ValueError('ERROR: Currently only available for orthogonal lattice, a,b,c=',a,b,c)
    ndivs = get_num_division(nsys0,width)
    if np.dot(ndivs,ndivs) > 10000000:
        raise ValueError('ERROR: ndivs too large.',ndivs)
    else:
        print(' # of divisions = {0:d} {1:d} {2:d}'.format(*ndivs))

    pdist = np.zeros(ndivs,dtype=float)
    for f in files:
        print(' Reading {0:s}...'.format(f))
        nsys = NAPSystem(fname=f,specorder=specorder)
        pdist += get_prob_dist(ndivs,nsys,sid,sgm)
    pdist = normalize_pdist(pdist,nsys0,sid)
    write_CHGCAR(nsys0,pdist)
    print(' Wrote CHGCAR')
    
