#!/usr/bin/env python
"""
Qn-value analysis. Qn means that a central polyhedron connects to n neighboring polyhedra.

Usage:
  qn.py [options] FILES [FILES...]

Options:
  -h, --help  Show this message and exit.
  --rcut RC   Cutoff radius between center atom and corner atoms of a polyhedron. [default: 3.0]
  --center CENTER
              Species of the center of polyhedra. [default: None]
  --corner CORNER
              Species of the corners of polyhedra. [default: None]
  --nmax NMAX
              Maximum number of n in Qn. [default: 4]
  --out4fp
              Whether or not to write a file for fp.
"""
import os
import sys
from datetime import datetime

from docopt import docopt
import numpy as np
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "210624"

def calc_Qn(nsys,center=None,corner=None,rcut=3.0,nmax=4):
    """
    Calculate Qn distribution of the given structure whose center and corner species
    are provided.

    Input:
    ------
    nsys : NAPSystem
        The system to be analyzed.
    center, corner : string
        Names of species of center or corner atoms.

    Return:
    -------
    qn : numpy array of float.
        Probability of Qn values of the system.
    """
    if center == None or corner == None:
        raise ValueError('center and corner species must be provided.')

    #...Scan over corner atoms and count those having bridging two centers.
    symbols = nsys.get_symbols()
    bridging = [ False for i in range(len(symbols)) ]
    for ia in range(len(nsys)):
        si = symbols[ia]
        if si != corner: continue
        n = 0
        for ja,dij in nsys.neighbors_of(ia,distance=True,rcut=rcut):
            sj = symbols[ja]
            if sj != center: continue
            n += 1
        if n > 1:
            bridging[ia] = True

    #...Scan over center atoms and count corner atoms that is bridging.
    Qn = np.zeros(nmax+1,)
    for ia in range(len(nsys)):
        si = symbols[ia]
        if si != center: continue
        n = 0
        for ja,dij in nsys.neighbors_of(ia,distance=True,rcut=rcut):
            sj = symbols[ja]
            if sj != corner: continue
            if bridging[ja]:
                n += 1
        Qn[min(n,nmax)] += 1.0

    #...Normalize by the number of polyhedra, that is number of center atoms
    icenter = nsys.specorder.index(center)
    ncenter = nsys.num_atoms(sid=icenter+1)
    if ncenter == 0:
        raise ValueError('ncenter==0')
    Qn /= ncenter
    
    return Qn

if __name__ == "__main__":

    args = docopt(__doc__)
    center = args['--center']
    corner = args['--corner']
    rcut   = float(args['--rcut'])
    nmax   = int(args['--nmax'])
    fnames = args['FILES']
    out4fp = args['--out4fp']

    Qn = np.zeros(nmax+1,)
    for fname in fnames:
        nsys = nappy.io.read(fname)
        Qn += calc_Qn(nsys,center=center,corner=corner,
                      rcut=rcut,nmax=nmax)
    try:
        Qn /= len(fnames)
    except:
        raise

    print('Qn:')
    for n in range(nmax+1):
        print('  n={0:d} : {1:7.4f}'.format(n,Qn[n]))
    
    if out4fp:
        nperline = 6
        cmd = ' '.join(s for s in sys.argv)
        with open('out.qn','w') as f:
            f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('#  {0:s}\n'.format(cmd))
            f.write('#\n')
            #...Num of data, weight for the data
            f.write('  {0:6d}  {1:7.3f}\n'.format(nmax+1,1.0))
            j0 = 0
            while True:
                f.write('  '.join(' {0:7.4f}'.format(Qn[j]) for j in range(j0,j0+nperline) if j < nmax+1))
                f.write('\n')
                j0 += nperline
                if j0 >= nmax+1:
                    break
