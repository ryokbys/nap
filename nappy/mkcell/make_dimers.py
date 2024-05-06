#!/usr/bin/env python
"""
Make cells that contain dimers.

Usage:
  make_dimers.py [options] N DMIN DMAX

Options:
  -h, --help  Show this message and exit.
  --vacuum VAC
              Size of vacuum that is added to dimer. [default: 10.0]
  -o OUTFILE
              Output file name. Format is detected automatically. [default: POSCAR]
  --species SPCS
              Species. [default: H,H]
"""
from docopt import docopt
import numpy as np

from nappy.napsys import NAPSystem
from nappy.atom import Atom

__author__ = "RYO KOBAYASHI"
__version__ = "181004"

def make_dimer(distance,latconst,spcs):
    if distance > latconst/2:
        raise ValueError('Lattice size is too small for given distance.')
    s = NAPSystem(specorder=spcs)
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    pos=[(0.0, 0.0, 0.0),
         (distance/latconst, 0.0, 0.0)]
    #...Atom 1
    atom = Atom()
    p = pos[0]
    atom.set_pos(p[0],p[1],p[2])
    atom.set_symbol(spcs[0])
    s.add_atom(atom)
    #...Atom 2
    atom = Atom()
    p = pos[1]
    atom.set_pos(p[0],p[1],p[2])
    atom.set_symbol(spcs[1])
    s.add_atom(atom)
    return s

if __name__ == "__main__":

    args = docopt(__doc__)
    n = int(args['N'])
    dmin = float(args['DMIN'])
    dmax = float(args['DMAX'])
    vac = float(args['--vacuum'])
    ofname= args['-o']
    spcs = args['--species'].split(',')
    
    if len(spcs) < 1:
        raise ValueError('Species is wrong. spcs=',spcs)
    elif len(spcs) < 2:
        spcs += spcs
        
    if n <= 1:
        raise ValueError('N must be larger than 1.')
    
    dd = (dmax-dmin)/(n-1)
    latconst = max(dmax+vac,2.0*dmax)
    print('Lattice constant = {0:7.3f}'.format(latconst))
    for i in range(n):
        d = dmin + dd*i
        struct = make_dimer(d,latconst,spcs)
        fname = ofname+'_{0:03d}'.format(i+1)
        struct.write(fname)
        print('{0:s}:  distance = {1:6.3f}'.format(fname,d))
