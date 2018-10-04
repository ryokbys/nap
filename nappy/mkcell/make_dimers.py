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
  --specorder SPECORDER
              Species. [default: H]
"""
from __future__ import print_function

from docopt import docopt
import numpy as np

from nappy.napsys import NAPSystem
from nappy.atom import Atom

__author__ = "RYO KOBAYASHI"
__version__ = "181004"

def make_dimer(distance,latconst,specorder):
    if distance > latconst/2:
        raise ValueError('Lattice size is too small for given distance.')
    s = NAPSystem(specorder=specorder)
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    positions=[(0.0, 0.0, 0.0),
               (distance/latconst, 0.0, 0.0)]
    for p in positions:
        atom = Atom()
        atom.set_pos(p[0],p[1],p[2])
        atom.set_symbol(specorder[0])
        s.add_atom(atom)
    return s

if __name__ == "__main__":

    args = docopt(__doc__)
    n = int(args['N'])
    dmin = float(args['DMIN'])
    dmax = float(args['DMAX'])
    vac = float(args['--vacuum'])
    ofname= args['-o']
    specorder = args['--specorder'].split(',')

    if n <= 1:
        raise ValueError('N must be larger than 1.')
    dd = (dmax-dmin)/(n-1)
    latconst = max(dmax+vac,2.0*dmax)
    print('Lattice cconstant = {0:7.3f}'.format(latconst))
    for i in range(n):
        d = dmin + dd*(i+1)
        struct = make_dimer(d,latconst,specorder)
        fname = ofname+'_{0:03d}'.format(i+1)
        struct.write(fname)
        print('{0:s}:  distance = {1:6.3f}'.format(fname,d))
