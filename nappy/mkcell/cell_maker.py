#!/usr/bin/env python
"""
Make typical crystalline structures of conventional cell.

Usage:
  cell_maker.py (sc|bcc|bcc110|bcc111|fcc|fcc110|hcp|diamond|nacl) [options]

Options:
  -h, --help  Show this help message and exit.
  -s, --size=SIZE
              Number of copies of the unit cell, in comma-separated format, nx,ny,nz. [default: 1,1,1]
  -l, --lattice-constant=LATCONST
             Lattice constant of an axis. [default: 5.427]
  -o OUTFILE
             Output file name. Format is detected automatically. [default: POSCAR]
"""
from __future__ import print_function

import sys
import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem

_default_specorder=['Si']


def make_sc(latconst=1.0):
    """
    Make a cell of simple cubic structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    symbol = _default_specorder[0]
    symbols = [ symbol ]
    poss = [[0.00, 0.00, 0.00]]
    vels = [[0., 0., 0.]]
    frcs = [[0., 0., 0.]]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_bcc(latconst=1.0,specorder=_default_specorder):
    """
    Make a cell of bcc structure with z along [001].
    """
    s= NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.50, 0.50, 0.50]]
    symbol = _default_specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s

def make_bcc110(latconst=1.0):
    """                                                  
    Make a cell of bcc structure with z along [110].
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.414, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.414 ])
    s.set_lattice(latconst,a1,a2,a3)
    symbol = _default_specorder[0]
    symbols = [ symbol, symbol, symbol, symbol]
    poss = [[0.00, 0.00, 0.00],
            [0.00, 0.50, 0.50],
            [0.50, 0.50, 0.00],
            [0.50, 0.00, 0.50]]
    vels = [ [0., 0., 0.] for i in range(4) ]
    frcs = [ [0., 0., 0.] for i in range(4) ]
    s.add_atoms(symbols, poss, vels, frcs)
    return s

def make_bcc111(latconst=1.0):
    """
    Make a cell of bcc structure with z along [111].
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.414, 0.0, 0.0 ])
    a2= np.array([ 0.0, 2.449, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.732 ])
    s.set_lattice(latconst,a1,a2,a3)
    symbol = _default_specorder[0]
    poss=[[0.00, 0.00, 0.00],
          [0.00, 0.00, 0.50],
          [0.00, 0.333, 0.167],
          [0.00, 0.333, 0.667],
          [0.00, 0.667, 0.333],
          [0.00, 0.667, 0.833],
          [0.50, 0.167, 0.333],
          [0.50, 0.167, 0.833],
          [0.50, 0.50, 0.00],
          [0.50, 0.50, 0.50],
          [0.50, 0.833, 0.167],
          [0.50, 0.833, 0.667]]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s

def make_fcc(latconst=1.0,specorder=_default_specorder):
    """
    Make a cell of fcc structure.
    """
    s= NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.50, 0.50, 0.00],
            [0.50, 0.00, 0.50],
            [0.00, 0.50, 0.50]]
    symbol = specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_fcc110(latconst=1.0,specorder=_default_specorder):
    """
    Make a cell of fcc structure with z along [110].
    """
    s= NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.414, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.414 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.00, 0.50, 0.00],
            [0.00, 0.00, 0.50],
            [0.00, 0.50, 0.50],
            [0.50, 0.25, 0.25],
            [0.50, 0.25, 0.75],
            [0.50, 0.75, 0.25],
            [0.50, 0.75, 0.75]]
    symbol = specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s

def make_honeycomb(latconst=1.0):
    """
    Make a cell of 2D honeycomb structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.5, 0.0 ])
    a3= np.array([ 0.0, 0.0, np.sqrt(3.0) ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.50, 0.00],
            [0.50, 0.50, 1./6],
            [0.50, 0.50, 0.50],
            [0.00, 0.50, 0.5 +1.0/6] ]
    symbol = _default_specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_diamond(latconst=1.0):
    """
    Make a cell of diamond structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.50, 0.50, 0.00],
            [0.50, 0.00, 0.50],
            [0.00, 0.50, 0.50],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75]]
    symbol = _default_specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_hcp(latconst=1.0):
    """
    Make a cell of hcp structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([-0.5, np.sqrt(3.0)/2, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.633 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [1.0/3, 2.0/3, 0.50] ]
    symbol = _default_specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s

def make_graphene(latconst=2.467,size=(1,1,1)):
    """
    Make graphene.
    """
    napsys = make_honeycomb(latconst=latconst)
    napsys.repeat(*size)
    napsys.add_vacuum(2.*latconst, 0.0, 10.*latconst*np.sqrt(3))
    
    return napsys

def make_2D_triangle(latconst=3.8,size=(1,1,1)):
    """
    Make 2D triangle lattice on x-z plane. 
    Note that it is not x-y plane.
    """
    specorder = ['Ar']
    s = NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0,  0.0, 0.0 ])
    a2= np.array([ 0.0, 10.0, 0.0 ])
    a3= np.array([ 0.0,  0.0, np.sqrt(3.0) ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.50, 0.00],
            [0.50, 0.50, 0.50]]
    symbol = _default_specorder[0]
    symbols = [ symbol for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    
    s.repeat(*size)
    s.add_vacuum(2.*latconst, 0.0, 10.*latconst*np.sqrt(3))
    return s

def make_nacl(latconst=1.0):
    specorder = ['Na','Cl']
    s = NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.50, 0.00, 0.00],
            [0.00, 0.50, 0.00],
            [0.00, 0.00, 0.50],
            [0.50, 0.50, 0.00],
            [0.50, 0.00, 0.50],
            [0.00, 0.50, 0.50],
            [0.50, 0.50, 0.50],]
    symbols = ['Na','Cl','Cl','Cl','Na','Na','Na','Cl']
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


#=======================================================================
if __name__ == "__main__":

    args= docopt(__doc__)

    # nx= int(args['--nx'])
    # ny= int(args['--ny'])
    # nz= int(args['--nz'])
    nx,ny,nz = [ int(x) for x in args['--size'].split(',') ]
    latconst= float(args['--lattice-constant'])
    ofname= args['-o']

    struct= None
    if args['sc']:
        struct= make_sc(latconst)
    elif args['bcc']:
        struct= make_bcc(latconst)
    elif args['bcc110']:
        struct= make_bcc110(latconst)
    elif args['bcc111']:
        struct= make_bcc111(latconst)
    elif args['fcc']:
        struct= make_fcc(latconst)
    elif args['fcc110']:
        struct= make_fcc110(latconst)
    elif args['hcp']:
        struct= make_hcp(latconst)
    elif args['diamond']:
        struct= make_diamond(latconst)
    elif args['nacl']:
        struct = make_nacl(latconst)

    if struct is None:
        print("Something wrong: structure is not created...")
        sys.exit()

    struct.repeat(nx,ny,nz)
    
    struct.write(ofname)
