#!/usr/bin/env python
"""
Make typical crystalline structures of conventional cell.

Usage:
  cell_maker.py (sc|bcc|fcc|hcp|dia|nacl|zb|wz) [options]

Options:
  -h, --help  Show this help message and exit.
  -s, --size=SIZE
              Number of copies of the unit cell, in comma-separated format, nx,ny,nz. [default: 1,1,1]
  -l, --lattice-constant=LATCONST
             Lattice constant of an axis. [default: 5.427]
  -o OUTFILE
             Output file name. Format is detected automatically. [default: POSCAR]
  --orientation ORIENTATION
             Orientation of z-axis in Miller index. [default: 0,0,1]
  --celltype CELLTYPE
             Conventional or primitive. [default: conventional]
  --specorder SPECORDER
             Species order. [default: None]
"""
from __future__ import print_function

import sys
import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem
import nappy

_default_specorder=['Si']


def make_sc(latconst=1.0,specorder=None):
    """
    Make a cell of simple cubic structure.
    """
    if specorder is None:
        raise ValueError('specorder must be given.')
    s= NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    symbols = [ specorder[0] ]
    poss = [[0.00, 0.00, 0.00]]
    vels = [[0., 0., 0.]]
    frcs = [[0., 0., 0.]]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_bcc(latconst=1.0,specorder=None):
    """
    Make a cell of bcc structure with z along [001].
    """
    if specorder is None:
        raise ValueError('specorder must be given.')
    s= NAPSystem(specorder=specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    poss = [[0.00, 0.00, 0.00],
            [0.50, 0.50, 0.50]]
    symbol = specorder[0]
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

def make_fcc(latconst=1.0,specorder=None):
    """
    Make a cell of fcc structure.
    """
    if specorder is None:
        specorder = ['Al']
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


def make_fcc110(latconst=1.0,specorder=None):
    """
    Make a cell of fcc structure with z along [110].
    """
    if specorder is None:
        specorder = ['Al']
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

def make_nacl(latconst=1.0,specorder=None):
    if specorder is None:
        specorder = ['Na','Cl']
    if len(specorder) < 2:
        specorder = ['Na','Cl']
        print('Since len(specorder) < 2, specorder is reset to ',specorder)
    
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

def make_zincblend(latconst=1.0,specorder=None):
    """
    Make a cell of diamond structure.
    """
    if specorder is None:
        specorder = ['Ga','N']
    if len(specorder) < 2:
        specorder = ['Ga','N']
        print('Since len(specorder) < 2, specorder is reset to ',specorder)
    
    s= NAPSystem(specorder=specorder)
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
    symbols = [ specorder[0] if i<4 else specorder[1] for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


def make_wurtzite(latconst=1.0,specorder=None,celltype='conventional'):
    """
    Make a cell of wurtzite structure.

    - celltype: conventional or primitive
    """
    if specorder is None:
        specorder = ['Ga','N']
    if len(specorder) < 2:
        specorder = ['Ga','N']
        print('Since len(specorder) < 2, specorder is reset to ',specorder)
    
    s = NAPSystem(specorder=specorder)
    if celltype[0] == 'c':
        #...conventional cell
        a1= np.array([ 1.00, 0.00, 0.00 ])
        a2= np.array([ 0.00, np.sqrt(3.0), 0.00 ])
        a3= np.array([ 0.00, 0.00, 1.633 ])
        s.set_lattice(latconst,a1,a2,a3)
        poss = [[0.00,  0.00,  0.00],
                [0.50,  0.50,  0.00],
                [0.50,  0.5/3, 0.50],
                [0.00,  0.5/3+0.5, 0.50],
                [0.50,  0.5/3,  0.125],
                [0.00,  0.5/3+0.5, 0.125],
                [0.00,  0.00,  0.625],
                [0.50,  0.50,  0.625],]
        symbols = [ specorder[0] if i<4 else specorder[1] for i in range(len(poss)) ]
    elif cenlltype[0] == 'p':
        #...primitive cell
        a1= np.array([ 1.0, 0.0, 0.0 ])
        a2= np.array([-0.5, np.sqrt(3.0)/2, 0.0 ])
        a3= np.array([ 0.0, 0.0, 1.633 ])
        s.set_lattice(latconst,a1,a2,a3)
        poss = [[0.00,  0.00,  0.00],
                [2.0/3, 1.0/3, 0.125],
                [2.0/3, 1.0/3, 0.50],
                [0.00,  0.00,  0.625],]
        symbols = [ specorder[0] if i<2 else specorder[1] for i in range(len(poss)) ]
    vels = [ [0., 0., 0.] for i in range(len(poss)) ]
    frcs = [ [0., 0., 0.] for i in range(len(poss)) ]
    s.add_atoms(symbols,poss,vels,frcs)
    return s


#=======================================================================
if __name__ == "__main__":

    args= docopt(__doc__)
    print(args)
    # nx= int(args['--nx'])
    # ny= int(args['--ny'])
    # nz= int(args['--nz'])
    nx,ny,nz = [ int(x) for x in args['--size'].split(',') ]
    latconst= float(args['--lattice-constant'])
    ofname= args['-o']
    orient = [ int(x) for x in args['--orientation'].split(',')]
    celltype = args['--celltype']
    specorder = [ x for x in args['--specorder'].split(',')]
    if specorder[0] == 'None':
        specorder = None

    struct= None
    print('specorder = ',specorder)
    if args['sc']:
        struct= make_sc(latconst,specorder)
    elif args['bcc']:
        if orient == [0,0,1]:
            struct= make_bcc(latconst,specorder)
        elif orient == [1,1,0]:
            struct= make_bcc110(latconst)
        elif orient == [1,1,1]:
            struct= make_bcc111(latconst)
        else:
            raise ValueError('The orientation is not available: ',orient)
    elif args['fcc']:
        if orient == [0,0,1]:
            struct= make_fcc(latconst)
        elif orient == [1,1,0]:
            struct= make_fcc110(latconst)
        else:
            raise ValueError('The orientation is not available: ',orient)
    elif args['hcp']:
        struct= make_hcp(latconst)
    elif args['dia']:
        struct= make_diamond(latconst)
    elif args['nacl']:
        struct = make_nacl(latconst)
    elif args['zb']:
        struct = make_zincblend(latconst,specorder=specorder)
    elif args['wz']:
        struct = make_wurtzite(latconst,celltype=celltype,specorder=specorder)

    if struct is None:
        print("Something wrong: structure is not created...")
        sys.exit()

    struct.repeat(nx,ny,nz)
    
    nappy.io.write(struct,fname=ofname)
