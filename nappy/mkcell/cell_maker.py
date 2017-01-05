#!/usr/bin/env python
"""
Make typical crystalline structures of conventional cell.

Usage:
  cell_maker.py (sc|bcc|fcc|hcp|diamond) [options]

Options:
  -h, --help  Show this help message and exit.
  --nx=NX    Num of copies of a unit cell in x-direction. [default: 1]
  --ny=NY    Num of copies of a unit cell in y-direction. [default: 1]
  --nz=NZ    Num of copies of a unit cell in z-direction. [default: 1]
  -l, --lattice-constant=LATCONST
             Lattice constant of an axis. [default: 5.472]
  -o OUTFILE
             Output file name. Format is detected automatically. [default: POSCAR]
"""

import os,sys
import numpy as np
from docopt import docopt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

from napsys import NAPSystem
from atom import Atom

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
    p=[0.00, 0.00, 0.00]
    atom= Atom()
    atom.set_pos(p[0],p[1],p[2])
    atom.set_symbol(_default_specorder[0])
    s.add_atom(atom)
    return s


def make_bcc(latconst=1.0):
    """
    Make a cell of bcc structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    positions=[(0.00, 0.00, 0.00),
               (0.50, 0.50, 0.50)]
    for p in positions:
        atom= Atom()
        atom.set_pos(p[0],p[1],p[2])
        atom.set_symbol(_default_specorder[0])
        s.add_atom(atom)
    return s


def make_fcc(latconst=1.0,size=(1,1,1)):
    """
    Make a cell of fcc structure.
    """
    s= NAPSystem(specorder=_default_specorder)
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    positions=[(0.00, 0.00, 0.00),
               (0.50, 0.50, 0.00),
               (0.50, 0.00, 0.50),
               (0.00, 0.50, 0.50)]
    for p in positions:
        atom= Atom()
        atom.set_pos(p[0],p[1],p[2])
        atom.set_symbol(_default_specorder[0])
        s.add_atom(atom)
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
    positions=[(0.00, 0.00, 0.00),
               (0.50, 0.50, 0.00),
               (0.50, 0.00, 0.50),
               (0.00, 0.50, 0.50),
               (0.25, 0.25, 0.25),
               (0.75, 0.75, 0.25),
               (0.75, 0.25, 0.75),
               (0.25, 0.75, 0.75)]
    for p in positions:
        atom= Atom()
        atom.set_pos(p[0],p[1],p[2])
        atom.set_symbol(_default_specorder[0])
        s.add_atom(atom)
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
    positions=[(0.00, 0.00, 0.00),
               (1.0/3, 2.0/3, 0.50)]
    for p in positions:
        atom= Atom()
        atom.set_pos(p[0],p[1],p[2])
        atom.set_symbol(_default_specorder[0])
        s.add_atom(atom)
    return s


#=======================================================================
if __name__ == "__main__":

    args= docopt(__doc__)
    #print args

    nx= int(args['--nx'])
    ny= int(args['--ny'])
    nz= int(args['--nz'])
    latconst= float(args['--lattice-constant'])
    ofname= args['-o']

    struct= None
    if args['sc']:
        struct= make_sc(latconst)
    elif args['bcc']:
        struct= make_bcc(latconst)
    elif args['fcc']:
        struct= make_fcc(latconst)
    elif args['hcp']:
        struct= make_hcp(latconst)
    elif args['diamond']:
        struct= make_diamond(latconst)

    if struct is None:
        print "Something wrong: structure is not created..."
        sys.exit()

    struct.repeat(nx,ny,nz)
    
    struct.write(ofname)
