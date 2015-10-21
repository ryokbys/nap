#!/usr/bin/env python
"""
Make typical crystalline structures in Bravais lattice.

Usage:
  bravais_maker.py (sc|bcc|fcc|hcp|diamond) [options]

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

from pmdsys import PMDSystem
from atom import Atom

def make_sc(latconst=1.0):
    """
    Make a Bravais cell of simple cubic structure.
    """
    s= PMDSystem()
    #...lattice
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(latconst,a1,a2,a3)
    p(0.00, 0.00, 0.00)
    atom= Atom()
    atom.set_pos(p[0],p[1],p[2])
    atom.set_sid(1)
    s.add_atom(atom)
    return s

def make_bcc(latconst=1.0):
    """
    Make a Bravais cell of bcc structure.
    """
    s= PMDSystem()
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
        atom.set_sid(1)
        s.add_atom(atom)
    return s

def make_fcc(latconst=1.0):
    """
    Make a Bravais cell of fcc structure.
    """
    s= PMDSystem()
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
        atom.set_sid(1)
        s.add_atom(atom)
    return s

def make_diamond(latconst=1.0):
    """
    Make a Bravais cell of diamond structure.
    """
    s= PMDSystem()
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
        atom.set_sid(1)
        s.add_atom(atom)
    return s

def make_hcp(latconst=1.0):
    """
    Make a Bravais cell of hcp structure.
    """
    s= PMDSystem()
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
        atom.set_sid(1)
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

    if struct == None:
        print "Something wrong: structure is not created..."
        sys.exit()

    struct.expand(nx,ny,nz)
    
    struct.write(ofname)
