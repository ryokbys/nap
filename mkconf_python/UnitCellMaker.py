import numpy as np

from AtomSystem import AtomSystem
from Atom import Atom

def bccBravaisCell():
    s= AtomSystem()
    #...lattice
    alc= 1.0
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(alc,a1,a2,a3)
    #...atom1
    atom= Atom()
    atom.set_pos( 0.0, 0.0, 0.0 )
    atom.set_sid(1)
    s.add_atom(atom)
    #...atom2
    atom= Atom()
    atom.set_pos( 0.5, 0.5, 0.5 )
    atom.set_sid(1)
    s.add_atom(atom)
    return s

def fccBravaisCell():
    s= AtomSystem()
    #...lattice
    alc= 1.0
    a1= np.array([ 1.0, 0.0, 0.0 ])
    a2= np.array([ 0.0, 1.0, 0.0 ])
    a3= np.array([ 0.0, 0.0, 1.0 ])
    s.set_lattice(alc,a1,a2,a3)
    #...atom1
    atom= Atom()
    atom.set_pos( 0.0, 0.0, 0.0 )
    atom.set_sid(1)
    s.add_atom(atom)
    #...atom2
    atom= Atom()
    atom.set_pos( 0.5, 0.5, 0.0 )
    atom.set_sid(1)
    s.add_atom(atom)
    #...atom3
    atom= Atom()
    atom.set_pos( 0.5, 0.0, 0.5 )
    atom.set_sid(1)
    s.add_atom(atom)
    #...atom4
    atom= Atom()
    atom.set_pos( 0.0, 0.5, 0.5 )
    atom.set_sid(1)
    s.add_atom(atom)
    return s


