#!/usr/bin/env python
"""
Extract structures from NEB results and write to POSCAR files.

Usage:
  nebout2poscars.py [options]

Options:
  -h, --help  Show this message and exit.
"""
import os,sys
from docopt import docopt
import numpy as np
from ase.io import write

sys.path.append(__file__)
from espresso_in import read_espresso_in

__author__ = "RYO KOBAYASHI"
__version__ = ""

def get_num_images(fname=None):
    if not fname or not os.path.exists(fname):
        raise IOError('There is no '+fname)
    num_images = 0
    with open(fname,'r') as f:
        for l in f.readlines():
            if 'num_of_images' in l:
                num_images = int(l.split()[2])
                return num_images
    raise IOError('num_of_images not found')
    

def get_images(fname,natm,num_images):
    from copy import deepcopy
    if not os.path.exists(fname):
        raise IOError('There is no '+fname)

    list_pos = []
    mode = 'idle'
    pos = np.zeros((natm,3),dtype=float)
    with open(fname) as f:
        lines = f.readlines()
        il = 0
        img = 0
        while True:
            il += 2
            img += 1
            if img > num_images:
                break
            pos[:,:] = 0.0
            for ia in range(natm):
                dat = lines[il].split()
                pos[ia,:] = [ float(x) for x in dat[1:4] ]
                il += 1
            list_pos.append(deepcopy(pos))
    return list_pos


def get_nebout_structures():
    """
    Extract structures from NEB results in the current directory.
    There must be:
    * neb.dat
    * pw_X.in, where X=(1,...,N) and N is the number of images.
    * pw.crd in which the final coordinates are written.
    """
    from ase.atoms import Atoms

    fnebdat = 'neb.dat'
    fpwin = 'pw_1.in'
    fcrd = 'pw.crd'
    
    num_images = get_num_images(fnebdat)
    natm,cell,elems,pos,cell_unit,pos_unit = read_espresso_in(fpwin)
    list_pos = get_images(fcrd,natm,num_images)
    
    list_atoms = []
    for i in range(num_images):
        atoms = Atoms(symbols=elems,
                      cell=cell,
                      scaled_positions=list_pos[i],
                      pbc=[1,1,1])
        list_atoms.append(atoms)
    return list_atoms


if __name__ == "__main__":

    args = docopt(__doc__)

    lsatoms = get_nebout_structures()
    for i,atoms in enumerate(lsatoms):
        write('POSCAR-{0:d}'.format(i+1),
              atoms,
              format='vasp',direct=True,
              sort=False,vasp5=True)
