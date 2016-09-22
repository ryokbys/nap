#!/usr/bin/env python
"""
Convert Quantum Espresso output to fitpot data.

Usage:
  espresso_in.py [options] FILE

Options:
  -h, --help  Show this message and help.
  --specorder=SPECORDER
              Specify the order of species needed to convert to pos. [default: Al,Mg,Si]

"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

sys.path.append(os.environ['HOME']+'/src/nap/nappy')
from atom import Atom
from pmdsys import PMDSystem,unitvec_to_hi,cartessian_to_scaled

__author__ = "Ryo KOBAYASHI"
__version__ = "160629"


def read_espresso_in(fname):
    """
    Read cell info and atom coordinates from the Quantum Espresso input file.

    fname: str
        File name of the input file to be read.
    """

    try:
        f= open(fname,'r')
        lines = f.readlines()
        f.close()
    except:
        raise IOError('Could not open the file: '+fname)

    for il in range(len(lines)):
        if 'ATOMIC_POSITIONS' in lines[il]:
            il_pos = il
        elif 'CELL_PARAMETERS' in lines[il]:
            il_cell = il
        elif '&SYSTEM' in lines[il]:
            il_system = il

    # read number of atoms
    il = il_system
    while True:
        il += 1
        name = lines[il].split()[0]
        if name == 'nat':
            natm = int(lines[il].split()[-1])
            break
        elif name == '/':
            raise IOError('There is no info about number of atoms.')

    # read cell parameters
    cell_unit = lines[il_cell].split()[1]
    il = il_cell
    cell = np.zeros((3,3),dtype=float)
    ixyz = 0
    while True:
        il += 1
        l = lines[il].split()
        if len(l) == 0:
            break
        cell[:,ixyz] = [ float(x) for x in l ]
        ixyz += 1
        if ixyz > 2:
            break

    # read atom species and positions
    pos_unit = lines[il_pos].split()[1]
    il = il_pos + 1
    elems = []
    pos = np.zeros((natm,3),dtype=float)
    for ia in range(natm):
        l = lines[il+ia].split()
        elems.append(l[0])
        pos[ia,:] = [ float(x) for x in l[1:4] ]
    
    return natm,cell,elems,pos,cell_unit,pos_unit
    

if __name__ == '__main__':
    args = docopt(__doc__)
    infname = args['FILE']
    specorder= args['--specorder'].split(',')

    natm,cell,elems,pos,cell_unit,pos_unit = read_espresso_in(infname)
    print('natm = ',natm)
    print('cell:',cell)
    # print('elems:',elems)
    # print('pos:',pos)

