#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
Read XDATCAR and write POSCAR_##### files.

Usage:
  xdatcar.py [options] XDATCAR

Options:
  -h, --help  Show this message and exit.

"""
from __future__ import print_function

import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem
from nappy.atom import Atom

__author__ = "Ryo KOBAYASHI"
__version__ = "180409"

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class XDATCAR(object):
    """
    A class for read/write XDATCAR file.
    """
    
    def __init__(self,fname="XDATCAR"):
        self.nsyss = []
        self.read(fname=fname)

    def read(self,fname='XDATCAR'):
        f= open(fname,'r')
        while True:
            #....."Direct configuration= #"
            try:
                line = f.readline()
                #...Check configuration
                if not "Direct configuration=" in line:
                    alc,a1,a2,a3,species,num_atoms = self.read_header(f)
                    # print(species)
                    # print(num_atoms)
                    natm = 0
                    for n in num_atoms:
                        natm += n
                    continue
                elif len(line.split()) == 0:
                    break
            except:
                break
            #...Following lines: atom positions
            #...Make a NAPSystem from these information and register to nsyss
            nsys = NAPSystem(specorder=species)
            nsys.set_lattice(alc,a1,a2,a3)
            inc = 0
            for inum,ni in enumerate(num_atoms):
                for j in range(ni):
                    inc += 1
                    ai = Atom()
                    ai.set_id(inc)
                    ai.set_symbol(species[inum])
                    data= f.readline().split()
                    # print(data)
                    if not is_number(data[0]):
                        raise ValueError('Wrong format?')
                    x1,x2,x3 = [ float(x) for x in data[0:3] ]
                    ai.set_pos(x1,x2,x3)
                    ai.set_vel(0.0,0.0,0.0)
                    nsys.add_atom(ai)
            self.nsyss.append(nsys)
        f.close()
        return

    def read_header(self,f):
        #.....2nd line: multiplying factor
        alc= float(f.readline().split()[0])
        #.....3rd-5th lines: lattice vectors
        a1 = np.array([ float(x) for x in f.readline().split() ])
        a2 = np.array([ float(x) for x in f.readline().split() ])
        a3 = np.array([ float(x) for x in f.readline().split() ])
        #.....6th line: num of atoms
        data= f.readline().split()
        species = []
        if not data[0].isdigit(): # if it is not digit, read next line
            species = data
            data = f.readline().split()
        num_atoms= np.array([ int(n) for n in data ])
        return alc,a1,a2,a3,species,num_atoms

if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['XDATCAR']
    xdatcar = XDATCAR(fname=fname)

    nsyss = xdatcar.nsyss
    
    for i,nsys in enumerate(nsyss):
        foutname = 'POSCAR_{0:05d}'.format(i)
        print(' >>> '+foutname)
        nsys.write_POSCAR(fname=foutname)
        
    print('Done')
