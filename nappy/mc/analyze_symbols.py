#!/usr/bin/env python
"""
Analyze dat.symbols file written by lmc.py.

Usage:
  analyze_symbols.py [options] DATSYMBOLS

Options:
  -h, --help  Show this message and exit.
  --size=SIZE
             Size of the system. [default: 3x3x3]
  --latt-const=LATTCONST
             Lattice constant of the system. [default: 4.0448]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
from ase.lattice.cubic import FaceCenteredCubic

sys.path.append(os.path.dirname(__file__))
from mc import loads_symbols
from pairlist import make_pair_list

__author__ = "RYO KOBAYASHI"
__version__ = ""


def count_coordination(idx,spcs,nlspr,lspr,symbols):
     """
     Count number of neighbors of specis SPCS around
     an atom-IDX.
     """
     count = [ 0 for s in spcs ] 
     for jj in range(nlspr[idx]):
         ja = lspr[idx,jj]
         sj = symbols[ja]
         if sj in spcs:
             count[spcs.index(sj)] += 1
     return count

def main(args):
    fname = args['DATSYMBOLS']
    lattconst = float(args['--latt-const'])
    size = [ int(x) for x in args['--size'].split('x')]

    atoms = FaceCenteredCubic(symbol='Al',
                              latticeconstant=lattconst,
                              size=size)
    sites = atoms.get_scaled_positions()
    #...make permanent neighbor list
    nlspr,lspr = make_pair_list(atoms,rcut=3.0)

    print(' reading '+fname+'...')
    with open(fname,'r') as f:
        lines = f.readlines()
    print(' finished reading '+fname)

    out_crd = open('out.coordination','w')
    out_crd.write('# istp,  Mg-Mg,  Mg-Si, Mg-Vac,'+
                  '  Si-Mg,  Si-Si, Si-Vac,'+
                  '  Vac-Mg, Vac-Si, Vac-Vac\n')
    spcs = ['Mg','Si','Vac']
    print(' analyzing...')
    ppercent = 0
    for il,line in enumerate(lines):
        data = line.split()
        istp = int(data[0])
        counts = {'Mg':[ 0. for s in spcs ],
                  'Si':[ 0. for s in spcs ],
                  'Vac':[ 0. for s in spcs ]}
        symbols = loads_symbols(data[1])
        if len(symbols) != len(atoms):
            raise ValueError('len(symbols)!=len(atoms)')
        for ia,si in enumerate(symbols):
            if si in spcs:
                count = count_coordination(ia,spcs,
                                           nlspr,lspr,
                                           symbols)
                for i in range(len(count)):
                     counts[si][i] += count[i]
        for si in spcs:
            nsymb = symbols.count(si)
            counts[si][:] = [ float(x)/nsymb for x in counts[si] ]
        out_crd.write(' {0:6d}'.format(istp))
        for si in spcs:
            out_crd.write(' {0:7.2f}'.format(counts[si][0]))
            out_crd.write(' {0:7.2f}'.format(counts[si][1]))
            out_crd.write(' {0:7.2f}'.format(counts[si][2]))
        out_crd.write('\n')

        #...show progress in percent
        percent = int(float(il)/len(lines) *100)
        if percent != ppercent:
             sys.stdout.write('\r   {0:3d}%'.format(percent))
             ppercent = percent
    out_crd.close()
    sys.stdout.write('\n')
    

if __name__ == "__main__":

    args = docopt(__doc__)
    main(args)
