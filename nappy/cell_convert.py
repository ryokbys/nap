#!/usr/bin/env python
"""
Change unit vectors from the one of the given system.

Usage:
  cell_convert.py [options] INFILE

Options:
  -h,--help  Show this message and exit.
  --specorder=SPECORDER
             Set species order. [default: Al,Mg,Si]
"""
from __future__ import print_function

import os,sys,copy
from docopt import docopt
import numpy as np

sys.path.append(__file__)
from atom import Atom
from pmdsys import PMDSystem

__author__ = 'Ryo KOBAYASHI'
__version__ = '160510'

if __name__ == '__main__':
    
    args = docopt(__doc__,version=__version__)
    infile = args['INFILE']
    specorder = args['--specorder'].split(',')

    psys = PMDSystem(fname=infile,specorder=specorder)
    psys.assign_pbc()
    psys.a1 = psys.a1 *psys.alc
    psys.a2 = psys.a2 *psys.alc
    psys.a3 = psys.a3 *psys.alc
    psys.alc = 1.0
    print('a1  = ',psys.a1)
    print('a2  = ',psys.a2)
    print('a3  = ',psys.a3)
    
    pos = psys.get_real_positions()
    spos = psys.get_scaled_positions()
    for i in range(min(len(psys.atoms),100)):
        a = psys.atoms[i]
        print('{0:5d} {1:s}'.format(a.id,a.symbol)
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(spos[i,0],
                                                       spos[i,1],
                                                       spos[i,2])
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(pos[i,0],
                                                       pos[i,1],
                                                       pos[i,2]))
    
    # print(psys.get_scaled_positions())
    # print(psys.get_real_positions())
    sa1new = np.zeros(3,dtype=float)
    sa2new = np.zeros(3,dtype=float)
    sa3new = np.zeros(3,dtype=float)
    #tmp = raw_input('Input new a1 vector: ')
    #a1new[:] = [ float(x) for x in tmp.split(',') ]
    sa1new[:] = [ 0.5, 0.5, 0.0]
    #tmp = raw_input('Input new a2 vector: ')
    #a2new[:] = [ float(x) for x in tmp.split(',') ]
    sa2new[:] = [ 0.0, 1.0, 0.0 ]
    #tmp = raw_input('Input new a3 vector: ')
    #a3new[:] = [ float(x) for x in tmp.split(',') ]
    sa3new[:] = [ 0.5, 0.5, 1.0 ]
    hmat = psys.get_hmat()
    a1new = np.dot(hmat,sa1new)
    a2new = np.dot(hmat,sa2new)
    a3new = np.dot(hmat,sa3new)
    print('new a1 in hmat_orig =',sa1new)
    print('new a2 in hmat_orig =',sa2new)
    print('new a3 in hmat_orig =',sa3new)
    print('new a1 =',a1new)
    print('new a2 =',a2new)
    print('new a3 =',a3new)
    psnew = PMDSystem(specorder=specorder)
    psnew.set_lattice(psys.alc,a1new,a2new,a3new)

    # Expand the original system for the search of atoms to be included 
    # in the new system.
    # First, compute how much we have to expand the original system
    hi = np.linalg.inv(hmat)
    icsa1new = [0,0,0]
    icsa2new = [0,0,0]
    icsa3new = [0,0,0]
    for i in range(3):
        if sa1new[i] < 0.0:
            icsa1new[i] = int(sa1new[i]-1.0)
        else:
            icsa1new[i] = int(sa1new[i]+1.0)
        if sa2new[i] < 0.0: 
            icsa2new[i] = int(sa2new[i]-1.0) 
        else:
            icsa2new[i] = int(sa2new[i]+1.0)
        if sa3new[i] < 0.0:
            icsa3new[i] = int(sa3new[i]-1.0) 
        else:
            icsa3new[i] = int(sa3new[i]+1.0)
    print(icsa1new)
    print(icsa2new)
    print(icsa3new)
    for i in range(3):
        if icsa1new[i] == 0:
            raise RuntimeError('icsa1new[i] == 0')
        if icsa2new[i] == 0:
            raise RuntimeError('icsa2new[i] == 0')
        if icsa3new[i] == 0:
            raise RuntimeError('icsa3new[i] == 0')
    irange1 = (min(icsa1new[0],icsa2new[0],icsa3new[0]),
               max(icsa1new[0],icsa2new[0],icsa3new[0]))
    irange2 = (min(icsa1new[1],icsa2new[1],icsa3new[1]),
               max(icsa1new[1],icsa2new[1],icsa3new[1]))
    irange3 = (min(icsa1new[2],icsa2new[2],icsa3new[2]),
               max(icsa1new[2],icsa2new[2],icsa3new[2]))

    print('irange1: ',irange1)
    print('irange2: ',irange2)
    print('irange3: ',irange3)
    expos = []
    symbols = psys.get_symbols()
    print('symbols :',symbols)
    exsymbols = []
    print('Expanding the original system...')
    for n3 in range(min(0,irange3[0]),irange3[1]):
        for n2 in range(min(0,irange2[0]),irange2[1]):
            for n1 in range(min(0,irange1[0]),irange1[1]):
                for ia in range(len(spos)):
                    sposi = copy.deepcopy(spos[ia])
                    sposi[0] += n1
                    sposi[1] += n2
                    sposi[2] += n3
                    posi = np.dot(hmat,sposi)
                    symbol = symbols[ia]
                    # print(ia,n1,n2,n3,symbol,sposi)
                    expos.append(posi)
                    exsymbols.append(symbol)

    print('Extracting the atoms inside the new unit vectors...')
    hmat= psnew.get_hmat()
    hi = np.linalg.inv(hmat)
    for ia,posi in enumerate(expos):
        sposi = np.dot(hi,posi)
        if 0.0 <= sposi[0] < 1.0 and \
           0.0 <= sposi[1] < 1.0 and \
           0.0 <= sposi[2] < 1.0:
            atom = Atom()
            symbol = exsymbols[ia]
            print('{0:5d} {1:s}'.format(ia,symbol)
                  +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(sposi[0],
                                                           sposi[1],
                                                           sposi[2]))
            
            atom.set_symbol(symbol)
            atom.set_pos(sposi[0],sposi[1],sposi[2])
            psnew.add_atom(atom)
            
    tmp = None
    #tmp = raw_input('Input periodic shift vector if you want: ')
    tmp = ' 0.5, 0.0, 0.5'
    if tmp:
        shift = [ float(x) for x in tmp.split(',')]
        for a in psnew.atoms:
            a.pos[0] += shift[0]
            a.pos[1] += shift[1]
            a.pos[2] += shift[2]
        psnew.assign_pbc()
    psnew.write_POSCAR(infile+'.new')
    print('Check '+infile+'.new')
