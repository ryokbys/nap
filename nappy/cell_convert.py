#!/usr/bin/env python
"""
Change unit vectors of the system from the original a1, a2, and a3
to ax1, ax2, and ax3, each of them is written as (x1,x2,x3) that are
multiplying factors to a1, a2, and a3, respectively.
For example, if ax1 = (x1,x2,x3), the new a1 vector is given as:
  a1_new = x1*a1 +x2*a2 +x3*a3

Usage:
  cell_convert.py [options] INFILE

Options:
  -h,--help  Show this message and exit.
  --ax1 AX1  Factors to be multiplied to a1,a2,a3, comma separated. [default: 1.0,0.0,0.0]
  --ax2 AX2  Factors to be multiplied to a1,a2,a3, comma separated. [default: 0.0,1.0,0.0]
  --ax3 AX3  Factors to be multiplied to a1,a2,a3, comma separated. [default: 0.0,0.0,1.0]
"""
from __future__ import print_function

import os,sys,copy
from docopt import docopt
import numpy as np

#sys.path.append(__file__)
import nappy
from nappy.napsys import NAPSystem

__author__ = 'Ryo KOBAYASHI'
__version__ = '221117'

def convert_by_multiply(nsys,ax1,ax2,ax3):
    """
    Convert the input nsys by using ax1,ax2,ax3.
    See the definition in __doc__.
    """
    nsys.assign_pbc()
    #...Reset alc to 1.0
    nsys.a1 = nsys.a1 *nsys.alc
    nsys.a2 = nsys.a2 *nsys.alc
    nsys.a3 = nsys.a3 *nsys.alc
    nsys.alc = 1.0
    print(' a1  = ',nsys.a1)
    print(' a2  = ',nsys.a2)
    print(' a3  = ',nsys.a3)
    
    pos = nsys.get_real_positions()
    spos = nsys.get_scaled_positions()
    symbols = nsys.get_symbols()
    for i in range(min(len(nsys.atoms),10)):
        print(' {0:5d} {1:s}'.format(i,symbols[i])
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(spos[i,0],
                                                       spos[i,1],
                                                       spos[i,2])
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(pos[i,0],
                                                       pos[i,1],
                                                       pos[i,2]))

    
    hmat = nsys.get_hmat()
    hmati= nsys.get_hmat_inv()
    a1new = ax1[0]*nsys.a1 +ax1[1]*nsys.a2 +ax1[2]*nsys.a3
    a2new = ax2[0]*nsys.a1 +ax2[1]*nsys.a2 +ax2[2]*nsys.a3
    a3new = ax3[0]*nsys.a1 +ax3[1]*nsys.a2 +ax3[2]*nsys.a3
    sa1new = np.dot(hmati,a1new)
    sa2new = np.dot(hmati,a2new)
    sa3new = np.dot(hmati,a3new)
    print(' new a1 =',a1new)
    print(' new a2 =',a2new)
    print(' new a3 =',a3new)
    psnew = NAPSystem(specorder=nsys.specorder)
    psnew.set_lattice(nsys.alc,a1new,a2new,a3new)

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
    print(' icsa1new: ',icsa1new)
    print(' icsa2new: ',icsa2new)
    print(' icsa3new: ',icsa3new)
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

    print(' irange1: ',irange1)
    print(' irange2: ',irange2)
    print(' irange3: ',irange3)
    expos = []
    symbols = nsys.get_symbols()
    print(' symbols :',symbols)
    exsymbols = []
    print(' Expanding the original system...')
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

    print(' Extracting the atoms inside the new unit vectors...')
    hmat= psnew.get_hmat()
    hi = np.linalg.inv(hmat)
    for ia,posi in enumerate(expos):
        sposi = np.dot(hi,posi)
        if 0.0 <= sposi[0] < 1.0 and \
           0.0 <= sposi[1] < 1.0 and \
           0.0 <= sposi[2] < 1.0:
            symbol = exsymbols[ia]
            print(' {0:5d} {1:s}'.format(ia,symbol)
                  +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(sposi[0],
                                                           sposi[1],
                                                           sposi[2]))
            
            # atom.set_symbol(symbol)
            # atom.set_pos(sposi[0],sposi[1],sposi[2])
            # psnew.add_atom(atom)
            psnew.add_atoms([symbol],[sposi])
            
    # tmp = None
    # #tmp = raw_input('Input periodic shift vector if you want: ')
    # tmp = ' 0.5, 0.0, 0.5'
    # if tmp:
    #     shift = [ float(x) for x in tmp.split(',')]
    #     for a in psnew.atoms:
    #         a.pos[0] += shift[0]
    #         a.pos[1] += shift[1]
    #         a.pos[2] += shift[2]
    #     psnew.assign_pbc()
    # psnew.write_POSCAR(infile+'.new')
    # print('Check '+infile+'.new')

    psnew.assign_pbc()
    return psnew


def main():
    args = docopt(__doc__,version=__version__)
    infile = args['INFILE']
    ax1 = [ float(x) for x in args['--ax1'].split(',') ]
    ax2 = [ float(x) for x in args['--ax2'].split(',') ]
    ax3 = [ float(x) for x in args['--ax3'].split(',') ]
    if len(ax1) != 3:
        raise ValueError('len(ax1) != 3')
    if len(ax2) != 3:
        raise ValueError('len(ax2) != 3')
    if len(ax3) != 3:
        raise ValueError('len(ax3) != 3')

    nsys = nappy.io.read(infile)
    newsys = convert_by_multiply(nsys,ax1,ax2,ax3)
    nappy.io.write(newsys,infile+'_new')
    return None

if __name__ == '__main__':
    
    main()
