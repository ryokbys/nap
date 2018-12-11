#!/bin/env python
# -*- coding: utf-8 -*-
"""
Create POSCAR files in which some atoms are shifted so that
two stacking fault (SF) planes are formed.

usage:
    make-SF-POSCARs.py [options] POSCAR (DX DY)...

options:
    -h, --help   show this help message and exit
    -p, --plane=<xyz>
                 SF plane to be inserted [default: z]
    --zmin=<zmin>
                 minimum of the range to be shifted [default: 0.0]
    --zmax=<zmax>
                 maximum of the range to be shfited [default: 1.0]
    --ndiv=<n>   number of division in each direction [default: 10]
"""
from __future__ import print_function

__author__ = "Ryo KOBAYASHI"
__version__ = "0.1"


import numpy as np
import copy
from docopt import docopt
from math import sin,cos,pi,sqrt
from random import random

from poscar import POSCAR

def _shift_and_write(poscar0,sfplane,ndiv,zmin,zmax,dr):

    dax= {'x':0, 'y':1, 'z':2}
    iax= daxis[sfplane]
    #...flag atoms to be shifted
    shift_atom= []
    for posi in poscar0.pos:
        flag= False
        if zmin < posi[iax] < zmax:
            flag= True
        shift_atom.append(flag)

    n= 0
    for i in xrange(len(dr)):
        disp= dr[i]
        for iv in xrange(ndiv):
            towards= []
            m= 0
            for ixyz in range(3):
                towards.append(0.0)
                if ixyz != iax:
                    towards[ixyz]= disp[m]/ndiv*(iv+1)
                    m += 1
            poscar= shift(poscar0,shift_atom,towards,iax)
            n += 1
            outfname='POSCAR{0:03d}'.format(n)
            poscar.write(fname=outfname)


def shift(poscar0,shift_atom,towards,iax):
    """
    Shifts atoms indicated in *shift_atom* list in the plane perpendicular to iax
    by *towards* and returns a new poscar object.
    """

    poscar= copy.deepcopy(poscar0)
    for ia in xrange(len(shift_atom)):
        if shift_atom[ia]:
            for ixyz in xrange(3):
                if ixyz == iax:
                    continue
                poscar.pos[ia][ixyz] += towards[ixyz]
                poscar.pos[ia][ixyz]= _pbc(poscar.pos[ia][ixyz])
    return poscar

def _pbc(x):
    if x < 0.0:
        x += 1.0
    if x >= 1.0:
        x -= 1.0
    return x


############################################################ main

if __name__ == "__main__":

    args= docopt(__doc__)

    zmin= float(args['--zmin'])
    zmax= float(args['--zmax'])
    sfplane= args['--plane']
    ndiv= int(args['--ndiv'])
    _fname= args['POSCAR']
    dr= []
    for i in range(len(args['DX'])):
        dr.append((float(args['DX'][i]),float(args['DY'][i])))

    if not sfplane in ('x','y','z'):
        print(' --plane must be x, y, or z !!!')
        print(' Now --plane=',sfplane)
        exit

    if zmin < 0.0 or zmin > 1.0 or zmax < 0.0 or zmax > 1.0 or zmax < zmin:
        print(' Should be 0.0 <= zmin < zmax < 1.0 !!!')
        print(' zmin, zmax = ',zmin, zmax)
        exit

    poscar= POSCAR()
    poscar.read(fname=_fname)

    _shift_and_write(poscar,sfplane,ndiv,zmin,zmax,dr)
    
    
