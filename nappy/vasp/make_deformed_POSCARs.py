#!/bin/env python
# -*- coding: utf-8 -*-
"""
Make deformed POSCAR files from the non-deformed POSCAR file.

Usage:
    make_deformed_POSCARs.py isotropic [options] POSCAR
    make_deformed_POSCARs.py uniaxial [options] POSCAR
    make_deformed_POSCARs.py orthorhombic [options] POSCAR
    make_deformed_POSCARs.py monoclinic [options] POSCAR
    make_deformed_POSCARs.py random [options] POSCAR

Options:
    -h, --help  Show this help message and exit.
    -d DEV      Maximum value of each strain element. [default: 0.01]
    -n NDEV     Number of deviation in each strain element. [default: 2]
    -o OFFSET   Offset of sequential number in output files. [default: 0]
"""

import numpy as np
import copy
from docopt import docopt

from POSCAR import POSCAR

#======================================== subroutines and functions

def isotropic(poscar):
    inc= _offset
    afac0= copy.deepcopy(poscar.afac)
    
    da= 2*_dev/_ndev *afac0
    al0= afac0 *(1.0 -dev)
    for na in range(_ndev+1):
        poscar.afac= al0 +da*na
        inc += 1
        fname = _fname+'-{0:03d}'.format(inc)
        print fname
        poscar.write(fname=fname)
    # retore poscar.afac
    poscar.afac = afac0
    

def uniaxial(poscar):
    ho= copy.deepcopy(poscar.h)
    h = np.zeros((3,3),dtype=float)
    inc= _offset
    #...uniaxial deformation in x-direction for C11
    de11= 2*_dev/_ndev
    #print " de11 =",de11
    for ne11 in range(_ndev+1):
        e11= -_dev +ne11*de11
        #print "  e11 =",e11
        h[0,0]= ho[0,0]*(1.0+e11)
        h[0,1]= ho[0,1]
        h[0,2]= ho[0,2]
        h[1,0]= ho[1,0]
        h[1,1]= ho[1,1]
        h[1,2]= ho[1,2]
        h[2,0]= ho[2,0]
        h[2,1]= ho[2,1]
        h[2,2]= ho[2,2]
        poscar.h= h
        inc += 1
        fname = _fname+'-{0:03d}'.format(inc)
        print fname
        poscar.write(fname=fname)
    #...restore poscar.h
    poscar.h= ho


def orthorhombic(poscar):
    """
    Orthorhombic volume-conserving strain for (C11-C12).
    """
    ho= copy.deepcopy(poscar.h)
    h = np.zeros((3,3),dtype=float)
    inc= _offset

    dg11= _dev/_ndev
    for ng11 in range(_ndev+1):
        g11= ng11*dg11
        h[0,0]= ho[0,0] +g11
        h[0,1]= ho[0,1]
        h[0,2]= ho[0,2]
        h[1,0]= ho[1,0]
        h[1,1]= ho[1,1] -g11
        h[1,2]= ho[1,2]
        h[2,0]= ho[2,0]
        h[2,1]= ho[2,1]
        h[2,2]= ho[2,2] +g11**2 /(1.0-g11**2)
        poscar.h= h
        inc += 1
        fname = _fname+'-{0:03d}'.format(inc)
        print fname
        poscar.write(fname=fname)
    #...restore poscar.h
    poscar.h= ho


def monoclinic(poscar):
    """
    Monoclinic volume-conserving strain for C44.
    """
    ho= copy.deepcopy(poscar.h)
    h = np.zeros((3,3),dtype=float)
    inc= _offset

    dg12= _dev/_ndev
    for ng12 in range(_ndev+1):
        g12= ng12*dg12
        h[0,0]= ho[0,0]
        h[0,1]= ho[0,1] +g12/2
        h[0,2]= ho[0,2]
        h[1,0]= ho[1,0] +g12/2
        h[1,1]= ho[1,1]
        h[1,2]= ho[1,2]
        h[2,0]= ho[2,0]
        h[2,1]= ho[2,1]
        h[2,2]= ho[2,2] +g12**2 /(4.0-g12**2)
        poscar.h= h
        inc += 1
        fname = _fname+'-{0:03d}'.format(inc)
        print fname
        poscar.write(fname=fname)
    #...restore poscar.h
    poscar.h= ho


def random(poscar):
    """
    Random deformation of the simulation cell.
    """
    #...original a vectors
    ho= copy.deepcopy(poscar.h)
    h = np.zeros((3,3),dtype=float)
    
    inc= _offset
    
    de11= 2*_dev/_ndev
    de22e33= 2*_dev/_ndev
    dg12= 2*_dev/_ndev
    dg23g31= 2*_dev/_ndev
    for ne11 in range(_ndev+1):
        e11= -_dev +ne11*de11
        for ne22e33 in range(_ndev+1):
            e22e33= -_dev +ne22e33*de22e33
            for ng12 in range(_ndev+1):
                g12= -_dev +ng12*dg12
                for ng23g31 in range(_ndev+1):
                    g23g31= -_dev +ng23g31*dg23g31
                    h[0,0]= (1.0+e11)*ho[0,0] +g12         *ho[0,1] +g23g31      *ho[0,2]
                    h[0,1]= g12      *ho[0,0] +(1.0+e22e33)*ho[0,1] +g23g31      *ho[0,2]
                    h[0,2]= g23g31   *ho[0,0] +g23g31      *ho[0,1] +(1.0+e22e33)*ho[0,2]
                    h[1,0]= (1.0+e11)*ho[1,0] +g12         *ho[1,1] +g23g31      *ho[1,2]
                    h[1,1]= g12      *ho[1,0] +(1.0+e22e33)*ho[1,1] +g23g31      *ho[1,2]
                    h[1,2]= g23g31   *ho[1,0] +g23g31      *ho[1,1] +(1.0+e22e33)*ho[1,2]
                    h[2,0]= (1.0+e11)*ho[2,0] +g12         *ho[2,1] +g23g31      *ho[2,2]
                    h[2,1]= g12      *ho[2,0] +(1.0+e22e33)*ho[2,1] +g23g31      *ho[2,2]
                    h[2,2]= g23g31   *ho[2,0] +g23g31      *ho[2,1] +(1.0+e22e33)*ho[2,2]
                    poscar.h= h
                    inc += 1
                    fname = _fname+'-{0:03d}'.format(inc)
                    print fname
                    poscar.write(fname=fname)
    #...restore poscar.h
    poscar.h= ho

############################################################ main

if __name__ == "__main__":

    args = docopt(__doc__)
    _dev = float(args['-d'])
    _ndev= int(args['-n'])
    _offset = int(args['-o'])
    _fname= args['POSCAR']

    poscar= POSCAR()
    poscar.read(fname=_fname)
    
    if args['isotropic']:
        isotropic(poscar)
    elif args['uniaxial']:
        uniaxial(poscar)
    elif args['orthorhombic']:
        orthorhombic(poscar)
    elif args['monoclinic']:
        monoclinic(poscar)
    elif args['random']:
        random(poscar)

