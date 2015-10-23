#!/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------
# Make deformed POSCAR files from the non-deformed POSCAR file.
# Deformation is performed by multiplying deformation matrix to h-mat.
#-----------------------------------------------------------------------
#

import numpy as np
import copy,optparse

from POSCAR import POSCAR

#======================================== subroutines and functions

def deform_cubic(poscar):
    """
    Deformation of the simulation cell that has cubic symmetry.
    Only uniaxial deformation for C11 and C12,
    off diagonal deformation for C44 are considered.
    Be sure that you are applying this routine to the cubic cell.
    """
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
        poscar.write(fname=_fname+'-{0:03d}'.format(inc))

    #...orthorhombic volume-conserving strain for (C11-C12)
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
        poscar.write(fname=_fname+'-{0:03d}'.format(inc))

    #...monoclinic volume-conserving strain for C44
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
        poscar.write(fname=_fname+'-{0:03d}'.format(inc))
        
    #...restore poscar.h
    poscar.h= ho


def deform_random(poscar):
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
                    poscar.write(fname=_fname+'-{0:03d}'.format(inc))
    #...restore poscar.h
    poscar.h= ho

############################################################ main

if __name__ == "__main__":
    _usage= '%prog [options] [POSCAR]'
    parser= optparse.OptionParser(usage=_usage)
    parser.add_option("-d","--dev",dest="dev",type="float",default=0.01,
                      help="maximum value of each strain element. Default: 0.01.")
    parser.add_option("-n","--num-dev",dest="ndev",type="int",default=2,
                      help="number of devision in each strain element. Default: 2.")
    parser.add_option("-o","--offset",dest="offset",type="int",default=0,
                      help="offset of sequential number in output file."
                      +" Default: 0.")
    parser.add_option("-m","--mode",dest="mode",
                      type="string",default="random",
                      help="deformation mode setting."
                      +" random: random deformation of the cell,"
                      +" cubic: for the case of cubic symmetry. Default: random.")
    
    (options,args)= parser.parse_args()
    
    _dev= options.dev
    _ndev= options.ndev
    _offset= options.offset
    _mode= options.mode
    
    _fname= args[0]
    poscar= POSCAR()
    poscar.read(fname=_fname)
    
    if _mode in ("random","Random","RANDOM"):
        deform_random(poscar)
    elif _mode in ("cubic","Cubic","CUBIC"):
        deform_cubic(poscar)

