#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make randomly displaced POSCAR files from the non-displaced POSCAR file.

Usage:
    make-displaced-POSCARs.py [options] POSCAR

Options:
    -h, --help  Show this help message and exit.
    -m, --max-displacement=<maxdis>
                Maximum displacemnt in Angstrom. [default: 0.03]
    -n, --num-output=<nout>
                Number of output files in total. [default: 10]
    -o, --offset=<offset>
                Offset of the sequential number of output files. [default: 0]
    -i, --id-of-atom=<idatom>
                ID of an only atom to be displaced.
                If this is `-1`, all the atom except one will be displaced. [default: -1]
"""

import copy
from docopt import docopt
import numpy as np
from math import sin,cos,pi,sqrt
from random import random

from POSCAR import POSCAR

def get_displacement(r,theta,phi):
    """
    Returns displacements in Cartesian coordinate from the given polar coordinate.
    """
    dx= r*sin(theta) +r*cos(phi)
    dy= r*sin(theta) +r*sin(phi)
    dz= r*cos(theta)
    return np.array([dx,dy,dz])

def make_random_poscars(fname='POSCAR',maxdis=0.05,
                      nout=10,offset=0):
    """
    Make displaced POSCARs from non-displaced POSCAR file.

    *fname* = 'POSCAR'  (str)
        Name of the non-displaced original POSCAR file.
    *maxdis* = 0.05  (float)
        Maximum displacement.
    *nout* = 10  (int)
        Number of displaced POSCAR files to be created.
    *offset* = 0  (int)
        Offset of the sequential POSCAR file names.
    """
    
    poscar= POSCAR()
    poscar.read(fname=fname)
    
    #...original a vectors
    ho= poscar.h *poscar.afac
    
    hi= np.linalg.inv(ho)
    
    inc= offset
    
    for i in range(nout):
        for ia in range(1,len(poscar.pos)):
            r= maxdis*random()
            theta= pi*random()
            phi= 2.0*pi*random()
            dr= get_displacement(r,theta,phi)
            drh= np.dot(hi,dr)
            poscar.pos[ia] = poscar.pos[ia] +drh
        inc += 1
        poscar.write(fname=fname+'-{0:03d}'.format(inc))


def make_1disp_poscars(fname='POSCAR',maxdis=0.05,
                       offset=0,idatom=0):
    """
    Make displaced POSCARs, in which only one atom is displaced,
    from non-displaced POSCAR file.

    *fname* = 'POSCAR'  (str)
        Name of the non-displaced original POSCAR file.
    *maxdis* = 0.05  (float)
        Maximum displacement.
    *offset* = 0  (int)
        Offset of the sequential POSCAR file names.
    *idatom* = 0  (int)
        ID of an atom to be displaced.
    """

    nd= 3
    
    poscar= POSCAR()
    poscar.read(fname=fname)
    
    #...original a vectors
    ho= poscar.h *poscar.afac
    
    hi= np.linalg.inv(ho)
    
    inc= offset

    dd= np.zeros(3,dtype=float)
    dis= maxdis/nd
    for ixyz in range(3):
        dd[:]= 0.0
        for i in range(-nd,nd+1):
            if i == 0: continue
            dd[ixyz]= i*dis
            ddh= np.dot(hi,dd)
            poscar.pos[idatom]= poscar.pos[idatom] +ddh
            inc += 1
            poscar.write(fname=fname+'-{0:03d}'.format(inc))


############################################################ main

if __name__ == "__main__":

    args= docopt(__doc__)
    maxdis= float(args['--max-displacement'])
    nout= int(args['--num-output'])
    offset= int(args['--offset'])
    idatom= int(args['--id-of-atom'])
    fname= args['POSCAR']

    if idatom < 0:
        make_random_poscars(fname,maxdis,nout,offset)
    else:
        make_1disp_poscars(fname,maxdis,offset,idatom)
