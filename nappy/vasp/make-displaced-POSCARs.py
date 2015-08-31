#!/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------
# Make displaced POSCAR files from the non-displacedPOSCAR file.
#-----------------------------------------------------------------------
#

import numpy as np
import copy,optparse
from math import sin,cos,pi,sqrt
from random import random

from POSCAR import POSCAR


def get_displacement(r,theta,phi):
    dx= r*sin(theta) +r*cos(phi)
    dy= r*sin(theta) +r*sin(phi)
    dz= r*cos(theta)
    return np.array([dx,dy,dz])

############################################################ main

if __name__ == "__main__":
    usage= '%prog [options] [POSCAR]'

    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-m","--max-displacement",dest="maxdis",
                      type="float",default=0.05,
                      help="maximum displacement in Angstrom."
                      +" Default value is 0.05.")
    parser.add_option("-n","--num-output",dest="nout",
                      type="int",default=10,
                      help="number of output files in total."
                      +" Default value is 10.")
    parser.add_option("-o","--offset",dest="offset",
                      type="int",default=0,
                      help="offset of sequential number in output file.")
    (options,args)= parser.parse_args()

    maxdis= options.maxdis
    nout= options.nout
    offset= options.offset
    
    fname= args[0]
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

