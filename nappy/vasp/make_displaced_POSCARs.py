#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make randomly displaced POSCAR files from the non-displaced POSCAR file.

Usage:
  make-displaced-POSCARs.py [options] POSCAR

Options:
  -h, --help  Show this help message and exit.
  --random-type=RTYPE
              Specify random number type. [default: uniform]
  -d, --displacement=<dis>
              Maximum displacemnt in Angstrom in case of uniform random number,
              and variance in case of normal random number. [default: 0.03]
  -n, --num-output=<nout>
              Number of output files in total. [default: 10]
  -o, --offset=<offset>
              Offset of the sequential number of output files. [default: 0]
  -i, --id-of-atom=<idatom>
              ID of an only atom to be displaced.
              If this is `-1`, all the atom except one will be displaced. [default: -1]
  --xmin=XMIN
              Min of range in x where atoms are to be displaced. [default: 0.0]
  --xmax=XMAX
              Max of range in x where atoms are to be displaced. [default: 1.0]
  --ymin=YMIN
              Min of range in y where atoms are to be displaced. [default: 0.0]
  --ymax=YMAX
              Max of range in y where atoms are to be displaced. [default: 1.0]
  --zmin=ZMIN
              Min of range in z where atoms are to be displaced. [default: 0.0]
  --zmax=ZMAX
              Max of range in z where atoms are to be displaced. [default: 1.0]
"""

import copy 
from docopt import docopt
import numpy as np
from math import sin,cos,pi,sqrt
from random import random,gauss

from poscar import POSCAR

from ase.io import read,write

__author__ = "RYO KOBAYASHI"
__version__ = "160727"

def get_displacement(r,theta,phi):
    """
    Returns displacements in Cartesian coordinate from the given polar coordinate.
    """
    dx= r*sin(theta) +r*cos(phi)
    dy= r*sin(theta) +r*sin(phi)
    dz= r*cos(theta)
    return np.array([dx,dy,dz])

def myrand(rtype):
    if rtype == 'uniform':
        return random()
    elif rtype == 'gauss':
        return gauss(0.0,1.0)
    else:
        raise RuntimeError('No such rtype: '+rtype)


def make_random_poscars(atoms0,rtype='uniform',dis=0.03,
                        nout=10,offset=0,
                        xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,
                        zmin=0.0,zmax=1.0):
    """
    Make displaced POSCARs from non-displaced POSCAR file.

    *atoms0*
        ASE's atoms object whose atoms to be displaced.
    *dis* = 0.03  (float)
        Characteristic displacement.
    *nout* = 10  (int)
        Number of displaced POSCAR files to be created.
    *offset* = 0  (int)
        Offset of the sequential POSCAR file names.
    """
    
    # poscar= POSCAR()
    # poscar.read(fname=fname)
    # #...original a vectors
    # ho= poscar.h *poscar.afac
    # hi= np.linalg.inv(ho)
    ho = atoms0.get_cell()
    hi = np.linalg.inv(ho)
    
    inc= offset
    
    for i in range(nout):
        atoms = copy.deepcopy(atoms0)
        pos = atoms.get_positions()
        spos = atoms.get_scaled_positions()
        for ia in range(1,len(pos)):
            sp = spos[ia]
            if sp[0] < xmin or sp[0] > xmax \
               or sp[1] < ymin or sp[1] > ymax \
               or sp[2] < zmin or sp[2] > zmax:
                continue
            r= abs(dis* myrand(rtype))
            theta= pi*random()
            phi= 2.0*pi*random()
            dr= get_displacement(r,theta,phi)
            #drh= np.dot(hi,dr)
            pos[ia] = pos[ia] +dr
        inc += 1
        atoms.set_positions(pos)
        atoms.wrap()
        write(fname+'-{0:03d}'.format(inc),
              images=atoms,format="vasp",vasp5=True,
              direct=True,sort=False)


def make_1disp_poscars(atoms0,rtype='uniform',dis=0.03,
                       offset=0,idatom=0,
                       xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,
                       zmin=0.0,zmax=1.0):
    """
    Make displaced POSCARs, in which only one atom is displaced,
    from non-displaced POSCAR file.

    *atoms0*
        ASE's atoms object whose atom to be displaced.
    *dis* = 0.03  (float)
        Maximum displacement.
    *offset* = 0  (int)
        Offset of the sequential POSCAR file names.
    *idatom* = 0  (int)
        ID of an atom to be displaced.
    """

    nd= 3
    
    # poscar= POSCAR()
    # poscar.read(fname=fname)
    # #...original a vectors
    # ho= poscar.h *poscar.afac
    # hi= np.linalg.inv(ho)
    ho = atoms0.get_cell()
    hi = np.linalg.inv(ho)
    
    inc= offset

    dd= np.zeros(3,dtype=float)
    d= dis/nd
    for ixyz in range(3):
        dd[:]= 0.0
        for i in range(-nd,nd+1):
            atoms= copy.deepcopy(atoms0)
            if i == 0: continue
            dd[ixyz]= i*d
            #ddh= np.dot(hi,dd)
            pos = atoms.get_positions()
            #poscar.pos[idatom]= poscar.pos[idatom] +ddh
            pos[idatom] = pos[idatom] + dd
            inc += 1
            atoms.set_positions(pos)
            atoms.wrap()
            #poscar.write(fname=fname+'-{0:03d}'.format(inc))
            write(fname+'-{0:03d}'.format(inc),
                  images=atoms,format="vasp",vasp5=True,
                  direct=True,sort=False)
            

############################################################ main

if __name__ == "__main__":

    args= docopt(__doc__)
    rtype= args['--random-type']
    dis= float(args['--displacement'])
    nout= int(args['--num-output'])
    offset= int(args['--offset'])
    idatom= int(args['--id-of-atom'])
    fname= args['POSCAR']
    xmin = float(args['--xmin'])
    xmax = float(args['--xmax'])
    ymin = float(args['--ymin'])
    ymax = float(args['--ymax'])
    zmin = float(args['--zmin'])
    zmax = float(args['--zmax'])

    try:
        atoms = read(fname,format="vasp")
    except:
        raise IOError('Cannot read '+fname+'.')

    if idatom < 0:
        make_random_poscars(atoms,rtype,dis,nout,offset,
                            xmin,xmax,ymin,ymax,zmin,zmax)
    else:
        make_1disp_poscars(fname,rtype,dis,offset,idatom,
                           xmin,xmax,ymin,ymax,zmin,zmax)
