#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make deformed POSCAR files from the non-deformed POSCAR file.

Usage:
  make_deformed_POSCARs.py isotropic [options] POSCAR
  make_deformed_POSCARs.py uniaxial [options] POSCAR
  make_deformed_POSCARs.py shear [options] POSCAR

Options:
  -h, --help  Show this help message and exit.
  --fmax=FMAX 
              Maximum factor to be multiplied to the cell. [default: 2.0]
  --fmin=FMIN 
              Minimum factor to be multiplied to the cell. [default: 0.8]
  --dmin=DMIN 
              Minimum deviation in case of non-uniform interval. [default: 0.01]
  -n NDIV     Number of divisions in each strain element. [default: 20]
  --offset=OFFSET
              Offset of sequential number in output file. [default: 0]
  --uniform
              Use uniform interval deviation. [default: False]
  --specorder=SPECORDER
              Order of species. [default: Al,Mg,Si]
"""
from __future__ import print_function

import os,sys
import numpy as np
import copy
from docopt import docopt

sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/../')
from poscar import POSCAR
from pmdsys import PMDSystem

_length_digit = 3
_angle_digit = 3

#======================================== subroutines and functions

def get_nonuniform_factors(min_factor=0.8,max_factor=1.5,
                            delta_min=0.01,ndiv=20):
    """
    Create systems whose sizes are changed.
    Deviation changes from small to large w.r.t. 
    the distance from factor 1.0.
    
    min/max_factor
        Min/max limit of factor to be multiplied to the cell
    delta_min
        Minimum deviation.
    ndiv
        Number of division towards further direction.
    """
    
    if max_factor < min_factor:
        raise RuntimeError('Error: max_factor < min_factor !')
    
    maxdev= max(abs(1.0-min_factor),abs(max_factor-1.0))
    w= (maxdev-ndiv*delta_min)*6.0/ndiv/(ndiv+1)/(ndiv-1)
    factors=[1.0,]
    for i in range(1,ndiv+1):
        dlt= i*delta_min +w/2*i*(i+1)*(i-1)/3
        if min_factor > 1.0:
            factors.append(1.0 +dlt)
        elif max_factor < 1.0:
            factors.append(1.0 -dlt)
        else:
            if dlt < abs(max_factor-1.0)+delta_min:
                factors.append(1.0 +dlt)
            if dlt < abs(min_factor-1.0)+delta_min:
                factors.append(1.0 -dlt)
    factors.sort()
    return factors

def get_uniform_factors(min_factor=0.8, max_factor=1.5,ndiv=20):

    if max_factor < min_factor:
        raise RuntimeError('Error: max_factor < min_factor !')

    if ndiv < 2:
        raise ValueError('ndiv should be larger than 1.')
    dlt = (max_factor -min_factor)/(ndiv-1)
    factors = []
    for i in range(ndiv):
        factors.append(min_factor +dlt*i)
    return factors


def isotropic(psys0,factors=[],offset=0):
    psyslist = []
    for fac in factors:
        psys = copy.deepcopy(psys0)
        alc0 = psys0.alc
        alc = alc0 *fac
        psys.alc = alc
        psyslist.append(psys)
    return psyslist


def uniaxial(psys0,factors=[]):
    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)
    psyslist = []
    for fac in factors:
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[0,:] *= fac
        psys.set_hmat(hmat)
        psyslist.append(psys)

    if rb != ra:
        for fac in factors:
            psys = copy.deepcopy(psys0)
            hmat = psys0.get_hmat()
            hmat[1,:] *= fac
            psys.set_hmat(hmat)
            psyslist.append(psys)
    if rc != ra and rc != rb:
        for fac in factors:
            psys = copy.deepcopy(psys0)
            hmat = psys0.get_hmat()
            hmat[2,:] *= fac
            psys.set_hmat(hmat)
            psyslist.append(psys)
    return psyslist


def shear(psys0,factors=[]):
    aa,ba,ca = psys0.get_lattice_angles()
    raa = round(aa,_angle_digit)
    rba = round(ba,_angle_digit)
    rca = round(ca,_angle_digit)

    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)

    psyslist = []
    for fac in factors:
        dlt = (fac-1.0)*(a+b)/2/2
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[0,1] += dlt
        hmat[1,0] += dlt
        psys.set_hmat(hmat)
        psyslist.append(psys)

    if not (rba == rca and set((ra,rc)) == set((ra,rb))):
        for fac in factors:
            dlt = (fac-1.0)*(a+c)/2/2
            psys = copy.deepcopy(psys0)
            hmat = psys0.get_hmat()
            hmat[0,2] += dlt
            hmat[2,0] += dlt
            psys.set_hmat(hmat)
            psyslist.append(psys)
        
    if not (raa == rca and set((rb,rc)) == set((ra,rb))) and \
       not (raa == rba and set((rb,rc)) == set((rb,rc))) :
        for fac in factors:
            dlt = (fac-1.0)*(b+c)/2/2
            psys = copy.deepcopy(psys0)
            hmat = psys0.get_hmat()
            hmat[1,2] += dlt
            hmat[2,1] += dlt
            psys.set_hmat(hmat)
            psyslist.append(psys)

    return psyslist


############################################################ main

if __name__ == "__main__":

    args = docopt(__doc__)

    fmax = float(args['--fmax'])
    fmin = float(args['--fmin'])
    dmin = float(args['--dmin'])
    ndiv = int(args['-n'])
    offset = int(args['--offset'])
    uniform = args['--uniform']
    specorder = args['--specorder'].split(',')
    fname= args['POSCAR']

    psys0 = PMDSystem(fname,specorder=specorder)
    factors = []
    if uniform:
        factors = get_uniform_factors(min_factor=fmin,max_factor=fmax,
                                      ndiv=ndiv)
    else:
        factors = get_nonuniform_factors(min_factor=fmin,max_factor=fmax,
                                         delta_min=dmin,ndiv=ndiv)

    
    if args['isotropic']:
        psyslist = isotropic(psys0,factors)
    elif args['uniaxial']:
        psyslist = uniaxial(psys0,factors)
    elif args['shear']:
        psyslist = shear(psys0,factors)

    num = offset
    for p in psyslist:
        num += 1
        outfname = fname+'-{0:05d}'.format(num)
        p.write_POSCAR(outfname)


