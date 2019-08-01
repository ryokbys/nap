#!/usr/bin/env python
"""
Deform the original system by adding some strain.

Usage:
  deform.py isotropic [options] INFILE OUTFILE
  deform.py uniaxial  [options] INFILE OUTFILE
  deform.py shear     [options] INFILE OUTFILE

Options:
  -h, --help  Show this help message and exit.
  --in-format=INFMT
              Format of the input file. [default: None]
  --out-format=OUTFMT
              Format of the output file. [default: pmd]
  --specorder=SPECORDER
              Order of species. [default: Al,Mg,Si]
  --strain=STRAIN
              Strain to be applied. [default: 0.01]
"""
from __future__ import print_function

import os,sys,copy
import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem, parse_filename

def isotropic(psys0,strain=0.01):
    fac = 1.0 + strain
    psys = copy.deepcopy(psys0)
    alc0 = psys0.alc
    alc = alc0 *fac
    psys.alc = alc
    return psys


def uniaxial(psys0,strain=0.01):
    fac = 1.0 + strain
    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)
    
    psyslist = []

    psys = copy.deepcopy(psys0)
    hmat = psys0.get_hmat()
    hmat[0,:] *= fac
    psys.set_hmat(hmat)
    psyslist.append(psys)

    if rb != ra:
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[1,:] *= fac
        psys.set_hmat(hmat)
        psyslist.append(psys)
    if rc != ra and rc != rb:
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[2,:] *= fac
        psys.set_hmat(hmat)
        psyslist.append(psys)
    return psyslist


def shear(psys0,strain=0.01):
    fac = 1.0 + strain
    aa,ba,ca = psys0.get_lattice_angles()
    raa = round(aa,_angle_digit)
    rba = round(ba,_angle_digit)
    rca = round(ca,_angle_digit)

    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)

    psyslist = []
    dlt = (fac-1.0)*(a+b)/2/2
    psys = copy.deepcopy(psys0)
    hmat = psys0.get_hmat()
    hmat[0,1] += dlt
    hmat[1,0] += dlt
    psys.set_hmat(hmat)
    psyslist.append(psys)

    if not (rba == rca and set((ra,rc)) == set((ra,rb))):
        dlt = (fac-1.0)*(a+c)/2/2
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[0,2] += dlt
        hmat[2,0] += dlt
        psys.set_hmat(hmat)
        psyslist.append(psys)
        
    if not (raa == rca and set((rb,rc)) == set((ra,rb))) and \
       not (raa == rba and set((rb,rc)) == set((rb,rc))) :
        dlt = (fac-1.0)*(b+c)/2/2
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[1,2] += dlt
        hmat[2,1] += dlt
        psys.set_hmat(hmat)
        psyslist.append(psys)
    return psyslist


if __name__ == "__main__":

    args= docopt(__doc__)

    infmt= args['--in-format']
    outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']
    specorder= args['--specorder'].split(',')
    strain= float(args['--strain'])

    print('args:')
    print(args)

    psys0= NAPSystem(fname=infname,ffmt=infmt,specorder=specorder)

    if not outfmt == None:
        outfmt= parse_filename(outfname)

    if args['isotropic']:
        psys = isotropic(psys0,strain)
    elif args['uniaxial']:
        psys = uniaxial(psys0,strain)[0]
    elif args['shear']:
        psys = shear(psys0,strain)[0]

    if outfmt == 'pmd':
        psys.write_pmd(outfname)
    elif outfmt == 'smd':
        psys.write_pmd(outfname)
    elif outfmt == 'akr':
        psys.write_akr(outfname)
    elif outfmt == 'POSCAR':
        psys.write_POSCAR(outfname)
    elif outfmt == 'dump':
        psys.write_dump(outfname)
    elif outfmt == 'xsf':
        psys.write_xsf(outfname)
    else:
        print('Cannot detect output file format.')
    
