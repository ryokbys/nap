#!/usr/bin/env python
"""
Deform the original system by adding some strain.

Usage:
  deform.py isotropic [options] INFILE OUTFILE
  deform.py uniaxial  [options] INFILE OUTFILE
  deform.py shear     [options] INFILE OUTFILE

Options:
  -h, --help  Show this help message and exit.
  --specorder=SPECORDER
              Order of species. [default: None]
  --factor=FACTOR
              Strain factor to be applied. [default: 0.95]
"""
from __future__ import print_function

import os,sys,copy
import numpy as np
from docopt import docopt

#from nappy.napsys import NAPSystem, parse_filename
#from nappy.io import read, write, parse_filename
import nappy

def isotropic(psys0,factor=1.0):
    psys = copy.deepcopy(psys0)
    alc0 = psys0.alc
    alc = alc0 *factor
    psys.alc = alc
    return psys


def uniaxial(psys0,factor=1.0):
    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)
    
    psyslist = []

    psys = copy.deepcopy(psys0)
    hmat = psys0.get_hmat()
    hmat[0,:] *= factor
    psys.set_hmat(hmat)
    psyslist.append(psys)

    if rb != ra:
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[1,:] *= factor
        psys.set_hmat(hmat)
        psyslist.append(psys)
    if rc != ra and rc != rb:
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[2,:] *= factor
        psys.set_hmat(hmat)
        psyslist.append(psys)
    return psyslist


def shear(psys0,factor=1.0):

    aa,ba,ca = psys0.get_lattice_angles()
    raa = round(aa,_angle_digit)
    rba = round(ba,_angle_digit)
    rca = round(ca,_angle_digit)

    a,b,c = psys0.get_lattice_lengths()
    ra = round(a,_length_digit)
    rb = round(b,_length_digit)
    rc = round(c,_length_digit)

    psyslist = []
    dlt = (factor-1.0)*(a+b)/2/2
    psys = copy.deepcopy(psys0)
    hmat = psys0.get_hmat()
    hmat[0,1] += dlt
    hmat[1,0] += dlt
    psys.set_hmat(hmat)
    psyslist.append(psys)

    if not (rba == rca and set((ra,rc)) == set((ra,rb))):
        dlt = (factor-1.0)*(a+c)/2/2
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[0,2] += dlt
        hmat[2,0] += dlt
        psys.set_hmat(hmat)
        psyslist.append(psys)
        
    if not (raa == rca and set((rb,rc)) == set((ra,rb))) and \
       not (raa == rba and set((rb,rc)) == set((rb,rc))) :
        dlt = (factor-1.0)*(b+c)/2/2
        psys = copy.deepcopy(psys0)
        hmat = psys0.get_hmat()
        hmat[1,2] += dlt
        hmat[2,1] += dlt
        psys.set_hmat(hmat)
        psyslist.append(psys)
    return psyslist

def main():
    import os,sys
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])))

    # infmt= args['--in-format']
    # outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']
    specorder= args['--specorder'].split(',')
    if specorder[0] == 'None':
        specorder = None
    factor= float(args['--factor'])

    psys0 = nappy.io.read(fname=infname,specorder=specorder)

    if args['isotropic']:
        psys = isotropic(psys0,factor)
    elif args['uniaxial']:
        psys = uniaxial(psys0,factor)[0]
    elif args['shear']:
        psys = shear(psys0,factor)[0]

    nappy.io.write(psys,fname=outfname)

    return None
    

if __name__ == "__main__":

    main()
