#!/usr/bin/env python
"""
Extract averaged volumen and lattice parameters from dump files.

Usage:
  dumps2vol.py [options] DUMPS [DUMPS...]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

from nappy.napsys import NAPSystem

__author__ = "RYO KOBAYASHI"
__version__ = "rev190906"

def nsys2lat(nsys):
    a,b,c = nsys.get_lattice_lengths()
    alpha,beta,gamma = nsys.get_lattice_angles()
    alpha = alpha /np.pi *180.0
    beta  = beta  /np.pi *180.0
    gamma = gamma /np.pi *180.0
    return a,b,c,alpha,beta,gamma

def main(args):

    dumps = args['DUMPS']
    nsum = 0
    volsum = 0.0
    asum= bsum= csum= 0.0
    alpsum= betsum= gmmsum= 0.0
    for dump in dumps:
        try:
            nsys = NAPSystem(fname=dump)
            volsum += nsys.volume()
            a,b,c,alpha,beta,gamma = nsys2lat(nsys)
            asum += a
            bsum += b
            csum += c
            alpsum += alpha
            betsum += beta
            gmmsum += gamma
            nsum += 1
        except:
            print('Failed {0:s} '.format(dump))
            pass

    if nsum < 1:
        raise ValueError('Something went wrong...')
    vol = volsum /nsum
    a = asum/nsum
    b = bsum/nsum
    c = csum/nsum
    alpha = alpsum/nsum
    beta  = betsum/nsum
    gamma = gmmsum/nsum

    with open('data.pmd.vol','w') as f:
        f.write('{0:15.3f}\n'.format(vol))

    with open('data.pmd.lat','w') as f:
        f.write(' {0:10.3f} {1:10.3f} {2:10.3f}'.format(a,b,c)
                +' {0:10.3f} {1:10.3f} {2:10.3f}\n'.format(alpha,beta,gamma))

    print('Wrote data.pmd.vol data.pmd.lat')

    
if __name__ == "__main__":

    args = docopt(__doc__)

    main(args)
