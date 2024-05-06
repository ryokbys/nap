#!/usr/bin/env python
"""
Extract average volume and lattice parameters from vasprun.xml.

Usage:
  vasprun2vol.py [options] VASPRUN.XML

Options:
  -h, --help    Show this message and exit.
  --skip NSKIP  Skip first NSKIP steps from the statistics. [default: 0]
"""
import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "rev190906"

def extract_vol_lat(vasprunfile='vasprun.xml',nskip=0):

    with open(vasprunfile,'r') as f:
        nvol = 0
        vsum = 0.0
        nlat = 0
        mode = None
        lsum = {'a':0.0, 'b':0.0, 'c':0.0,
                'alpha':0.0, 'beta':0.0, 'gamma':0.0}
        for il,line in enumerate(f.readlines()):
            data = line.split()
            if len(data) == 0:
                mode = None
                continue
            if 'volume' in line:
                mode = None
                nvol += 1
                if nvol <= nskip:
                    continue
                vsum += float(data[2])
            elif '<varray name="basis" >' in line:
                mode = 'basis'
                ib = 0
                vecs = np.zeros((3,3))
            elif mode == 'basis' and '</varray>' in line:
                mode = None
                nlat += 1
                if nlat <= nskip:
                    continue
                a = np.linalg.norm(vecs[0])
                b = np.linalg.norm(vecs[1])
                c = np.linalg.norm(vecs[2])
                alpha = np.arccos(np.dot(vecs[1],vecs[2])/b/c) /np.pi *180.0
                beta  = np.arccos(np.dot(vecs[0],vecs[2])/a/c) /np.pi *180.0
                gamma = np.arccos(np.dot(vecs[0],vecs[1])/a/b) /np.pi *180.0
                lsum['a'] += a
                lsum['b'] += b
                lsum['c'] += c
                lsum['alpha'] += alpha
                lsum['beta'] += beta
                lsum['gamma'] += gamma
            elif mode == 'basis' and '<v>' in line:
                vecs[ib,:] = [ float(x) for x in data[1:4] ]
                ib += 1
            else:
                mode = None
    #...average
    vol = 0.0
    if nvol > nskip:
        vol = vsum /(nvol-nskip)
    else:
        raise ValueError('There is no vol data...')
    lat = {}
    if nlat > nskip:
        for k,v in lsum.items():
            lat[k] = v /(nlat-nskip)
    return vol, lat

def main(args):

    fname = args['VASPRUN.XML']
    nskip = int(args['--skip'])

    vol, lat = extract_vol_lat(fname,nskip)

    with open('data.ref.vol','w') as f:
        f.write('{0:15.3f}\n'.format(vol))

    lat_key_order = ('a','b','c','alpha','beta','gamma')
    with open('data.ref.lat','w') as f:
        for k in lat_key_order:
            f.write(' {0:10.3f}'.format(lat[k]))
        f.write('\n')

    print('Wrote data.ref.vol data.ref.lat')

if __name__ == "__main__":

    args = docopt(__doc__)

    main(args)
