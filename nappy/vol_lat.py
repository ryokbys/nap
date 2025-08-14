#!/usr/bin/env python
"""
Extract averaged volume and lattice parameters from files.

Usage:
  {0:s} [options] FILES [FILES...]

Options:
  -h, --help  Show this message and exit.
  --skip NSKIP
              Skip first NSKIP steps from the statistics.
              If this is -1, vol and lat of the final step are taken. [default: 0]
  --out4fp    Flag to write out in general fp.py format. [default: Fault]
  --prefix PREFIX
              Prefix for output files. [default: data.pmd]
  --format FORMAT
              Input file format. [default: extxyz]
"""
import os,sys
from docopt import docopt
import numpy as np

import nappy
from nappy.common import get_key

__author__ = "Ryo KOBAYASHI"
__version__ = "250814"

def nsys2lat(nsys):
    a,b,c = nsys.get_lattice_lengths()
    alpha,beta,gamma = nsys.get_lattice_angles()
    alpha = alpha /np.pi *180.0
    beta  = beta  /np.pi *180.0
    gamma = gamma /np.pi *180.0
    return a,b,c,alpha,beta,gamma

def main():

    from nappy.util import header_msg
    print(header_msg(sys.argv))

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])), version=__version__)

    files = args['FILES']
    if len(files) > 1:
        files.sort(key=get_key, reverse=True)
    nskip = int(args['--skip'])
    #del files[:nskip]
    prefix = args['--prefix']
    out4fp = args['--out4fp']
    fmt = args['--format']
    if fmt == 'None':
        fmt = None

    nsum = 0
    volsum = 0.0
    asum= bsum= csum= 0.0
    alpsum= betsum= gmmsum= 0.0
    for i,fi in enumerate(files):
        assert os.path.exists(fi), f"[Error] File, {fi}, does not exist !!!"
        # if not os.path.exists(fi):
        #     print(f"[Error] File, {fi}, does not exist !!!")
        #     sys.exit()
        print(' File =',fi)
        nsys = nappy.io.read(fname=fi, format=fmt)
        if type(nsys) is list:
            for nsysi in nsys:
                volsum += nsysi.get_volume()
                a,b,c,alpha,beta,gamma = nsys2lat(nsysi)
                asum += a
                bsum += b
                csum += c
                alpsum += alpha
                betsum += beta
                gmmsum += gamma
                nsum += 1
        else:
            volsum += nsysi.get_volume()
            a,b,c,alpha,beta,gamma = nsys2lat(nsysi)
            asum += a
            bsum += b
            csum += c
            alpsum += alpha
            betsum += beta
            gmmsum += gamma
            nsum += 1

    assert nsum > 0, 'NSUM==0, which should not happe!'
    vol = volsum /nsum
    a = asum/nsum
    b = bsum/nsum
    c = csum/nsum
    alpha = alpsum/nsum
    beta  = betsum/nsum
    gamma = gmmsum/nsum

    #...Regardless of prefix, write out.vol and out.lat
    with open('out.vol','w') as f:
        f.write('{0:15.3f}\n'.format(vol))
    with open('out.lat','w') as f:
        f.write(' {0:10.3f} {1:10.3f} {2:10.3f}'.format(a,b,c)
                +' {0:10.3f} {1:10.3f} {2:10.3f}\n'.format(alpha,beta,gamma))
    #...Format of output (named by prefix) depends on out4fp
    if out4fp:
        with open(prefix+'.vol','w') as f:
            f.write('# Volume\n')
            f.write('     1     1.0\n')
            f.write('{0:15.3f}\n'.format(vol))
        with open(prefix+'.lat','w') as f:
            f.write('# Lattice parameters\n')
            f.write('     6     1.0\n')
            f.write(' {0:10.3f} {1:10.3f} {2:10.3f}'.format(a,b,c)
                    +' {0:10.3f} {1:10.3f} {2:10.3f}\n'.format(alpha,beta,gamma))
    else:
        with open(prefix+'.vol','w') as f:
            f.write('{0:15.3f}\n'.format(vol))
        with open(prefix+'.lat','w') as f:
            f.write(' {0:10.3f} {1:10.3f} {2:10.3f}'.format(a,b,c)
                    +' {0:10.3f} {1:10.3f} {2:10.3f}\n'.format(alpha,beta,gamma))

    print(' Wrote {0:s}.vol {0:s}.lat'.format(prefix))


if __name__ == "__main__":
    main()
