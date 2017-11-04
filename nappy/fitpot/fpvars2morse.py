#!/usr/bin/env python
"""
Convert fitpot parameters to in.params.Morse for pmd.

Usage:
  fp2morse.py [options] FITPOT_VAR_FILE

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = ""

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['FITPOT_VAR_FILE']

    ds = []
    alps = []
    rs = []
    with open(infname,'r') as f:
        lines = f.readlines()
        ndat = int(lines[0].split()[0])
        n = 0
        for i in range(1,len(lines)):
            n += 1
            if n > ndat: break
            ds.append(float(lines[3*(i-1)+1].split()[0]))
            n += 1
            if n > ndat: break
            alps.append(float(lines[3*(i-1)+2].split()[0]))
            n += 1
            if n > ndat: break
            rs.append(float(lines[3*(i-1)+3].split()[0]))

    # print(' ds  = ',ds)
    # print(' alps= ',alps)
    # print(' rs  = ',rs)
    with open('in.params.Morse','w') as f:
        f.write('# is, js,  D,      alpha,  rmin\n')
        for l in range(len(ds)):
            f.write(' {0:3d} {1:3d}'.format(1,l+2))
            f.write(' {0:7.3f} {1:7.3f} {2:7.3f}\n'.format(ds[l],alps[l],rs[l]))
            
    print(' Check in.params.Morse.')
    
