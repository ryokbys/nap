#!/usr/bin/env python
"""
Create `in.params.linreg` from `in.params.desc`.

Usage:
  desc2linreg.py [options]

Options:
  -h, --help  Show this message and exit.
  -i IFNAME   Input file name. [default: in.params.desc]
  -o OFNAME   Output file name. [default: in.params.linreg]
  -m MAX      Max value of the weights that are randomly determined. [default: 1.0e+10]
"""
import numpy as np
from docopt import docopt
import random
import nappy.nn.desc as desc

__author__ = "RYO KOBAYASHI"
__version__ = "180807"


def create_linreg_input(nsf,desc,wmax,fname='in.params.linreg'):
    rc2 = 2.0
    rc3 = 1.0
    for d in descs:
        sftype = d[0]
        if 'angular' in sftype:
            rc3 = np.ceil(max(rc3,d[4]))
        elif sftype in ('gauss','cosine'):
            rc2 = np.ceil(max(rc2,d[3]))
    
    with open(fname,'w') as f:
        f.write(' {0:4d}  {1:7.3f}  {2:7.3f}\n'.format(nsf,rc2,rc3))
        for i in range(nsf):
            w = (random.random() -0.5)*2.0 *min(1.0e-2,wmax)
            f.write('  {0:12.4e}  {1:12.4e}  {2:12.4e}\n'.format(w,-wmax,wmax))
    print(' Wrote {0:s}'.format(fname))
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['-i']
    outfname = args['-o']
    wmax = float(args['-m'])

    nsp,nsf,descs,r_inner = desc.read_desc(infname)

    create_linreg_input(nsf,descs,wmax,outfname)
    
    
