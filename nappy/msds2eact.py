#!/usr/bin/env python
"""
Compute activation energy from temperature vs diffusion coefficient.
DIRS should be like '300K' and the digits before 'K' is used as temperature.

Usage:
  msds2eact.py [options] DIRS [DIRS...]

Options:
  -h, --help  Show this message and exit.
  -o OUT      Output file name of D vs T. [default: out.D-T]
  -t DT       Time interval between successive data points in fs. [default: 1.0]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

from nappy.msd2diff import read_out_msd, msd2diff

__author__ = "RYO KOBAYASHI"
__version__ = "181008"

_kB = 8.6173303e-5

if __name__ == "__main__":

    args = docopt(__doc__)
    dirs = args['DIRS']
    dt = float(args['-t'])
    outfname = args['-o']

    #...Sort dirs list in numerical order
    dirs.sort(cmp=lambda x,y: cmp(int(x.replace('K','')), int(y.replace('K',''))),
              reverse=True)

    Ts = []
    Ds = []
    for d in dirs:
        T = d.replace('K','')
        print(' T = ',T)
        ts,msds = read_out_msd(d+'/out.msd',dt)
        dcoeff = msd2diff(ts,msds)
        Ts.append(float(T))
        Ds.append(dcoeff)

    with open(outfname,'w') as f:
        f.write('# temperature [K], diffusion coefficent [m^2/sec]\n')
        for T,D in zip(Ts,Ds):
            f.write('    {0:8.2f}     {1:12.4e}\n'.format(T,D))
    print(' Wrote {0:s}'.format(outfname))
    
    Tinvs = [ 1.0/T for T in Ts ]
    logDs = [ np.log(D) for D in Ds ]
    A = np.array([Tinvs, np.ones(len(Tinvs))])
    A = A.T
    a,b = np.linalg.lstsq(A, logDs, rcond=None)[0]
    
    Eact = a *(-_kB)
    print(' Activation energy = {0:.3f} eV'.format(Eact))
