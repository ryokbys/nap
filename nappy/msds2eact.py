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
  --dim DIM   Spatial dimension of diffusion. [default: 3]
  --offset OFFSET
              Offset of the data. [default: 0]
  --plot
              Plot E_act vs 1/T graph (optional) or write a gnuplot script. [default: False]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
from scipy import stats

from nappy.msd2diff import read_out_msd, msd2D

__author__ = "RYO KOBAYASHI"
__version__ = "181008"

_kB = 8.6173303e-5

def make_gnuplot_file(outDT,Eact,D0):
    txt = """set format y '10^{{%L}}'
set xl '1000/T (1/K)'
set yl 'D [cm^2/sec]'
set log y
f(x) = exp({0:.3f} -{1:.3f} /8.617e-5 *(x/1000))
plot '{2:s}' us (1000.0/$1):2:3 w yerr lc 'blue' pt 7 t 'data', f(x) t 'fitted' lc 'blue'
""".format(D0,Eact,outDT)
    with open('plot_D-T.gp','w') as f:
        f.write(txt)
    print(' Wrote plot_D-T.gp')
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    dirs = args['DIRS']
    dt = float(args['-t'])
    outfname = args['-o']
    offset = int(args['--offset'])
    plot = args['--plot']
    dim = int(args['--dim'])

    #...Sort dirs list in numerical order
    dirs.sort(cmp=lambda x,y: cmp(int(x.replace('K','')), int(y.replace('K',''))),
              reverse=True)

    Ts = np.zeros(len(dirs))
    Ds = np.zeros(len(dirs))
    Dstds = np.zeros(len(dirs))
    fac = 1.0e-16 /1.0e-15 #...A^2/fs to cm^2/s
    for i,d in enumerate(dirs):
        T = d.replace('K','')
        ts,msds = read_out_msd(d+'/out.msd',dt,offset)
        D,b,Dstd = msd2D(ts,msds,fac,dim=dim)
        print(' T,D = {0:5d}K, {1:12.4e} +/- {2:12.4e} [cm^2/s]'.format(int(T),D,Dstd))
        Ts[i] = float(T)
        Ds[i] = D
        Dstds[i] = Dstd

    with open(outfname,'w') as f:
        f.write('# T [K],       D [cm^2/sec],    sgm(D) [cm^2/sec]\n')
        for T,D,Dstd in zip(Ts,Ds,Dstds):
            f.write('    {0:8.2f}  {1:12.4e}  {2:12.4e}\n'.format(T,D,Dstd))
    print(' Wrote {0:s}'.format(outfname))
    
    Tinvs = np.array([ 1.0/T for T in Ts ])
    logDs = np.array([ np.log(D) for D in Ds ])
    a,b,r,p,stderr = stats.linregress(Tinvs,logDs)
    
    Eact = a *(-_kB)
    Eaerr = stderr *_kB
    print(' Ea = {0:.3e} +/- {1:.3e} [eV]'.format(Eact,Eaerr))
    print(' D0 = {0:.4e} [cm^2/s]'.format(np.exp(b)))

    make_gnuplot_file(outfname,Eact,b)
    
    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
        plt.xlabel('1000/T [1/K]')
        plt.ylabel('log(D [cm^2/sec])')
        plt.plot(Tinvs*1000,logDs,'bo',label='data')
        fvals = np.array([ b +a*(1.0/T) for T in Ts])
        plt.plot(Tinvs*1000,fvals,'r-',label='fitted')
        plt.savefig("graph_msds2eact.png", format='png',
                    dpi=300, bbox_inches='tight')
