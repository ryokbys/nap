#!/usr/bin/env python
"""
Compute diffusion coefficient from MSD data.

Usage:
  msd2diff.py [options] MSD_FILE

Options:
  -h, --help  Show this message and exit.
  -o, --offset OFFSET
              Offset of given data. [default: 0]
  -t DT       Time interval between successive data points in fs. [default: 1.0]
  --plot      Plot a fitted graph. [default: False]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "181008"


def read_out_msd(fname='out.msd',dt=1.0,offset=0):
    with open(fname,'r') as f:
        lines = f.readlines()
    ts = []
    msds = []
    n0 = 0
    msd0 = 0.0
    for il,line in enumerate(lines):
        if line[0] == '#':
            continue
        data = line.split()
        if il < offset:
            n0 = int(data[0])
            msd0 = float(data[1])
            continue
        n = int(data[0])
        msd = float(data[1])
        ts.append((n-n0)*dt)
        msds.append(msd-msd0)
    return np.array(ts),np.array(msds)

def msd2diff(ts,msds,fac):
    """
    Compute diffusion coefficient from time [fs] vs MSD [m^2] data 
    by solving least square problem using numpy.
    Return diffusion coefficient multiplied by FAC.
    """
    A= np.array([ts, np.ones(len(ts))])
    A = A.T
    xvar = np.var(A[:,0])
    p,res,_,_ = np.linalg.lstsq(A,msds,rcond=None)
    a = p[0]
    b = p[1]
    # fac = 1.0e-16 /1.e-15
    a = a *fac
    b = b *fac
    # print(res[0],xvar,np.mean(A[:,0]),len(ts))
    std = np.sqrt(res[0]/len(ts)/xvar) *fac 
    return a,b,std

if __name__ == "__main__":

    args = docopt(__doc__)
    dt = float(args['-t'])
    fname = args['MSD_FILE']
    offset = int(args['--offset'])
    plot = args['--plot']

    ts,msds = read_out_msd(fname,dt,offset)
    #...Assuming input MSD unit in A^2/fs and output in cm^2/s
    fac = 1.0e-16 /1.0e-15
    #...Least square
    a,b,std = msd2diff(ts,msds,fac)
    print(' Diffusion coefficient = {0:12.4e}'.format(a)+
          ' +/- {0:12.4e} [cm^2/s]'.format(std))

    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
        #...Original time unit == fs
        unit = 'fs'
        tfac = 1.0
        if ts[-1] > 1.0e+5: #...if max t > 100ps, time unit in ps
            unit = 'ps'
            tfac = 1.0e-3
        plt.xlabel('Time ({0:s})'.format(unit))
        plt.ylabel('MSD (A^2/{0:s})'.format(unit))
        fvals = np.array([ (t*a+b)/fac for t in ts ])
        plt.plot(ts*tfac,msds/tfac,'b-',label='MSD data')
        plt.plot(ts*tfac,fvals/tfac,'r-',label='Fitted curve')
        plt.savefig("graph_msd2diff.png", format='png',
                    dpi=300, bbox_inches='tight')
        print(' Wrote graph_msd2diff.png')
        
