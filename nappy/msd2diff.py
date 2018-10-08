#!/usr/bin/env python
"""
Compute diffusion coefficient from MSD data.

Usage:
  msd2diff.py [options] MSD_FILE

Options:
  -h, --help  Show this message and exit.
  -t DT       Time interval between successive data points in fs. [default: 1.0]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "181008"


def read_out_msd(fname='out.msd',dt=1.0):
    with open(fname,'r') as f:
        lines = f.readlines()
    ts = []
    msds = []
    for il,line in enumerate(lines):
        if line[0] == '#':
            continue
        data = line.split()
        n = int(data[0])
        msd = float(data[1])
        ts.append(n*dt)
        msds.append(msd)
    return ts,msds

def msd2diff(ts,msds):
    """
    Compute diffusion coefficient from time [fs] vs MSD [m^2] data 
    by solving least square problem using numpy.
    Return diffusion coefficient in [m^2/sec].
    """
    A= np.array([ts, np.ones(len(ts))])
    A = A.T
    a,b = np.linalg.lstsq(A,msds,rcond=None)[0]
    dcoeff = a *1.0e-20 /1.e-15
    return dcoeff

if __name__ == "__main__":

    args = docopt(__doc__)
    dt = float(args['-t'])
    fname = args['MSD_FILE']

    ts,msds = read_out_msd(fname,dt)

    #...Least square
    dcoeff = msd2diff(ts,msds)
    print('Diffusion coefficient = {0:12.4e} [m^2/s]'.format(dcoeff))
