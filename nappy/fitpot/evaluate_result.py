#!/usr/bin/env python
"""
Evaluate the fitpot result.

Usage:
  evaluate_result.py histogram [options] DATAFILE

Options:
  -h,--help  Show this message and exit.
  -w,--width=WIDTH
             Width of the bin in histogram in eV. [default: 0.0005]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np


__author__ = "RYO KOBAYASHI"
__version__ = ""


def histogram(fname="",width=0.0005):
    
    with open(fname,'r') as f:
        lines = f.readlines()
        data = np.zeros((len(lines),2),dtype=float)
        for il,l in enumerate(lines):
            d = l.split()
            data[il,:] = [ float(d[0]),float(d[1]) ]

    diffmax = 0.0
    for i in range(len(data)):
        diff = abs(data[i,0]-data[i,1])
        diffmax = max(diff,diffmax)
    ndim = int(diffmax/width)+1
    hg = np.zeros((ndim),dtype=int)
    for i in range(len(data)):
        diff = abs(data[i,0]-data[i,1])
        ihg = int(diff/width)
        hg[ihg] += 1
    print("# histogram of "+fname
          +", width={0:10.5f}".format(width))
    ld = len(data)
    for i in range(len(hg)):
        print("{0:10.5f} {1:6d} {2:8.2f}".format(i*width+0.5*width,
                                                 hg[i],
                                                 float(hg[i])/ld))

if __name__ == "__main__":

    args = docopt(__doc__)
    dfile = args['DATAFILE']
    width = float(args['--width'])
    hist = args['histogram']

    if hist:
        histogram(dfile,width)
