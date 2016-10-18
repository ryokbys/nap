#!/usr/bin/env python
"""
Take staggered averages of all the columns.

Usage:
  staggered_average.py [options] FILE

Options:
  -h, --help  Show this message and exit.
  --width WIDTH
              Averaging window width. [default: 100]
  --offset OFFSET
              Offset of staggered average. [default: 10]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = ""

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['FILE']
    width = int(args['--width'])
    offset = int(args['--offset'])

    #...1st analyze data
    f = open(fname,'r')
    lines = f.readlines()
    nrow = len(lines)
    ncol = len(lines[-1].split())
    ndata = nrow/offset
    f.close()

    print(' nrow  = ',nrow)
    print(' ncol  = ',ncol)
    print(' ndata = ',ndata)
    data = np.zeros((ndata,ncol),dtype=float)
    
    f = open(fname,'r')
    fo = open(fname+'.ave','w')
    #...write out 1st line
    il = 0
    for line in lines:
        if line[0] == '#':
            fo.write(line)
            continue
        # if il == 0:
        #     fo.write(line)
        for i in range(ndata):
            if i*offset <= il < i*offset + width:
                data[i,:] += [ float(x) for x in line.split() ]
        il += 1
        if il >= nrow-width:
            lastdata = (il-width)/offset
            break
    f.close()

    data[:,:] /= width
    for i in range(min(lastdata,ndata)):
        fo.write(' {0:10.1f}'.format(data[i,0]))
        for ii in range(1,len(data[i])):
            fo.write(' {0:7.2f}'.format(data[i,ii]))
        fo.write('\n')
    fo.close()
    print(' check '+fname+'.ave')
