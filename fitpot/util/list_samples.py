#!/usr/bin/env python
"""

Usage:
  list_samples.py [options] DIR

Options:
  -h,--help  Show this message and exit.

"""
from __future__ import print_function

import os,sys
from glob import glob
from docopt import docopt

def uniq(arr):
    newarr = []
    for a in arr:
        if not a in newarr:
            newarr.append(a)
    return newarr


if __name__ == "__main__":

    args = docopt(__doc__)

    dname = args['DIR']
    
    smpls = glob(dname+"/smpl_*")
    for i in range(len(smpls)):
        smpls[i] = smpls[i].split('/')[-1][5:-6]
    uniq_smpls = uniq(smpls)
    uniq_smpls.sort()
    print('sample_error   {0:d}'.format(len(uniq_smpls)))
    for s in uniq_smpls:
        print('{0:>20s}  0.001  0.1'.format(s))


