#!/usr/bin/env python
"""
Count number of samples of types.

Usage:
  count_samples.py [options] DIR

Options:
  -h,--help  Show this message and exit.

"""
from __future__ import print_function

import os,sys
from docopt import docopt
from glob import glob


if __name__ == "__main__":
    
    args = docopt(__doc__)
    smpldir = args['DIR']
    
    samples = glob(smpldir+'/smpl_*')
    print('Number of samples = {}'.format(len(samples)))
    
    for i in range(len(samples)):
        samples[i] = samples[i][:-5]

    uniqs = []
    num_uniqs= {}
    for s in samples:
        if not s in uniqs:
            uniqs.append(s)
            num_uniqs[s] = 1
        else:
            num_uniqs[s] += 1
    uniqs.sort()
    for s in uniqs:
        print('{0:s}: {1:5d}'.format(s.replace(smpldir+'/',''),num_uniqs[s]))
