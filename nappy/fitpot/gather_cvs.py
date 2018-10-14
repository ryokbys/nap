#!/usr/bin/env python
"""
Gather CV data from given directories.

Usage:
  gather_cvs.py [options] DIRS [DIRS...]

Options:
  -h, --help  Show this message and exit.
  --prefix PREFIX
              Prefix to be removed from given dirname. [default: ]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "181002"

def read_outcv(fname='out.cv'):
    with open(fname,'r') as f:
        lines = f.readlines()
    trns = []
    tsts = []
    wgts = []
    for il,line in enumerate(lines):
        if 'CV_' in line:
            data = line.split()
            trn = float(data[1])
            tst = float(data[2])
            wgt = float(data[3])
            trns.append(trn)
            tsts.append(tst)
            wgts.append(wgt)
    return trns,tsts,wgts

def get_stats(arr):
    nparr = np.array(arr)
    mean = nparr.mean()
    med = np.median(nparr)
    dmax = nparr.max()
    dmin = nparr.min()
    return mean,med,dmax,dmin

if __name__ == "__main__":

    args = docopt(__doc__)
    dirs = args['DIRS']
    prefix = args['--prefix']

    f_all = open('out.cvs.all','w')
    f_trn = open('out.cvs.trn','w')
    f_tst = open('out.cvs.tst','w')
    f_wgt = open('out.cvs.wgt','w')
    for d in dirs:
        trns,tsts,wgts = read_outcv(d+'/out.cv')
        dval = d.replace(prefix,'')
        for i in range(len(trns)):
            f_all.write('{0:s} {1:13.5e} '.format(dval,trns[i])
                        +'{0:13.5e} {1:13.5e}\n'.format(tsts[i],wgts[i]))
        trn_stats = get_stats(trns)
        tst_stats = get_stats(tsts)
        wgt_stats = get_stats(wgts)
        f_trn.write('{0:s} {1:13.5e} '.format(dval,trn_stats[0])
                    +'{0:13.5e} {1:13.5e} {2:13.5e}\n'.format(trn_stats[1],
                                                              trn_stats[2],
                                                              trn_stats[3]))
        f_tst.write('{0:s} {1:13.5e} '.format(dval,tst_stats[0])
                    +'{0:13.5e} {1:13.5e} {2:13.5e}\n'.format(tst_stats[1],
                                                              tst_stats[2],
                                                              tst_stats[3]))
        f_wgt.write('{0:s} {1:13.5e} '.format(dval,wgt_stats[0])
                    +'{0:13.5e} {1:13.5e} {2:13.5e}\n'.format(wgt_stats[1],
                                                              wgt_stats[2],
                                                              wgt_stats[3]))

    f_all.close()
    f_trn.close()
    f_tst.close()
    f_wgt.close()
    os.system('sort -g -k1,1 -o {0:s} {0:s}'.format('out.cvs.all'))
    os.system('sort -g -k1,1 -o {0:s} {0:s}'.format('out.cvs.trn'))
    os.system('sort -g -k1,1 -o {0:s} {0:s}'.format('out.cvs.tst'))
    os.system('sort -g -k1,1 -o {0:s} {0:s}'.format('out.cvs.wgt'))
    print('Check the following output files:')
    print('  - out.cvs.all')
    print('  - out.cvs.trn')
    print('  - out.cvs.tst')
    print('  - out.cvs.wgt')
