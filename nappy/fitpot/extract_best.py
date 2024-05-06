#!/usr/bin/env python
"""
Extract best values from ``out.fitpot`` file.

Usage:
  extract_best.py [options] OUT_FITPOT

Options:
  -h, --help  Show this message and exit.
  --which WHICH
              Which best is used: training or test. [default: test]
"""
import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "180906"

def extract_best(fname='out.fitpot',which_best='test'):
    """
    Extract best values of iter,RMSE(energy),RMSE(force) from out.fitpot.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    fvals = []
    ergs = []
    frcs = []
    for il,line in enumerate(lines):
        if 'iter' in line and 'ftrn' in line:
            data = line.split()
            it = int(data[1])
            if 'niter' in line:
                ftrn = float(data[3])
                ftst = float(data[4])
            else:
                ftrn = float(data[2])
                ftst = float(data[3])
            fvals.append((it,ftrn,ftst))
        elif 'ENERGY:' in line:
            data = line.split()
            if data[0] == '#':
                continue
            it = int(data[1])
            etrn = float(data[3])
            etst = float(data[5])
            ergs.append((it,etrn,etst))
        elif 'FORCE:' in line:
            if data[0] == '#':
                continue
            data = line.split()
            it = int(data[1])
            try:
                frctrn = float(data[3])
            except:
                frctrn = 10.0
            try:
                frctst = float(data[5])
            except:
                frctst = 10.0
            frcs.append((it,frctrn,frctst))

    fmin = 1.0e+30
    for i in range(len(fvals)):
        it,ftrn,ftst = fvals[i]
        if which_best == 'test':
            f = ftst
        else:
            f = ftrn
        if f < fmin:
            fitmin = it
            ftrnmin = ftrn
            ftstmin = ftst
            if which_best == 'test':
                fmin = ftstmin
            else:
                fmin = ftrnmin
                
    emin = 1.0e+30
    for i in range(len(ergs)):
        it,etrn,etst = ergs[i]
        if which_best == 'test':
            e = etst
        else:
            e = etrn
        if e < emin:
            eitmin = it
            etrnmin = etrn
            etstmin = etst
            if which_best == 'test':
                emin = etstmin
            else:
                emin = etrnmin

    frcmin = 1.0e+30
    for i in range(len(frcs)):
        it,frctrn,frctst = frcs[i]
        if which_best == 'test':
            frc = frctst
        else:
            frc = frctrn
        if frc < frcmin:
            frcitmin = it
            frctrnmin = frctrn
            frctstmin = frctst
            if which_best == 'test':
                frcmin = frctstmin
            else:
                frcmin = frctrnmin

    return fitmin,ftrnmin,ftstmin, eitmin,etrnmin,etstmin,frcitmin,frctrnmin,frctstmin
    

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['OUT_FITPOT']
    which_best = args['--which']

    fitmin,ftrnmin,ftstmin, eitmin,etrnmin,etstmin,frcitmin,frctrnmin,frctstmin = extract_best(fname,which_best)

    print('Min fvalue      = {0:6d}  {1:12.4e}  {2:12.4e}'.format(fitmin,ftrnmin,ftstmin))
    print('Min RMSE energy = {0:6d}  {1:12.4e}  {2:12.4e}'.format(eitmin,etrnmin,etstmin))
    print('Min RMSE force  = {0:6d}  {1:12.4e}  {2:12.4e}'.format(frcitmin,frctrnmin,frctstmin))
    
