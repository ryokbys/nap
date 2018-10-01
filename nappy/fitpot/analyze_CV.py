#!/usr/bin/env python
"""
Analyze cross validation (CV) results in CV_* directories.
Output file of `fitpot` is assumed to be `out.fitpot`.

Usage:
  analyze_CV.py [options] CV_DIRS [CV_DIRS...]

Options:
  -h, --help  Show this message and exit.
  --which WHICH
              Which value is to be used: iter, energy or force. [default: energy]
  --threshold THRESHOLD
              Threshold for the outlier data. [default: 1.0]
"""
from __future__ import print_function

from docopt import docopt
import numpy as np
from nappy.fitpot.extract_best import extract_best

__author__ = "RYO KOBAYASHI"
__version__ = "180906"

def extract_varfname(fname='out.fitpot'):
    with open(fname,'r') as f:
        for line in f.readlines():
            if 'param_file' in line:
                paramfname = line.split()[1]
                return paramfname
    raise EOFError('{0:s} does not include the line contains param_file.'.format(fname))

def weight_norm(fname='in.vars.fitpot'):
    with open(fname,'r') as f:
        lines = f.readlines()
        nwgts = int(lines[0].split()[0])
    wgts = np.zeros(nwgts,dtype=float)
    for il in range(1,nwgts+1):
        data = lines[il].split()
        wgts[il-1] = float(data[0])
    return np.linalg.norm(wgts)

if __name__ == "__main__":

    args = docopt(__doc__)
    cv_dirs = args['CV_DIRS']
    which = args['--which']
    threshold = float(args['--threshold'])

    cv_dirs.sort()
    data_trn = []
    data_tst = []
    wgt_norm = []
    print('Directory, min(training), min(test), wgt_norm')
    for d in cv_dirs:
        fitmin,ftrnmin,ftstmin,eitmin,etrnmin,etstmin,frcitmin,frctrnmin,frctstmin = extract_best(d+'/out.fitpot')
        wgtfname = extract_varfname(d+'/out.fitpot')
        if which == 'iter':
            wnrm = weight_norm(d+'/'+wgtfname+'.{0:d}'.format(fitmin))
            print('  {0:s}  {1:12.4e}  {2:12.4e}  {3:12.4e}'.format(d,ftrnmin,ftstmin,wnrm))
            if ftrnmin > threshold or ftstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(ftrnmin)
                data_tst.append(ftstmin)
                wgt_norm.append(wnrm)
        elif which == 'energy':
            wnrm = weight_norm(d+'/'+wgtfname+'.{0:d}'.format(eitmin))
            print('  {0:s}  {1:12.5f}  {2:12.5f}  {3:12.4e}'.format(d,etrnmin,etstmin,wnrm))
            if etrnmin > threshold or etstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(etrnmin)
                data_tst.append(etstmin)
                wgt_norm.append(wnrm)
        elif which == 'force':
            wnrm = weight_norm(d+'/'+wgtfname+'.{0:d}'.format(frcitmin))
            print('  {0:s}  {1:12.5f}  {2:12.5f}  {3:12.4e}'.format(d,frctrnmin,frctstmin,wnrm))
            if frctrnmin > threshold or frctstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(frctrnmin)
                data_tst.append(frctstmin)
                wgt_norm.append(wnrm)

    trnarray = np.array(data_trn)
    tstarray = np.array(data_tst)
    wgtarray = np.array(wgt_norm)
    print('Mean, std, max, and min:')
    dmean = trnarray.mean()
    dstd = trnarray.std()
    dmax = trnarray.max()
    dmin = trnarray.min()
    print('  Training= {0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f}'.format(dmean,dstd,
                                                                       dmax,dmin))
    dmean = tstarray.mean()
    dstd  = tstarray.std()
    dmax  = tstarray.max()
    dmin  = tstarray.min()
    print('  Test=     {0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f}'.format(dmean,dstd,
                                                                       dmax,dmin))
    
    dmean = wgtarray.mean()
    dstd  = wgtarray.std()
    dmax  = wgtarray.max()
    dmin  = wgtarray.min()
    print('  Weight=     {0:12.4e} {1:12.4e} {2:12.4e} {3:12.4e}'.format(dmean,dstd,
                                                                         dmax,dmin))
