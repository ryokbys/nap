#!/usr/bin/env python
"""
Analyze cross validation (CV) results in CV_* directories.

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



if __name__ == "__main__":

    args = docopt(__doc__)
    cv_dirs = args['CV_DIRS']
    which = args['--which']
    threshold = float(args['--threshold'])

    cv_dirs.sort()
    data_trn = []
    data_tst = []
    for d in cv_dirs:
        fitmin,ftrnmin,ftstmin,eitmin,etrnmin,etstmin,frcitmin,frctrnmin,frctstmin = extract_best(d+'/out.fitpot')
        if which == 'iter':
            print('  {0:s}  {1:12.4e}  {2:12.4e}'.format(d,ftrnmin,ftstmin))
            if ftrnmin > threshold or ftstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(ftrnmin)
                data_tst.append(ftstmin)
        elif which == 'energy':
            print('  {0:s}  {1:12.5f}  {2:12.5f}'.format(d,etrnmin,etstmin))
            if etrnmin > threshold or etstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(etrnmin)
                data_tst.append(etstmin)
        elif which == 'force':
            print('  {0:s}  {1:12.5f}  {2:12.5f}'.format(d,frctrnmin,frctstmin))
            if frctrnmin > threshold or frctstmin > threshold:
                print('Not to register this data because it seems to be an outlier.')
            else:
                data_trn.append(frctrnmin)
                data_tst.append(frctstmin)

    trnarray = np.array(data_trn)
    tstarray = np.array(data_tst)
    print('Mean, std, max, and min:')
    dmean = trnarray.mean()
    dstd = trnarray.std()
    dmax = trnarray.max()
    dmin = trnarray.min()
    print('Training= {0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f}'.format(dmean,dstd,
                                                                     dmax,dmin))
    dmean = tstarray.mean()
    dstd  = tstarray.std()
    dmax  = tstarray.max()
    dmin  = tstarray.min()
    print('Test=     {0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f}'.format(dmean,dstd,
                                                                     dmax,dmin))
    
