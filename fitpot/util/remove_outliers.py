#!/usr/bin/env python
"""
Remove/move outlier data within given data set from STDIN.

Usage:
  remove_outliers.py [options] DIRS [DIRS...]

Options:
  -h, --help    Show this message and exit.
  -t THRESHOLD  Threshold energy value (eV/atom) over which
                data will be removed. [default: 1.0]
"""
from __future__ import print_function

import os
from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "170117"

not_used_dir = 'not_used'

if __name__ == "__main__":

    args = docopt(__doc__)

    dirs = args['DIRS']
    threshold = float(args['-t'])

    os.system('mkdir -p '+not_used_dir)

    ergs = []
    for d in dirs:
        with open(d+'/erg.ref') as f:
            erg = float(f.readline())
        with open(d+'/frc.ref') as f:
            natm = int(f.readline())
        ergs.append(erg/natm)

    emin = min(ergs)
    print('emin = {0:12.7f}'.format(emin))
    for i,d in enumerate(dirs):
        e = ergs[i]
        ediff = e -emin
        if ediff > threshold:
            cmd = 'mv '+d+' '+not_used_dir+'/'
            print(cmd+'; because ediff = {0:12.7f} > threshold'.format(ediff))
            os.system(cmd)
