#!/usr/bin/env python
"""
Analyze energies and forces (stresses?) of specified samples.

Usage:
  analyze_samples.py [options] DIRS [DIRS...]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "170830"


if __name__ == "__main__":

    args = docopt(__doc__)

    dirs = args['DIRS']

    ergs = []
    frcs = []
    for d in dirs:
        with open(d+'/frc.ref','r') as g:
            lines = g.readlines()
            natm = int(lines[0].split()[0])
            for i in range(1,natm+1):
                data = lines[i].split()
                frc = [ float(x) for x in data[0:3] ]
        with open(d+'/erg.ref','r') as f:
            erg = float(f.readline().split()[0])
            erg /= natm
        ergs.append(erg)
        frcs.append(frc)

    e_ave = 0.0
    e_var = 0.0
    e_min = 1e+30
    e_max = -1e+30
    f_ave = 0.0
    f_var = 0.0
    f_min = 1e+30
    f_max = -1e+30
    for i in range(len(ergs)):
        erg = ergs[i]
        frc = frcs[i]
        e_ave += erg
        e_var += erg*erg
        e_min = min(e_min,erg)
        e_max = max(e_max,erg)
        for j in range(3):
            f_ave += frc[j]
            f_var += frc[j]*frc[j]
            f_min = min(f_min,frc[j])
            f_max = max(f_max,frc[j])
    e_ave /= len(ergs)
    e_var = e_var/len(ergs) -e_ave**2
    f_ave /= len(frcs)*3
    f_var = f_var/(len(frcs)*3) -f_ave**2

    print('Energy per atom:')
    print('  Average:            {0:8.4f}'.format(e_ave))
    print('  Standard deviation: {0:8.4f}'.format(np.sqrt(e_var)))
    print('  Minimum:            {0:8.4f}'.format(e_min))
    print('  Maximum:            {0:8.4f}'.format(e_max))
    print('  Energy range:       {0:8.4f}'.format(e_max-e_min))

    print('Force component:')
    print('  Average:            {0:8.4f}'.format(f_ave))
    print('  Standard deviation: {0:8.4f}'.format(np.sqrt(f_var)))
    print('  Minimum:            {0:8.4f}'.format(f_min))
    print('  Maximum:            {0:8.4f}'.format(f_max))
    print('  Froce range:        {0:8.4f}'.format(f_max-f_min))
