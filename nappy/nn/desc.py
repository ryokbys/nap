#!/usr/bin/env python
"""
Packages for descriptors.

Usage:
  desc.py [options]

Options:
  -h, --help  Show this message and exit.
  -i IFNAME   Input file name. [default: in.params.desc]
  -o OFNAME   Output file name. [default: in.params.desc]
"""
from __future__ import print_function

from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "180807"

def read_desc(fname='in.params.desc'):
    with open(fname,'r') as f:
        lines = f.readlines()
    r_inner = {}
    for line in lines:
        data = line.split()
        if data[0] == '!' or data[0] == '#':
            if len(data) == 1:
                continue
            if len(data) == 5 and data[1] == 'inner_cutoff:':
                isp = int(data[2])
                rin = float(data[3])
                rout = float(data[4])
                r_inner[isp] = (rin,rout)
        elif len(data) == 2:
            nsp = int(data[0])
            nsf = int(data[1])
            descs = []
        elif len(data) > 4:
            itype = int(data[0])
            isp = int(data[1])
            jsp = int(data[2])
            if itype == 1: # Gaussian
                rc = float(data[3])
                eta = float(data[4])
                rs = float(data[5])
                descs.append(('gauss',isp,jsp,rc,eta,rs))
            elif itype == 101: # angular
                ksp = int(data[3])
                rc = float(data[4])
                almbd = float(data[5])
                descs.append(('angular',isp,jsp,ksp,rc,almbd))
    return nsp,nsf,descs,r_inner

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['-i']
    outfname = args['-o']

    nsp,nsf,descs,r_inner = read_desc(infname)

    
    
