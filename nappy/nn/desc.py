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
            elif itype == 2: # cosine
                rc = float(data[3])
                xi = float(data[4])
                descs.append(('cosine',isp,jsp,rc,xi))
            elif itype == 101: # angular
                ksp = int(data[3])
                rc = float(data[4])
                almbd = float(data[5])
                descs.append(('angular',isp,jsp,ksp,rc,almbd))
            elif itype == 102: # angular2
                ksp = int(data[3])
                rc = float(data[4])
                almbd = float(data[5])
                descs.append(('angular2',isp,jsp,ksp,rc,almbd))
            elif itype == 103: # angular3
                ksp = int(data[3])
                rc = float(data[4])
                a1 = float(data[5])
                descs.append(('angular3',isp,jsp,ksp,rc,a1))
            elif itype == 104: # angular4
                ksp = int(data[3])
                rc = float(data[4])
                a1 = float(data[5])
                descs.append(('angular4',isp,jsp,ksp,rc,a1))
    return nsp,nsf,descs,r_inner

def write_desc(nsp,nsf,descs,r_inner,fname='in.params.desc'):
    with open(fname,'w') as f:
        if len(r_inner) != 0:
            for k,v in r_inner.items():
                isp = k
                rin = v[0]
                rout = v[1]
                f.write('#  inner_cutoff:  {0:d} {1:6.2f} {2:6.2f}\n'.format(isp,
                                                                             rin,
                                                                             rout))
        f.write(' {0:3d}  {1:5d}\n'.format(nsp,nsf))
        for d in descs:
            sftype = d[0]
            isp = d[1]
            jsp = d[2]
            if sftype == 'gauss':
                rc = d[3]
                eta = d[4]
                rs = d[5]
                f.write('    1 {0:2d} {1:2d}   '.format(isp,jsp))
                f.write('  {0:6.2f}  {1:9.4f} {2:8.4f}\n'.format(rc,eta,rs))
            elif sftype == 'cosine':
                rc = d[3]
                xi = d[4]
                f.write('    2 {0:2d} {1:2d}   '.format(isp,jsp))
                f.write('  {0:6.2f}  {1:9.4f}\n'.format(rc,xi))
            elif sftype == 'angular':
                ksp = d[3]
                rc = d[4]
                almbd = d[5]
                f.write('  101 {0:2d} {1:2d} {2:2d}'.format(isp,jsp,ksp))
                f.write('  {0:6.2f}  {1:9.4f}\n'.format(rc,almbd))
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['-i']
    outfname = args['-o']

    nsp,nsf,descs,r_inner = read_desc(infname)

    
    
