#!/usr/bin/env python
"""
Reduce the number of descriptors in `in.params.desc` by looking at `in.params.linreg`.

Usage:
  reduce_desc.py [options] IN_DESC IN_LINREG

Options:
  -h, --help  Show this message and exit.
  -t TOL      Tolerance of weight value to be removed. [default: 1e-8]
"""
from __future__ import print_function

from docopt import docopt
import nappy.nn.desc as desc

__author__ = "RYO KOBAYASHI"
__version__ = "180807"

def read_linreg(fname='in.params.linreg'):
    with open(fname,'r') as f:
        lines = f.readlines()
    nsf = 0
    rc2 = 0.0
    rc3 = 0.0
    weights = []
    for line in lines:
        if line[0] == '!' or line[0] == '#':
            continue
        data = line.split()
        if nsf == 0:
            nsf = int(data[0])
            rc2 = float(data[1])
            rc3 = float(data[2])
        else:
            weights.append([ float(x) for x in data ])
    return nsf,rc2,rc3,weights

def write_linreg(nsf1,rc2,rc3,weights,fname='in.params.linreg.new'):
    with open(fname,'w') as f:
        f.write(' {0:6d}  {1:8.3f}  {2:8.3f}\n'.format(nsf1,rc2,rc3))
        for w in weights:
            f.write('  {0:12.4e}  {1:12.4e}  {2:12.4e}\n'.format(w[0],w[1],w[2]))
    return

def get_reduced_descs(nsp,nsf,descs,weights,wtol):
    descs_new = []
    weights_new = []
    for i,d in enumerate(descs):
        w = weights[i][0]
        if abs(w) < wtol:
            continue
        descs_new.append(d)
        weights_new.append(weights[i])
    return descs_new,weights_new

if __name__ == "__main__":

    args = docopt(__doc__)
    f_desc = args['IN_DESC']
    f_linreg = args['IN_LINREG']
    of_desc = f_desc+'.new'
    of_linreg = f_linreg+'.new'
    wtol = float(args['-t'])

    nsp,nsf,descs,r_inner = desc.read_desc(f_desc)
    nsf1,rc2,rc3,weights = read_linreg(f_linreg)
    if nsf1 != nsf:
        raise ValueError('nsf1 != nsf, which should not happen.')
    descs_new,weights_new = get_reduced_descs(nsp,nsf,descs,weights,wtol)
    print(' Number of chosen descriptors  = {0:d}'.format(len(descs_new)))
    print(' Number of removed descriptors = {0:d}'.format(len(descs)-len(descs_new)))
    desc.write_desc(nsp,len(descs_new),descs_new,r_inner,fname=of_desc)
    write_linreg(len(weights_new),rc2,rc3,weights_new,fname=of_linreg)
    
