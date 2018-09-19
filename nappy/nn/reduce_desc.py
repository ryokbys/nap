#!/usr/bin/env python
"""
Reduce the number of descriptors in `in.params.desc` by looking at 
`in.params.linreg` in case of linreg or `in.vars.fitpot` in case of NN2.

Usage:
  reduce_desc.py [options] IN_DESC IN_PARAMS

Options:
  -h, --help  Show this message and exit.
  -t TOL      Tolerance of weight value to be removed. [default: 1e-8]
  -n NUMS     Numbers of nodes in case of NN2, such as "10,5" for 10-5-1 network. [default: 10,5]
"""
from __future__ import print_function

from docopt import docopt
import nappy.nn.desc as desc

__author__ = "RYO KOBAYASHI"
__version__ = "180807"

def read_vars(fname='in.vars.fitpot'):
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

def write_vars(nsf1,rc2,rc3,weights,fname='in.vars.fitpot.new'):
    with open(fname,'w') as f:
        f.write(' {0:6d}  {1:8.3f}  {2:8.3f}\n'.format(nsf1,rc2,rc3))
        for w in weights:
            f.write('  {0:12.4e}  {1:12.4e}  {2:12.4e}\n'.format(w[0],w[1],w[2]))
    return

def get_reduced_descs_linreg(nsp,nsf,descs,weights,wtol):
    descs_new = []
    weights_new = []
    for i,d in enumerate(descs):
        w = weights[i][0]
        if abs(w) < wtol:
            continue
        descs_new.append(d)
        weights_new.append(weights[i])
    return descs_new,weights_new
        
def get_reduced_descs_nn(nsp,nums,descs,weights,wtol):
    descs_new = []
    weights_new = []
    inc = 0
    for i,d in enumerate(descs):
        w_has_value = False
        for j in range(nums[1]):
            w = weights[inc][0]
            inc += 1
            if abs(w) > wtol:
                w_has_value = True
        if w_has_value:
            descs_new.append(d)
            n = i*nums[1]
            for k in range(nums[1]):
                weights_new.append(weights[n+k])
    #...Add weights between 2nd and output layers
    n = nums[1]*nums[0]
    for i in range(nums[1]):
        weights_new.append(weights[n])
        n += 1
    return descs_new,weights_new
        

if __name__ == "__main__":

    args = docopt(__doc__)
    f_desc = args['IN_DESC']
    f_params = args['IN_PARAMS']
    of_desc = f_desc+'.new'
    of_params = f_params+'.new'
    wtol = float(args['-t'])
    nums = [ int(x) for x in args['-n'].split(',') ]
    nums.append(1)

    nsp,nsf,descs,r_inner = desc.read_desc(f_desc)
    nsf1,rc2,rc3,weights = read_vars(f_params)
    if 'linreg' in f_params:
        if nsf != nsf1:
            raise ValueError('nsf != nsf1, which should not happen in case of '+f_params)
        descs_new,weights_new = get_reduced_descs_linreg(nsp,nsf,descs,weights,wtol)
    elif 'vars' in f_params:
        if nsf != nums[0]:
            raise ValueError('nsf != nums[0], which should not happen in case of '+f_params)
        if len(nums) < 3:
            raise ValueError('nums option should be wrong.')
        descs_new,weights_new = get_reduced_descs_nn(nsp,nums,descs,weights,wtol)
        
    print(' Number of chosen descriptors  = {0:d}'.format(len(descs_new)))
    print(' Number of removed descriptors = {0:d}'.format(len(descs)-len(descs_new)))
    desc.write_desc(nsp,len(descs_new),descs_new,r_inner,fname=of_desc)
    write_vars(len(weights_new),rc2,rc3,weights_new,fname=of_params)
    
    
