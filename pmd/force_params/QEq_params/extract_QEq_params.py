#!/usr/bin/env python
"""
Extract atomic parameters for QEq potential.

Usage:
  extract_bvs_params.py [options] DATA_FILE NAME [NAME...]

Options:
  -h, --help   Show this message and exit.
"""
from __future__ import print_function

from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "180112"

out_Coulomb= 'in.params.Coulomb'

def read_data_file(fname):
    params = {}
    with open(fname,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '#':
                continue
            data = line.split()
            idx = int(data[0])
            name = data[1]
            ie1 = float(data[2])
            ie2 = float(data[3])
            ea = float(data[4])
            rad = float(data[5])
            en = float(data[6])
            params[name] = [idx,name,ie1,ie2,ea,rad,en]
    return params

def anum_to_range(atomic_number):
    """
    Calculate and return the lower and upper limits of the charge of given atomic number.
    """
    nstates = (0,2,8,18,32,50)
    if atomic_number > 86:
        raise ValueError('Atomic number greater than 86 is not available.')
    elif atomic_number <= sum_array(nstates,1):
        n = 1
    elif atomic_number <= sum_array(nstates,2):
        n = 2
    elif atomic_number <= sum_array(nstates,3):
        n = 3
    elif atomic_number <= sum_array(nstates,4):
        n = 4
    elif atomic_number <= sum_array(nstates,5):
        n = 5
    else:
        raise ValueError('Atomic number is something wrong: ',atomic_number)

    freedom = (0,2,6,10,14,18,22)
    nval = atomic_number - sum_array(nstates,n-1)
    if nval < sum_array(freedom,1):
        l = 1
    elif nval < sum_array(freedom,2):
        l = 2
    elif nval < sum_array(freedom,3):
        l = 3
    elif nval < sum_array(freedom,4):
        l = 4
    else:
        l = 5
    if not l <= n:
        raise ValueError('not l<=n')

    print('anum,n,l,nval=',atomic_number,n,l,nval)
    nseat = sum_array(nstates,n) -sum_array(nstates,n-1)
    nseatopen = nseat - nval
    for il in range(l+1,n+1):
        nseatopen -= freedom[il]
        
    print('nseat,nseatopen=',nseat,nseatopen)
    qlow = -float(nseatopen)
    qup = float(min(nval, freedom[l]+freedom[l-1]))
    return qlow,qup

def sum_array(array,n):
    if len(array) < n+1:
        raise ValueError('len(array) < n')
    s = 0
    for i in range(n+1):
        s += array[i]
    return s

def write_Coulomb_params(fname,params,specorder):
    with open(fname,'w') as f:
        #...declare it is 'variable_charge' Coulomb
        f.write(' variable_charge \n')
        n = 0
        e0 = 0.0
        for k in specorder:
            n += 1
            p = params[k]
            anum = p[0]
            name = p[1]
            ie = p[2]
            ea = -p[4]
            xi = (ie+ea)/2
            J  = (ie-ea)
            qlow,qup = anum_to_range(anum)
            f.write('{0:4d} {1:5s}'.format(n,name)
                    +' {0:9.3f} {1:9.3f} {2:9.4f}'.format(xi,J,e0)
                    +' {0:5.1f} {1:5.1f}\n'.format(qlow,qup))
            

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['DATA_FILE']
    specorder = [ name for name in args['NAME'] ]

    params = read_data_file(fname)
    write_Coulomb_params(out_Coulomb,params,specorder)
