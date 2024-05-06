#!/usr/bin/env python
"""
Extract necessary BVS parameters from a file that includes a lot of parameters.

Usage:
  extract_bvs_params.py [options] ALL_PARAMS_FILE IDS [IDS...]

Options:
  -h, --help   Show this message and exit.
"""
from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "170327"

out_Morse = 'in.params.Morse'
out_Coulomb= 'in.params.Coulomb'

def read_params_file(fname):
    params = {}
    with open(fname,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '#':
                continue
            data = line.split()
            idx = int(data[0])
            name = data[1]
            vid = int(data[2])
            d0 = float(data[3])
            rmin = float(data[4])
            alpha = float(data[5])
            crad = float(data[6])
            qnum = int(data[7])
            params[idx] = [name,vid,d0,rmin,alpha,crad,qnum]
    return params

def write_Morse_params(fname,params,id_order):
    with open(fname,'w') as f:
        n = 1
        for k in id_order:
            p = params[k]
            n += 1
            f.write('{0:4d} {1:4d}'.format(1,n)
                    +' {0:8.5f} {1:8.5f} {2:8.5f}\n'.format(p[2],p[4],p[3]))


def write_screened_Coulomb_params(fname,params,id_order):
    with open(fname,'w') as f:
        #...declare it is 'screened' Coulomb
        f.write(' screened_bvs \n')
        n = 1
        #...species #1 is oxygen in case of BVS_Adams
        f.write('{0:4d} {1:5s}'.format(n,'O')
                +' {0:4d} {1:6.2f} {2:4d}\n'.format(2, 0.66, 2))
        for k in id_order:
            n += 1
            p = params[k]
            f.write('{0:4d} {1:5s}'.format(n,p[0])
                    +' {0:4d} {1:6.2f} {2:4d}\n'.format(p[1], p[5], p[6]))
        #...Write interactions
        f.write(' interactions \n')
        f.write(' {0:3d} {1:3d}\n'.format(1,1))
        nsp = len(id_order)
        for i in range(2,2+nsp):
            for j in range(i,2+nsp):
                f.write(' {0:3d} {1:3d}\n'.format(i,j))
            

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['ALL_PARAMS_FILE']
    id_order = [ int(i) for i in args['IDS'] ]

    params = read_params_file(fname)
    write_Morse_params(out_Morse,params,id_order)
    write_screened_Coulomb_params(out_Coulomb,params,id_order)
