#!/usr/bin/env python
"""
Calculate elastic constants C_ij by least square fitting
using 24 Green-Lagrange deformation tensors, which is the same
scheme as materialsproject.org does.

PREPARE mode creates directories and POSCARs to be computed.
ANALYZE mode reads the stress values obtained by some code.

Usage:
    elastic_constants.py prepare [options] POSCAR
    elastic_constants.py analyze STRSFILE

Options:
  -h, --help  Shows this message and exit.
"""
from __future__ import print_function

import os
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit

from ase.io import read

__author__  = 'Ryo KOBAYASHI'
__version__ = '170308'
__licence__ = 'MIT'

_confname = 'conf.elastic_constants.json'

#...constants
_outfname = 'out.elasticity'
_prefix = 'elasticity_'
_dlt1s = (-0.01, -0.005, 0.005, 0.01)
_dlt2s = (-0.06, -0.03, 0.03, 0.06)

def quad_func(x,a,b):
    return a *x**2 +b

def get_erg(fname):
    with open(fname) as f:
        l = f.readline()
        erg = float(l.split()[0])
    return erg

def get_deformations():

    fmats = []
    for dlt1 in _dlt1s:
        for t in range(3):
            fmat = np.identity(3,dtype=float)
            fmat[t,t] += dlt1
            fmats.append(fmat)

    for dlt2 in _dlt2s:
        for t in range(3):
            fmat = np.identity(3,dtype=float)
            if t==0:
                fmat[1,2] = dlt2
            elif t==1:
                fmat[0,2] = dlt2
            elif t==2:
                fmat[0,1] = dlt2
            fmats.append(fmat)
    return fmats

def prepare(infname='POSCAR'):

    #...original system
    orig_atoms = read(infname,format='vasp')

    #...get 24 deformations
    fmats = get_deformations()

    #...deform original system and save to _prefix_##/POSCAR
    for i,fmat in enumerate(fmats):
        atoms = orig_atoms.copy()
        dname = _prefix +"{0:02d}".format(i)
        os.system('mkdir -p {}'.format(dname))
        print(dname)
        cell0 = atoms.get_cell()
        cell = np.dot(cell0,fmat.T)
        atoms.set_cell(cell,scale_atoms=True)
        atoms.write(dname+'/POSCAR',format='vasp',vasp5=True,direct=True,
                    sort=False)
    
    print('prepare done')
    print('')
    print('Perform VASP calculations in these directories '
          +'and run the following command,')
    print('  python elasticity.py fit')

def cdote(strns,*params):
    """
    Compute C_ij*e_j to get s_i.
    """
    ctnsr = params2ctnsr(params)
    
    strss = np.zeros((len(strns),6),dtype=float)
    for i,strn in enumerate(strns):
        strs = np.dot(ctnsr,strn)
        strss[i] = strs
    return strss.flatten()


def params2ctnsr(params):
    """
    Create C_ij tensor from flat params vector assuring symmetry of C_ij.
    """
    ctnsr = np.zeros((6,6),dtype=float)
    # n = 0
    # for i in range(6):
    #     for j in range(i,6):
    #         ctnsr[i,j] = params[n]
    #         ctnsr[j,i] = ctnsr[i,j]
    #         n += 1
    ctnsr[0,0] = params[0]
    ctnsr[1,1] = params[1]
    ctnsr[2,2] = params[2]
    ctnsr[1,2] = params[3]
    ctnsr[0,2] = params[4]
    ctnsr[0,1] = params[5]
    ctnsr[3,3] = params[6]
    ctnsr[4,4] = params[7]
    ctnsr[5,5] = params[8]
    ctnsr[3,5] = params[9]
    ctnsr[0,4] = params[10]
    ctnsr[1,4] = params[11]
    ctnsr[2,4] = params[12]

    for i in range(6):
        for j in range(i,6):
            ctnsr[j,i] = ctnsr[i,j]
    return ctnsr


def analyze(strsfname):

    #...get 24 deformations
    fmats = get_deformations()
    strns = []
    for i,fmat in enumerate(fmats):
        emat = np.zeros((3,3),dtype=float)
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3))
        strn = np.zeros(6,dtype=float)
        strn[0] = emat[0,0]
        strn[1] = emat[1,1]
        strn[2] = emat[2,2]
        strn[3] = (emat[1,2]+emat[2,1])
        strn[4] = (emat[0,2]+emat[2,0])
        strn[5] = (emat[0,1]+emat[1,0])
        strns.append(strn)

    #...get 24 stress values
    strss = np.zeros((len(fmats),6),dtype=float)
    for i in range(len(fmats)):
        dname = _prefix +"{0:02d}".format(i)
        #print(dname)
        try:
            with open(dname+'/'+strsfname,'r') as f:
                data = f.readline().split()
                strss[i] = np.array([ float(d) for d in data ])
        except:
            raise
        print(fmats[i])
        print(strss[i])
    strss = strss.flatten()

    #...parameters 21 elements
    #params = np.zeros(21,dtype=float)
    #...parameters 13 elements
    params = np.zeros(13,dtype=float)

    #...fit
    opt,covar = curve_fit(cdote,strns,strss,p0=params)
    ctnsr = params2ctnsr(opt)

    print('C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')
    


if __name__ == '__main__':

    args= docopt(__doc__,version=__version__)

    if args['prepare']:
        infname = args['POSCAR']
        prepare(infname)
    elif args['analyze']:
        strsfname = args['STRSFILE']
        analyze(strsfname)
