#!/usr/bin/env python
"""
Calculate elastic constants C_ij by least square fitting
using Green-Lagrange deformation tensors, which is the same
scheme as materialsproject.org does.

PREPARE mode creates directories and POSCARs to be computed.
ANALYZE mode reads the stress values obtained by some code.

Usage:
    elasticity.py prepare [options] POSCAR
    elasticity.py analyze [options] STRSFILE

Options:
  -h, --help  Shows this message and exit.
  --delta1 DELTA1
              Maximum strain value of diagonal elements. [default: 0.002]
  --delta2 DELTA2
              Maximum strain value of off-diagonal elements. [default: 0.01]
"""
from __future__ import print_function

import os
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit, fmin_cg, check_grad

from ase.io import read

__author__  = 'Ryo KOBAYASHI'
__version__ = '170308'
__licence__ = 'MIT'

_confname = 'conf.elastic_constants.json'

#...constants
_outfname = 'out.elasticity'
_prefix = 'elasticity_'
# _dlt1s = (-0.01, -0.005, 0.005, 0.01)
# _dlt2s = (-0.06, -0.03, 0.03, 0.06)

def quad_func(x,a,b):
    return a *x**2 +b

def get_erg(fname):
    with open(fname) as f:
        l = f.readline()
        erg = float(l.split()[0])
    return erg

def get_deformations(dlt1max,dlt2max):

    dlt1s = [-dlt1max, -dlt1max/2, dlt1max/2, dlt1max]
    dlt2s = [-dlt2max, -dlt2max/2, dlt2max/2, dlt2max]

    fmats = []
    for dlt1 in dlt1s:
        for t in range(3):
            fmat = np.identity(3,dtype=float)
            fmat[t,t] += dlt1
            fmats.append(fmat)

    # elems = ((1,2),(0,2),(0,1),(2,1),(2,0),(1,0))
    elems = ((1,2),(0,2),(0,1))
    for dlt2 in dlt2s:
        for e in elems:
            fmat = np.identity(3,dtype=float)
            fmat[e] = dlt2
            fmats.append(fmat)
        # for t in range(3):
        #     fmat = np.identity(3,dtype=float)
        #     if t==0:
        #         fmat[1,2] = dlt2
        #         # fmat[2,1] = dlt2
        #     elif t==1:
        #         fmat[0,2] = dlt2
        #         # fmat[2,0] = dlt2
        #     elif t==2:
        #         fmat[0,1] = dlt2
        #         # fmat[1,0] = dlt2
        #     fmats.append(fmat)
    return fmats

def prepare(infname='POSCAR',dlt1max=0.01,dlt2max=0.06):

    #...original system
    orig_atoms = read(infname,format='vasp')

    #...get deformations
    fmats = get_deformations(dlt1max,dlt2max)

    #...deform original system and save to _prefix_##/POSCAR
    for i,fmat in enumerate(fmats):
        atoms = orig_atoms.copy()
        dname = _prefix +"{0:02d}".format(i)
        os.system('mkdir -p {}'.format(dname))
        print(dname)
        cell0 = atoms.get_cell()
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3)) +np.identity(3)
        cell = np.dot(emat,cell0)
        # cell = np.dot(cell0,fmat.T)
        atoms.set_cell(cell,scale_atoms=True)
        atoms.write(dname+'/POSCAR',format='vasp',vasp5=True,direct=True,
                    sort=False)
    
    print('prepare done')
    print('')
    print('Perform VASP calculations in these directories '
          +'and run the following command,')
    print('  python elasticity.py analyze str.ref')

def cdote(strns,*params):
    """
    Compute C_ij*e_j to get s_i.
    """
    ctnsr = params2ctnsr(params)
    
    strss = np.zeros((len(strns),6),dtype=float)
    for i,strn in enumerate(strns):
        strs = np.dot(ctnsr,strn)
        #print('i,strs= ',i,strs)
        strss[i] = strs
    return strss.flatten()

def func(x,*args):
    """
    Objective function to be minimized in least square fitting.

    .. math::

      L= \sum_s 1/N_s \sum_i 1/6 [ S_i^s -\sum_j C_{ij} e_j^s]^2

    
    Options
    -------
    x : N-dimensional array
        Variable array to be optimized.
    args : 2-dimensional array of parameter arrays
        args[0] is the array of arrays of strains given by this script.
        args[1] is the array of arrays of stresses obtained by the external program.
        These arrays should not be flattened.

    Returns
    -------
    val : float
        Objective function value (scalar).
    """
    ctnsr = params2ctnsr(x)
    # print(args)
    strns = args[0]
    strss0 = args[1]
    val = 0.0
    n = 0
    for i,strn in enumerate(strns):
        strs = np.dot(ctnsr,strn)
        strs0 = strss0[i]
        for j in range(len(strs)):
            val += (strs[j]-strs0[j])**2
            n += 1
    val /= n
    print('val = ',val)
    return val

def dfunc(x,*args):
    """
    Derivative of the objective function to be minimized in least square fitting
    with respect to the parameters.
    
    Options
    -------
    x : N-dimensional float
        Variable array to be optimized.
    args : 2-dimensional array of parameter arrays
        args[0] is the array of arrays of strains given by this script.
        args[1] is the array of arrays of stresses obtained by the external program.
        These arrays should not be flattened.

    Returns
    -------
    df : N-dimensional float
        Derivative of the objective function value.
    """
    ctnsr = params2ctnsr(x)
    strns = args[0]
    strss0 = args[1]
    print('ctnsr=',ctnsr)
    
    # residue vector, (S_i -\sum_j C_{ij} e_j)
    residue = np.zeros(6,dtype=float)
    for i,strn in enumerate(strns):
        strs = np.dot(ctnsr,strn)
        # print('i=',i)
        # print('strn =',strn)
        # print('strs =',strs)
        strs0= strss0[i]
        # print('strs0=',strs0)
        for j in range(len(strs)):
            residue[j] += strs0[j] -strs[j]
        # print('i,strs0-strs =',i,strs0-strs)
        
    print('residue = ',residue)
    df = np.zeros(len(x),dtype=float)
    n = 0
    for i in range(6):
        for j in range(i,6):
            if i == j:
                for istrn,strn in enumerate(strns):
                    df[n] += 2.0*(-strn[j])*residue[i]
            else:
                for istrn,strn in enumerate(strns):
                    df[n] += 2.0*(-strn[j])*residue[i] \
                             +2.0*(-strn[i])*residue[j]
            n += 1
    # for i in range(len(df)):
    #     print('i,df[i] = ',i,df[i])
    return df
    
def params2ctnsr(params):
    """
    Create C_ij tensor from flat params vector assuring symmetry of C_ij.
    """
    ctnsr = np.zeros((6,6),dtype=float)
    n = 0
    for i in range(6):
        for j in range(i,6):
            ctnsr[i,j] = params[n]
            n += 1
    # ctnsr[0,0] = params[0]
    # ctnsr[1,1] = params[1]
    # ctnsr[2,2] = params[2]
    # ctnsr[1,2] = params[3]
    # ctnsr[0,2] = params[4]
    # ctnsr[0,1] = params[5]
    # ctnsr[3,3] = params[6]
    # ctnsr[4,4] = params[7]
    # ctnsr[5,5] = params[8]
    # ctnsr[3,5] = params[9]
    # ctnsr[0,4] = params[10]
    # ctnsr[1,4] = params[11]
    # ctnsr[2,4] = params[12]
    for i in range(6-1):
        for j in range(i+1,6):
            ctnsr[j,i] = ctnsr[i,j]
    return ctnsr


def analyze(strsfname,dlt1max=0.01,dlt2max=0.06):

    #...get deformations
    fmats = get_deformations(dlt1max,dlt2max)
    strns = []
    for i,fmat in enumerate(fmats):
        emat = np.zeros((3,3),dtype=float)
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3))
        # print('fmat.T*fmat=',np.dot(fmat.T,fmat))
        # print('identity',np.identity(3))
        # print('fmats=',fmats[i])
        # print('emat=',emat)
        strn = np.zeros(6,dtype=float)
        strn[0] = emat[0,0]
        strn[1] = emat[1,1]
        strn[2] = emat[2,2]
        strn[3] = emat[1,2] *2.0
        strn[4] = emat[0,2] *2.0
        strn[5] = emat[0,1] *2.0
        # print('i,strn=',i,strn)
        strns.append(strn)

    #...get stress values from external calculations
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
        # print('strs=',strss[i])
        # print(strns[i])

    #...parameters 21 elements
    params = np.zeros(21,dtype=float)
    #...parameters 13 elements
    #params = np.zeros(13,dtype=float)

    # params[0] = 510.0
    # params[1] = 200.0
    # params[2] = 200.0
    # params[6] = 510.0
    # params[7] = 200.0
    # params[11]= 510.0
    # params[15]= 150.0
    # params[18]= 150.0
    # params[20]= 150.0

    # print(check_grad(func,dfunc, params, strns,strss, epsilon=1.0))
    # return

    #....CG
    # opt = fmin_cg(func,params,fprime=dfunc,args=(strns,strss),full_output=False)
    
    #...fit
    strss = strss.flatten()
    opt,covar = curve_fit(cdote,strns,strss,p0=params)

    ctnsr = params2ctnsr(opt)
    # perr = np.sqrt(np.diag(covar))
    # print('std dev = ',perr)

    print('C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')

    # strss2 = np.zeros((len(fmats),6),dtype=float)
    # for i,strn in enumerate(strns):
    #     strs = np.dot(ctnsr,strn)
    #     # print('i,strs = {0:3d}'.format(i),end='')
    #     # for s in strs:
    #     #     print(' {0:8.2f}'.format(s),end='')
    #     # print('')
    #     strss2[i,:] = strs[:]
    # strss2 = strss2.flatten()
    # #strss = strss.flatten()
    # n= 0
    # for i in range(len(strss)):
    #     n += 1
    #     print('i,strs2,strs,dstrs = {0:4d}'.format(i)
    #           +' {0:10.4f} {1:10.4f} {2:10.4f}'.format(strss2[i],
    #                                                    strss[i],
    #                                                    abs(strss2[i]-strss[i])))
    #     if n== 6:
    #         print('')
    #         n = 0


if __name__ == '__main__':

    args= docopt(__doc__,version=__version__)
    
    dlt1max = float(args['--delta1'])
    dlt2max = float(args['--delta2'])

    if args['prepare']:
        infname = args['POSCAR']
        prepare(infname,dlt1max=dlt1max,dlt2max=dlt2max)
    elif args['analyze']:
        strsfname = args['STRSFILE']
        analyze(strsfname,dlt1max=dlt1max,dlt2max=dlt2max)
