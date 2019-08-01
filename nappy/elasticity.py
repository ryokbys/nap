#!/usr/bin/env python
"""
Calculate elastic constants C_ij by least square fitting
using Green-Lagrange deformation tensors, which is the same
scheme as materialsproject.org does.

PREPARE mode creates directories and POSCARs to be computed.
ANALYZE mode reads the stress values obtained by some code.

Usage:
    elasticity.py prepare [options] POSCAR
    elasticity.py analyze [options] POSCAR STRSFILE

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
from scipy.optimize import curve_fit
import spglib

from ase.io import read
from nappy.napsys import NAPSystem

__author__  = 'Ryo KOBAYASHI'
__version__ = '170308'
__licence__ = 'MIT'

_confname = 'conf.elastic_constants.json'

#...constants
_outfname = 'out.elasticity'
_prefix = 'elasticity_'

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
            
    #elems = ((1,2),(0,2),(0,1),(2,1),(2,0),(1,0))
    elems = ((1,2),(0,2),(0,1))
    for dlt2 in dlt2s:
        for e in elems:
            fmat = np.identity(3,dtype=float)
            fmat[e] = dlt2
            fmats.append(fmat)
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
    print('After performing VASP or pmd calculations in these directories, '
          +'run the following command:')
    print('  $ python elasticity.py analyze str.ref')
    print('or')
    print('  $ python elasticity.py analyze strs.pmd')
    print('')

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
        strs0= strss0[i]
        for j in range(len(strs)):
            residue[j] += strs0[j] -strs[j]
        
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
    for i in range(6-1):
        for j in range(i+1,6):
            ctnsr[j,i] = ctnsr[i,j]
    return ctnsr


def analyze(infname,strsfname,dlt1max=0.01,dlt2max=0.06):

    #...original system
    atoms0 = read(infname,format='vasp')
    
    #...get deformations
    fmats = get_deformations(dlt1max,dlt2max)
    strns = []
    for i,fmat in enumerate(fmats):
        emat = np.zeros((3,3),dtype=float)
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3))
        strn = np.zeros(6,dtype=float)
        strn[0] = emat[0,0]
        strn[1] = emat[1,1]
        strn[2] = emat[2,2]
        strn[3] = emat[1,2] *2.0
        strn[4] = emat[0,2] *2.0
        strn[5] = emat[0,1] *2.0
        strns.append(strn)

    #...get stress values from external calculations
    strss = np.zeros((len(fmats),6),dtype=float)
    for i in range(len(fmats)):
        dname = _prefix +"{0:02d}".format(i)
        try:
            with open(dname+'/'+strsfname,'r') as f:
                data = f.readline().split()
                strss[i] = np.array([ float(d) for d in data ])
        except:
            raise

    #...parameters 21 elements
    params = np.zeros(21,dtype=float)
    #...parameters 13 elements
    #params = np.zeros(13,dtype=float)
    
    #...fit
    strss = strss.flatten()
    opt,covar = curve_fit(cdote,strns,strss,p0=params)

    ctnsr = params2ctnsr(opt)
    # perr = np.sqrt(np.diag(covar))
    # print('std dev = ',perr)

    print('delta1_max = {0:8.5f}'.format(dlt1max))
    print('delta2_max = {0:8.5f}'.format(dlt2max))
    print('')

    print('C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')

    cij = reduce_cij(atoms0,ctnsr)
    print('C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')
    some_moduli(cij)
    return cij

def some_moduli(cij):
    sij = np.linalg.inv(cij)
    c112233 = cij[0,0]+cij[1,1]+cij[2,2]
    c122331 = cij[0,1]+cij[0,2]+cij[1,2]
    c445566 = cij[3,3]+cij[4,4]+cij[5,5]
    s112233 = sij[0,0]+sij[1,1]+sij[2,2]
    s122331 = sij[0,1]+sij[0,2]+sij[1,2]
    s445566 = sij[3,3]+sij[4,4]+sij[5,5]
    kv = (c112233 +2.0*c122331)/9
    kr = 1.0/(s112233 +2.0*(s122331))
    gv = (c112233 -c122331 +3.0*c445566)/15
    gr = 15.0 /(4.0*s112233 -4.0*s122331 +3.0*s445566)
    kvrh = (kv+kr)/2
    gvrh = (gv+gr)/2
    prto2 = (3.0*kvrh -2.0*gvrh)/(6.0*kvrh +2.0*gvrh)

    print('')
    # print(' Definition of the following values, see ' \
    #     +'https://materialsproject.org/wiki/index.php/Elasticity_calculations')
    # print(' K_V   = {0:10.3f} GPa'.format(kv))
    # print(' K_R   = {0:10.3f} GPa'.format(kr))
    # print(' G_V   = {0:10.3f} GPa'.format(gr))
    # print(' G_R   = {0:10.3f} GPa'.format(gv))
    print(' Bulk modulus    = {0:10.3f} GPa'.format(kvrh))
    print(' shear modulus   = {0:10.3f} GPa'.format(gvrh))
    print(' Poisson\'s ratio = {0:10.3f}'.format(prto2))
    

def reduce_cij(atoms0,cij0,eps=1.e-4):
    """
    Reduce number of independent Cij according to the crystal system of original cell.
    It is not Cij=Cji.
    """

    cij = cij0
    symdata = spglib.get_symmetry_dataset(atoms0)
    napsys = NAPSystem(ase_atoms=atoms0)
    sgnum = symdata['number']
    print('Spacegroup number = ',sgnum,' ==> ',end='')
    
    a,b,c = napsys.get_lattice_lengths()
    alpha,beta,gamma = napsys.get_lattice_angles()

    aeqb = abs(a-b) < eps*min(a,b)
    beqc = abs(b-c) < eps*min(b,c)
    ceqa = abs(c-a) < eps*min(c,a)
    
    aleqpi2 = abs(alpha-np.pi/2) < eps*np.pi/2
    bteqpi2 = abs(beta -np.pi/2) < eps*np.pi/2
    gmeqpi2 = abs(gamma-np.pi/2) < eps*np.pi/2

    if 0 < sgnum <= 2:  # Triclinic
        print('Triclinic')
        pass
    elif sgnum <= 15:  # Monoclinic
        print('Monoclinic')
        pass
    elif sgnum <= 74:  # Orthorhombic
        print('Orthorhombic')
        pass
    elif sgnum <= 142:  # Tetragonal
        print('Tetragonal')
        pass
    elif sgnum <= 194:  # Hexagonal
        print('Hexagonal')
        print('Number of independent C_ij elements are reduced to 6.')
        print('C_66 should be 1/2(C_11-C_12) but this is not enforced now.')
        if not aleqpi2:
            c22 = (cij[1,1] +cij[2,2])/2
            c13 = (cij[0,1] +cij[0,2])/2
            c55 = (cij[4,4] +cij[5,5])/2
            cij[1,1] = cij[2,2] = c22
            cij[0,1] = cij[0,2] = c13
            cij[4,4] = cij[5,5] = c55
        elif not bteqpi2:
            c11 = (cij[0,0] +cij[2,2])/2
            c12 = (cij[0,1] +cij[1,2])/2
            c44 = (cij[3,3] +cij[5,5])/2
            cij[0,0] = cij[2,2] = c11
            cij[0,1] = cij[1,2] = c12
            cij[3,3] = cij[5,5] = c44
        elif not gmeqpi2:
            c11 = (cij[0,0] +cij[1,1])/2
            c12 = (cij[0,2] +cij[1,2])/2
            c44 = (cij[3,3] +cij[4,4])/2
            cij[0,0] = cij[1,1] = c11
            cij[0,2] = cij[1,2] = c12
            cij[3,3] = cij[4,4] = c44
        
    elif sgnum <= 230:  # Cubic
        print('Cubic')
        print('Number of independent C_ij elements are reduced to 3.')
        c11 = (cij[0,0] +cij[1,1] +cij[2,2])/3
        c12 = (cij[0,1] +cij[0,2] +cij[1,2])/3
        c44 = (cij[3,3] +cij[4,4] +cij[5,5])/3
        cij[0,0] = cij[1,1] = cij[2,2] = c11
        cij[0,1] = cij[0,2] = cij[1,2] = c12
        cij[3,3] = cij[4,4] = cij[5,5] = c44
    else:
        raise ValueError('Invalid space group number, ',sgnum)
    # Just symmetrize Cij
    for i in range(6):
        for j in range(i,6):
            cij[j,i] = cij[i,j]

    return cij

if __name__ == '__main__':

    args= docopt(__doc__,version=__version__)
    
    dlt1max = float(args['--delta1'])
    dlt2max = float(args['--delta2'])

    infname = args['POSCAR']
    
    if args['prepare']:
        prepare(infname,dlt1max=dlt1max,dlt2max=dlt2max)
    elif args['analyze']:
        strsfname = args['STRSFILE']
        analyze(infname,strsfname,dlt1max=dlt1max,dlt2max=dlt2max)
