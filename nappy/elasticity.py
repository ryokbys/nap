#!/usr/bin/env python
"""
Calculate elastic constants C_ij by least square fitting
using Green-Lagrange deformation tensors, which is the same
scheme as materialsproject.org does.

PREPARE mode creates directories and input files (as the same name as INFILE) to be computed. INFILE can be POSCAR or pmdini.
ANALYZE mode reads the stress values obtained by some code. 
Stress information is read from pmdfin file.
COMPUTE mode performs MD or MS calculations using pmd as a backend.

Usage:
    elasticity.py prepare [options] INFILE
    elasticity.py analyze [options] INFILE

Options:
  -h, --help  Shows this message and exit.
  --delta1 DELTA1
              Maximum strain value of diagonal elements. Available only with prepare. [default: 0.002]
  --delta2 DELTA2
              Maximum strain value of off-diagonal elements. Available only with prepare. [default: 0.01]
  --out4fp    Whether or not output Cij data for fp.py.
  --out4fp-name FNAME
              Filename of the output of Cij for fp.py. [default: data.pmd.Cij]
"""
import os,sys
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit, least_squares
import spglib
import yaml

import ase.io
import nappy
from nappy.napsys import NAPSystem
import copy

__author__  = 'Ryo KOBAYASHI'
__version__ = '230107'
__licence__ = 'MIT'

_confname = 'conf.elast.yaml'

#...constants
_prefix = 'elast_'

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
    #nsys0 = NAPSystem(fname=infname)
    nsys0 = nappy.io.read(infname)
    #orig_atoms = read(infname,format='vasp')
    orig_atoms = nsys0.to_ase_atoms()

    #...get deformations
    fmats = get_deformations(dlt1max,dlt2max)

    #...deform original system and save to _prefix_##/POSCAR
    for i,fmat in enumerate(fmats):
        atoms = orig_atoms.copy()
        #nsys = copy.deepcopy(nsys0)
        dname = _prefix +"{0:02d}".format(i)
        os.system('mkdir -p {}'.format(dname))
        print(dname)
        cell0 = atoms.get_cell()
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3)) +np.identity(3)
        cell = np.dot(emat,cell0)
        #print(i,emat,cell)
        # cell = np.dot(cell0,fmat.T)
        atoms.set_cell(cell,scale_atoms=True)
        if 'POSCAR' in infname:
            atoms.write(dname+'/'+infname,format='vasp',vasp5=True,direct=True,
                        sort=False)
        elif 'pmdini' in infname:
            nsys = nappy.io.from_ase(atoms)
            nappy.io.write(nsys,fname=dname+'/'+infname)
        else:
            raise ValueError('Available INFILEs: POSCAR and pmdini')
        with open(dname+'/README','w') as f:
            f.write('# Strain of the sample made by elasticity.py.\n\n')
            f.write('strain tensor:\n' )
            f.write('  {0:6.3f}  {1:6.3f}  {2:6.3f}\n'.format(*emat[0,0:3]))
            f.write('  {0:6.3f}  {1:6.3f}  {2:6.3f}\n'.format(*emat[1,0:3]))
            f.write('  {0:6.3f}  {1:6.3f}  {2:6.3f}\n'.format(*emat[2,0:3]))
            f.write('\nResulting cell matrix:\n')
            f.write('  {0:10.3f}  {1:10.3f}  {2:10.3f}\n'.format(*cell[0,0:3]))
            f.write('  {0:10.3f}  {1:10.3f}  {2:10.3f}\n'.format(*cell[1,0:3]))
            f.write('  {0:10.3f}  {1:10.3f}  {2:10.3f}\n'.format(*cell[2,0:3]))
    
    print('prepare done')
    print('')
    print('After performing VASP or pmd calculations in these directories, '
          +'run the following command:')
    print('  $ python elasticity.py analyze str.ref')
    print('or')
    print('  $ python elasticity.py analyze strs.pmd')
    print('')
    return None

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
            strsj = strs[j]
            strsj0= strs0[j]
            val += (strs[j]-strs0[j])**2
            n += 1
    val /= n
    # print('val = ',val)
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

def analyze(infname,dlt1max=0.01,dlt2max=0.06):

    #...original system
    nsys0 = nappy.io.read(infname)
    #atoms0 = ase.io.read(infname,format='vasp')
    atoms0 = nsys0.to_ase_atoms()
    #strs0 = nsys0.get_stress()

    #...get deformations
    fmats = get_deformations(dlt1max,dlt2max)
    print('\n F-matrices:')
    for i in range(len(fmats)):
        print('   {0:>3d}:'.format(i),fmats[i].flatten())
    strns = np.zeros((len(fmats),6),dtype=float)
    for i,fmat in enumerate(fmats):
        emat = np.zeros((3,3),dtype=float)
        emat = 0.5 *(np.dot(fmat.T,fmat) -np.identity(3))
        # strn = np.zeros(6,dtype=float)
        # strn[0] = emat[0,0]
        # strn[1] = emat[1,1]
        # strn[2] = emat[2,2]
        # strn[3] = emat[1,2] *2.0
        # strn[4] = emat[0,2] *2.0
        # strn[5] = emat[0,1] *2.0
        # strns.append(strn)
        strns[i,0] = emat[0,0]
        strns[i,1] = emat[1,1]
        strns[i,2] = emat[2,2]
        strns[i,3] = emat[1,2] *2.0
        strns[i,4] = emat[0,2] *2.0
        strns[i,5] = emat[0,1] *2.0
    print('\n strains:')
    for i in range(len(fmats)):
        print('   {0:>3d}:'.format(i),strns[i])
    
    #...get stress values from external calculations
    print('\n stresses:')
    strss = np.zeros((len(fmats),6),dtype=float)
    for i in range(len(fmats)):
        dname = _prefix +"{0:02d}".format(i)
        try:
            # with open(dname+'/'+strsfname,'r') as f:
            #     data = f.readline().split()
            #     strss[i] = np.array([ float(d) for d in data ])
            #     print('   {0:>3d}:'.format(i),strss[i])
            nsysi = nappy.io.read(dname+'/pmdfin')
            strsi = nsysi.get_stress()
            # Need to be careful about the definition of sign of stress value.
            # The tensile stress should be positive in order to get appropriate signs of Cij.
            # So multiplying minus sign here, because the sign of stress in pmdfin is opposite.
            strss[i] = -strsi
        except Exception as e:
            raise

    #...parameters 21 elements
    params = np.zeros(21,dtype=float)
    #...parameters 13 elements
    #params = np.zeros(13,dtype=float)
    
    #...optimize using curve_fit method, which may result in segmentation fault (230423)
    # strss = strss.flatten()
    # print('curve_fit...')
    # opt,covar = curve_fit(cdote,strns,strss,p0=params)
    
    #...Optimize using least_squares method, which would be more stable than curve_fit?
    args = [ strns, strss ]
    res = least_squares(func, params, args=args)
    opt = res.x
    
    print('params2ctnsr')
    ctnsr = params2ctnsr(opt)
    # perr = np.sqrt(np.diag(covar))
    # print('std dev = ',perr)

    print(' delta1_max = {0:8.5f}'.format(dlt1max))
    print(' delta2_max = {0:8.5f}'.format(dlt2max))
    print('')

    print(' C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')

    cij = reduce_cij(atoms0,ctnsr)
    print(' C_ij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:10.3f}'.format(ctnsr[i,j]),end='')
        print('')
    some_moduli(cij)
    return cij

def some_moduli(cij):
    sij = np.linalg.inv(cij)
    c123 = cij[0,0]+cij[1,1]+cij[2,2]
    c231 = cij[0,1]+cij[0,2]+cij[1,2]
    c456 = cij[3,3]+cij[4,4]+cij[5,5]
    s123 = sij[0,0]+sij[1,1]+sij[2,2]
    s231 = sij[0,1]+sij[0,2]+sij[1,2]
    s456 = sij[3,3]+sij[4,4]+sij[5,5]
    kv = (c123 +2.0*c231)/9
    kr = 1.0/(s123 +2.0*(s231))
    gv = (c123 -c231 +3.0*c456)/15
    gr = 15.0 /(4.0*s123 -4.0*s231 +3.0*s456)
    kvrh = (kv+kr)/2
    gvrh = (gv+gr)/2
    prto2 = (3.0*kvrh -2.0*gvrh)/(6.0*kvrh +2.0*gvrh)

    print('')
    print(f' Bulk modulus (K)     = {kvrh:10.3f} GPa')
    print(f' shear modulus (G)    = {gvrh:10.3f} GPa')
    print(f" Poisson's ratio (nu) = {prto2:10.3f}")
    print('')
    print(' Definition of elastic moduli, see ' \
        +'https://materialsproject.org/wiki/index.php/Elasticity_calculations')
    txt = """   c123 = c11 +c22 +c33  = {0:10.3f}
   c231 = c12 +c13 +c23  = {1:10.3f}
   c456 = c44 +c55 +c66  = {2:10.3f}
   s123 = s11 +s22 +s33  = {3:10.3f}
   s231 = s12 +s13 +s23  = {4:10.3f}
   s456 = s44 +s55 +s66  = {5:10.3f}
   Kv   = (c123 +2*c231)/9  = {6:10.3f}
   Kr   = 1.0 /(s123 +2*s231)  = {7:10.3f}
   Gv   = (c123 -c231 +3*c456)/15  = {8:10.3f}
   Gr   = 15.0 /(4.0*s123 -4.0*s231 +3.0*s456)  = {9:10.3f}
   K    = (Kv +Kr)/2
   G    = (Gv +Gr)/2
   nu   = (3*K -2*G)/(6*K +2*G)
""".format(c123,c231,c456,s123,s231,s456,kv,kr,gv,gr)
    print(txt)


def reduce_cij(atoms0,cij0,eps=1.e-4):
    """
    Reduce number of independent Cij according to the crystal system of original cell.
    It is not Cij=Cji.
    """

    cij = cij0
    symdata = spglib.get_symmetry_dataset(atoms0)
    #nsys = NAPSystem(ase_atoms=atoms0)
    nsys = nappy.io.from_ase(atoms0)
    sgnum = symdata['number']
    
    a,b,c = nsys.get_lattice_lengths()
    alpha,beta,gamma = nsys.get_lattice_angles()

    aeqb = abs(a-b) < eps*min(a,b)
    beqc = abs(b-c) < eps*min(b,c)
    ceqa = abs(c-a) < eps*min(c,a)
    
    aleqpi2 = abs(alpha-np.pi/2) < eps*np.pi/2
    bteqpi2 = abs(beta -np.pi/2) < eps*np.pi/2
    gmeqpi2 = abs(gamma-np.pi/2) < eps*np.pi/2

    print('Spacegroup number = ',sgnum,' ==> ',end='')
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

def out4fp_Cij(cij,fname='data.pmd.Cij'):
    """
    Output Cij data for fp.py.

    The order of the data is as following,
      11, 22, 33, 23, 13, 12, 44, 55, 66
    """
    from datetime import datetime
    with open(fname,'w') as f:
        nowstr = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        f.write('# Cij data computed by elasticity.py at {0:s}.\n'.format(nowstr))
        f.write('# The order of Cij is (11,22,33,23,13,12,44,55,66)\n')
        f.write('# datatype: independent\n')
        f.write('   9  1.0\n')
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[0,0]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[1,1]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[2,2]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[1,2]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[0,2]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[0,1]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[3,3]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[4,4]))
        f.write('   {0:12.4e}   {0:.1f}\n'.format(cij[5,5]))
    return None

def main():

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    # args= docopt(__doc__,version=__version__)
    dlt1max = float(args['--delta1'])
    dlt2max = float(args['--delta2'])
    out4fp = args['--out4fp']
    out4fp_name = args['--out4fp-name']

    infname = args['INFILE']
    
    if args['prepare']:
        prepare(infname,dlt1max=dlt1max,dlt2max=dlt2max)
        # overwrite conf.elast.yaml file
        conf = {'delta1':dlt1max, 'delta2':dlt2max}
        with open(_confname,'w') as f:
            f.write(yaml.dump(conf,default_flow_style=False))
    elif args['analyze']:
        try:
            with open(_confname,'r') as f:
                conf = yaml.safe_load(f)
        except Exception as e:
            raise
        print(' Loading ',_confname,' done\n   ==> ',conf)
        dlt1max = conf['delta1']
        dlt2max = conf['delta2']
        #strsfname = args['STRSFILE']
        cij = analyze(infname,dlt1max=dlt1max,dlt2max=dlt2max)
        #...Output Cij for fp.py
        if out4fp:
            out4fp_Cij(cij,fname=out4fp_name)
            print('\n Wrote ',out4fp_name)
        print('\n {0:s} analyze done'.format(os.path.basename(sys.argv[0])))

    return None

if __name__ == '__main__':

    main()
