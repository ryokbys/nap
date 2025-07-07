#!/usr/bin/env python
"""
Compute probability distribution of specified species.
The size of the system cell is determined as that of the system read 1st.

Usage:
  {0:s} [options] FILES [FILES...]

Options:
  -h, --help  Show this message and exit.
  -v          Verbose output for debugging.
  -w WIDTH    Mesh width in Ang. [default: 0.2]
  --sigma SIGMA
              Sigma value of the gaussian smearing. If it is positive,
              use this value in Ang. If it is negative, adopt |SIGMA|*WIDTH. [default: -1.0]
  --species SPC
              Secies name whose distributions are computed.
              If it is not given, all the atoms are to be taken into account. [default: None]
  --specorder SPECORDER
              Species order separated by comma. [default: None]
"""
import os
import sys
from docopt import docopt
import numpy as np
from scipy import stats

from nappy.common import get_key
#from nappy.napsys import NAPSystem
from nappy.io import read, write

from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "241126"

def get_num_division(nsys,width):
    ndivs = np.zeros(3,dtype=int)
    la,lb,lc = nsys.get_lattice_lengths()
    ndivs[0] = la/width +1
    ndivs[1] = lb/width +1
    ndivs[2] = lc/width +1
    return ndivs

def get_prob_dist_kde(ndivs, nsys, spc, sgm):
    tsgm2 = 2.0 * sgm * sgm
    tsgm2i = 1.0 / tsgm2
    la, lb, lc = nsys.get_lattice_lengths()
    ra = 1.0/ndivs[0]
    rb = 1.0/ndivs[1]
    rc = 1.0/ndivs[2]
    dv = la*ra *lb*rb *lc*rc
    prefactor = dv / (np.pi*tsgm2)**1.5
    poss = nsys.get_scaled_positions()
    symbols = np.array(nsys.get_symbols(), dtype=str)
    to_consider = symbols == spc
    xr = np.linspace(0.0, 1.0, ndivs[0])
    yr = np.linspace(0.0, 1.0, ndivs[1])
    zr = np.linspace(0.0, 1.0, ndivs[2])
    ic(xr.shape, yr.shape, zr.shape)
    X,Y,Z = np.meshgrid(xr,yr,zr,indexing='ij')
    grid_points= np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
    kde = stats.gaussian_kde(poss[to_consider].T)
    pdist0 = kde(grid_points)
    pdist = pdist0.reshape(ndivs)
    ic(grid_points.shape)
    ic(X.shape, Y.shape, Z.shape)
    ic(pdist0.shape, pdist.shape)
    return pdist

def get_prob_dist(ndivs, nsys, spc, sgm):
    tsgm2 = 2.0 * sgm * sgm
    tsgm2i = 1.0 / tsgm2
    pdist = np.zeros(ndivs, dtype=float)
    la, lb, lc = nsys.get_lattice_lengths()
    ra = 1.0/ndivs[0]
    rb = 1.0/ndivs[1]
    rc = 1.0/ndivs[2]
    dv = la*ra * lb*rb * lc*rc
    prefactor = dv / (np.pi*tsgm2)**1.5
    poss = nsys.get_scaled_positions()
    symbols = np.array(nsys.get_symbols(), dtype=str)
    to_consider = symbols == spc
    for ia in range(nsys.num_atoms()):
        if not to_consider[ia]:
            continue
        pi = poss[ia]
        for i in range(ndivs[0]):
            rai = ra*i -pi[0]
            rai = rai if abs(rai) < 0.5 else rai -1.0*np.sign(rai)
            da = la*rai
            if abs(da) > 2.5*sgm:
                continue
            for j in range(ndivs[1]):
                rbi = rb*j -pi[1]
                rbi = rbi if abs(rbi) < 0.5 else rbi -1.0*np.sign(rbi)
                db = lb*rbi
                if abs(db) > 2.5*sgm:
                    continue
                for k in range(ndivs[2]):
                    rci = rc*k -pi[2]
                    rci = rci if abs(rci) < 0.5 else rci -1.0*np.sign(rci)
                    dc = lc*rci
                    if abs(dc) > 2.5*sgm:
                        continue
                    dr2 = da*da +db*db +dc*dc
                    pdist[i,j,k] = pdist[i,j,k] \
                                   +prefactor*np.exp(-dr2 *tsgm2i)
    return pdist

def write_CHGCAR(nsys,pdist,fname='CHGCAR'):
    poscar = 'POSCAR_tmp'
    #nsys.write_POSCAR(fname=poscar)
    write(nsys,fname=poscar)
    with open(poscar,'r') as f:
        lines = f.readlines()
    with open(fname,'w') as f:
        for line in lines:
            f.write(line)
        f.write('\n')
        shape = np.shape(pdist)
        f.write(' {0:d} {1:d} {2:d}\n'.format(*shape))
        n = 0
        for k in range(shape[2]):
            for j in range(shape[1]):
                for i in range(shape[0]):
                    p = pdist[i,j,k]
                    f.write('  {0:12.4e}'.format(p))
                    if n % 5 == 4:
                        f.write('\n')
                    n += 1
    os.system('rm '+poscar)
    return None

def normalize_pdist(pdist, nsys, spc):
    print(' Normalize pdist...')
    s = np.sum(pdist)
    sid = nsys.species2sid(spc)
    assert sid > 0, "something is wrong with species2sid..."
    na = nsys.num_atoms(sid=sid)
    ic(na, s, sid)
    assert s > 1.0e-14, f"np.sum(pdist) is too small, {s}"
    pdist *= float(na)/s
    #print('   Max in pdist = ',np.max(pdist))
    #print('   Sum of pdist = ',np.sum(pdist))
    #print('   Num of atoms considered = ',na)
    return pdist

def main():
    from tqdm import tqdm

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    if args['-v']:
        ic.enable()
    files = args['FILES']
    spc = args['--species']
    if spc == 'None':
        spc = None
    width = float(args['-w'])
    sgm = float(args['--sigma'])
    if sgm < 0.0:
        sgm = abs(sgm)*width
    specorder = args['--specorder'].split(',')
    if 'None' in specorder:
        specorder = None
    if specorder == None and not (files[0].find('POSCAR') > -1 or
                                  files[0].find('pmd') > -1 ):
        raise ValueError('ERROR: specorder must be specified, unless files are POSCAR format.')

    #files.sort(key=get_key,reverse=True)

    print(f' Reading {files[0]}...')
    nsyss0 = read(fname=files[0], specorder=specorder)
    if type(nsyss0) is list:  # otherwise nsys0 is NAPSystem object
        nsys0 = nsyss0[0]
    else:
        nsys0 = nsyss0
    a, b, c = nsys0.get_lattice_angles()
    if abs(a-np.pi/2) > 180.0/np.pi or abs(b-np.pi/2) > 180.0/np.pi \
       or abs(c-np.pi/2) > 180.0/np.pi:
        raise ValueError('ERROR: Currently not available for non-orthogonal lattice,'
                         +' alpha,beta,gamma=', a, b, c)
    ndivs = get_num_division(nsys0, width)
    if np.dot(ndivs, ndivs) > 10000000:
        raise ValueError('ERROR: ndivs too large.', ndivs)
    print(' # of divisions = {0:d} {1:d} {2:d}'.format(*ndivs))

    pdist = np.zeros(ndivs, dtype=float)
    for i, f in enumerate(files):
        if i == 0:  # if i==0, file is already loaded
            nsyss = nsyss0
        else:
            print(f' Reading {f}...')
            nsyss = read(fname=f, specorder=specorder)
        if type(nsyss) is list:
            for nsys in tqdm(nsyss):
                pdist += get_prob_dist(ndivs, nsys, spc, sgm)
        else:
            pdist += get_prob_dist(ndivs, nsyss, spc, sgm)
    pdist = normalize_pdist(pdist, nsys0, spc)
    write_CHGCAR(nsys0, pdist)
    print(' Wrote CHGCAR')
    return None


if __name__ == "__main__":
    main()
