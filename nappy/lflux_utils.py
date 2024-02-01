#!/usr/bin/env python
"""
IO utilities that can be used in the project, lflux.

Usage:
  lflux_utils.py [options] INFILE

Options:
  -h, --help  Show this message and exit.
  --nproc NPROC  Number of processes. [default: 1]
  --noutlflux NOUT  Number of output data in out.lflux. [default: 200]
  --sysfile SYSFILE  File that contains the base system information. [default: dump_0]
"""

import os,sys
from docopt import docopt
import numpy as np
import copy
from tqdm import tqdm

__author__ = "RYO KOBAYASHI"
__version__ = "230115"

_ang2bohr = 1.0/0.529177
_bohr2ang = 1.0 /_ang2bohr


def read_gcube(fname):
    with open(fname,'r') as f:
        comment1 = f.readline()
        comment2 = f.readline()
        data = f.readline().split() # 3rd
        natm = int(data[0])
        orig = np.array([ float(d) for d in data[1:4] ])*_bohr2ang
        data = f.readline().split() # 4th
        n1 = int(data[0])
        da1 = np.array([ float(d) for d in data[1:4]]) *_bohr2ang
        data = f.readline().split() # 5th
        n2 = int(data[0])
        da2 = np.array([ float(d) for d in data[1:4]]) *_bohr2ang
        data = f.readline().split() # 6th
        n3 = int(data[0])
        da3 = np.array([ float(d) for d in data[1:4]]) *_bohr2ang
        nd = n1*n2*n3
        # Make mesh
        xs = [ da1[0]*(i1+0.5) for i1 in range(n1) ]
        ys = [ da2[1]*(i2+0.5) for i2 in range(n2) ]
        zs = [ da3[2]*(i3+0.5) for i3 in range(n3) ]
        X,Y,Z = np.meshgrid(xs,ys,zs)
        #X,Y,Z = np.mgrid[0.0:da[0]*ngx:da[0], 0.0:db[1]*ngy,db[1], 0.0:dc[2]*ngz,dc[2] ]
        # Skip atoms
        for ia in range(natm):
            f.readline()
        origdata = np.zeros(nd)
        ig = 0
        for line in f.readlines():
            data = [ float(d) for d in line.split() ]
            for i in range(len(data)):
                origdata[ig] = data[i]
                ig += 1
        voldata = np.zeros((n1,n2,n3))
        ig = 0
        for i1 in range(n1):
            for i2 in range(n2):
                for i3 in range(n3):
                    voldata[i1,i2,i3] = origdata[ig]
                    ig += 1
    return natm,n1,n2,n3,da1,da2,da3,orig,voldata

def write_gcube(n1,n2,n3,a1,a2,a3,fts,fname='out.lflux.cube',
                org=[0.0,0.0,0.0], pos=[]):
    """
    Write local flux data in Gaussian cube format.
    The shape of the given volumetric data, fts, should be 1-dimensional.
    """

    if len(fts.shape) == 1:  # flat 1D array
        if fts.shape[0] != n1*n2*n3:
            raise ValueError('Mismatch array length.')
    elif len(fts.shape) == 3:  # 3D array
        if fts.shape != (n1,n2,n3):
            raise ValueError('Mismatch array shape.')
        fts.reshape((n1*n2*n3,),order='C')
    if len(a1) != 3 or len(a2) != 3 or len(a3) != 3:
        raise ValueError('Lattice vectors a1,a2,a3 must be set appropriately.')
    if type(a1) != np.ndarray or type(a2) != np.ndarray or type(a3) != np.ndarray:
        a1 = np.array(a1)
        a2 = np.array(a2)
        a3 = np.array(a3)
    
    with open(fname,'w') as f:
        f.write('Gaussian cube for local flux data,\n')
        f.write('written from a jupyter notebook.\n')
        # f.write('  1   {0:12.5f}  {1:12.5f}  {2:12.5f}\n'.format(la/ngx/2*ang2bohr,   # Including one dummy atom
        #                                                          lb/ngy/2*ang2bohr,   # slightly shifted origin
        #                                                          lc/ngz/2*ang2bohr))
        f.write('  {0:d}   {1:12.5f}  {2:12.5f}  {3:12.5f}\n'.format(max(1,len(pos)),
                                                                     org[0]*_ang2bohr,
                                                                     org[1]*_ang2bohr,
                                                                     org[2]*_ang2bohr))  
        # In the official site of Gaussian cube format, the sign of the number of grid in each direction
        # has the meaning of unit; Bohr (positive) or Ang (negative), but VMD only works with positive values.
        da1 = a1/n1*_ang2bohr
        da2 = a2/n2*_ang2bohr
        da3 = a3/n3*_ang2bohr
        f.write(' {0:d}  {1:12.5f}  {2:12.5f}  {3:12.5f}\n'.format(n1,da1[0],da1[1],da1[2]))
        f.write(' {0:d}  {1:12.5f}  {2:12.5f}  {3:12.5f}\n'.format(n2,da2[0],da2[1],da2[2]))
        f.write(' {0:d}  {1:12.5f}  {2:12.5f}  {3:12.5f}\n'.format(n3,da3[0],da3[1],da3[2]))
        # Add at least one atom (dummy)
        if len(pos) < 1:
            f.write(' 1  1.000    0.000    0.000    0.000\n')
        else:
            for pi in pos:
                f.write('  1  1.000   {0:10.3f}  {1:10.3f}  {2:10.3f}\n'.format(*pi))
        # Volumetric data hereafter
        if len(fts) < 0:
            f.write('0.0\n')
        else:
            if len(fts.shape) == 1:
                vtxt = ''
                inc = 0
                for vd in fts:
                    vtxt += f' {vd:12.4e}'
                    inc += 1
                    if inc % 6 == 0:
                        vtxt += '\n'
                        inc = 0
                vtxt += '\n'
                f.write(vtxt)
    return

def read_lflux(fname='out.lflux',noutlflux=1000):
    with open(fname,'r') as f:
        inc = 0
        time = np.zeros(noutlflux+1)
        for il,line in tqdm(enumerate(f.readlines())):
            if line[0] == '#':
                if 'ngx,ngy,ngz,ng' in line:
                    print(line)
                    data = line.split()
                    ngx,ngy,ngz,ng = [ int(x) for x in data[7:11] ]
                    print('ngx,ngy,ngz,ng = ',ngx,ngy,ngz,ng)
                    dflux = np.zeros((noutlflux+1,ng))
                continue
            #print(line)
            data = line.split()
            if len(data) < 2:
                continue
            #print(data[0:5])
            istp = int(data[0])
            time[inc] = float(data[1])
            dflux[inc,0:ng] = [ float(x) for x in data[2:ng+2] ]
            inc += 1
            #print('#',end='')
    return ngx,ngy,ngz,ng,dflux,time

def func1(X,a,b):
    Y = a*X +b
    return Y

def fit_dflux_time(dflux,time):
    from scipy.optimize import curve_fit
    nout, ng = dflux.shape
    fts = np.zeros(ng)
    errs = np.zeros(ng)
    for ig in range(ng):
        flux = dflux[:,ig]
        popt,pcov = curve_fit(func1,time[int(nout/2):nout],flux[int(nout/2):nout])
        fts[ig] = popt[0]
        errs[ig] = pcov[0,0]
        if ig % int(ng/10) == 0:
            print('  {0:5d}/{1:5d}'.format(ig,ng))
    return fts,errs

def fit_flux_time_grid(ig,nout,flux,time):
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func1,time[int(nout/2):nout],flux[int(nout/2):nout])
    return popt[0], pcov[0,0]

def conv_lflux_to_gcube(infname='out.lflux',noutlflux=1000,outfname='out.lflux.cube',
                        a1=[], a2=[], a3=[], nproc=1):
    """
    Convert out.lflux file to Gaussian cube format.
    """
    if len(a1) != 3 or len(a2) != 3 or len(a3) != 3:
        raise ValueError('Lattice vectors a1,a2,a3 must be given appropriately.')
    
    from multiprocess import Pool
    if nproc > 0:
        pool = Pool(processes=nproc)
    else:
        pool = Pool()

    print(' Reading out.lflux...')
    nx,ny,nz,ncell,dflux,time = read_lflux(fname=infname, noutlflux=noutlflux)
    nout, ng = dflux.shape
    fts = np.zeros(ng)
    errs = np.zeros(ng)

    print(' Processing grids...')
    # inc = 0
    for igg in tqdm(range(0,ng,nproc)):
        prcs = []
        for ip in range(nproc):
            ig = igg +ip
            if ig >= ng: continue
            flux = dflux[:,ig]
            prcs.append(pool.apply_async(fit_flux_time_grid,
                                         (ig,nout,flux,time) ))
        results = [ res.get() for res in prcs ]
        for ip in range(nproc):
            ig = igg +ip
            if ig >= ng: continue
            ft, err = results[ip]
            fts[ig] = ft
            errs[ig] = err
        # if igg > int(ng/10*inc):
        #     print('  {0:5d}/{1:5d}'.format(igg,ng))
        #     inc += 1
    #fts, errs = fit_dflux_time(dflux,time)
    print(' Writing cube file...')
    write_gcube(nx,ny,nz,a1,a2,a3,fts,fname=outfname)
    return None


if __name__ == "__main__":
    import nappy
    
    args = docopt(__doc__)
    nproc = int(args['--nproc'])
    nout = int(args['--noutlflux'])
    sysfile = args['--sysfile']
    infile = args['INFILE']

    nsys = nappy.io.read(sysfile)
    a1,a2,a3 = nsys.get_lattice_vectors()
    outfile = infile +'.cube'
    conv_lflux_to_gcube(infname=infile,
                        noutlflux=nout,
                        outfname=outfile,
                        a1=a1, a2=a2, a3=a3,
                        nproc=nproc)
