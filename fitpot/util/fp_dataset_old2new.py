#!/usr/bin/env python
"""
Convert old-style fitpot format to new-style one.

Usage:
  {0:s} [options] DATASETDIR

Options:
  -h, --help  Show this message and exit.
"""
import os,sys
from docopt import docopt
import numpy as np
import glob
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = ""

def read_ref_data(dname):
    with open(dname+'/erg.ref','r') as f:
        erg = float(f.readline().split()[0])
    with open(dname+'/strs.ref','r') as f:
        line = f.readline()
        dat = line.split()
        strs = np.array([ float(d) for d in dat ])
    with open(dname+'/frc.ref','r') as f:
        lines = f.readlines()
    for il,line in enumerate(lines):
        if il == 0:
            natm = int(line.split()[0])
            frcs = np.zeros((natm,3),dtype=float)
        else:
            frcs[il-1,:] = [ float(x) for x in line.split() ]
    return erg, frcs, strs

def write_smpl(smplname,nsys,erg,frcs,strs):
    with open(smplname,'w') as f:
        f.write('#\n')
        f.write('# specorder:')
        for s in nsys.specorder:
            f.write(f'  {s}')
        f.write('\n')
        f.write(f'# energy:  {erg:0.5f}\n')
        f.write(f'# stress:')
        for s in strs:
            f.write(f'  {s:0.5f}')
        f.write('\n')
        f.write('# auxiliary_data:  fx  fy  fz\n')
        f.write('#\n')
        f.write('    1.000\n')
        h = nsys.get_hmat()
        f.write(f'  {h[0,0]:10.5f}  {h[1,0]:10.5f}  {h[2,0]:10.5f}  0.00  0.00  0.00\n')
        f.write(f'  {h[0,1]:10.5f}  {h[1,1]:10.5f}  {h[2,1]:10.5f}  0.00  0.00  0.00\n')
        f.write(f'  {h[0,2]:10.5f}  {h[1,2]:10.5f}  {h[2,2]:10.5f}  0.00  0.00  0.00\n')
        f.write(f'  {len(nsys):d}\n')
        sposs = nsys.get_scaled_positions()
        svels = nsys.get_scaled_velocities()
        tags = nsys.get_tags()
        for i in range(len(nsys)):
            ti = tags[i]
            spi = sposs[i]
            fi = frcs[i]
            f.write( f'  {ti:16.14f}  {spi[0]:22.14e} {spi[1]:22.14e} {spi[2]:22.14e}'
                    +f'  0.00 0.00 0.00'
                    +f'  {fi[0]:10.4f} {fi[1]:10.4f} {fi[2]:10.4f}\n')
        return None

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    datasetdir = args['DATASETDIR']
    dirs = glob.glob(datasetdir+'/smpl_*')
    new_datasetdir = datasetdir.rstrip('/')+'_new'
    os.mkdir(new_datasetdir)
    for d in dirs:
        print('.',end='',flush=True)
        nsys = nappy.io.read(d+'/pos',format='pmd')
        erg, frcs, strs = read_ref_data(d)
        fname = d.replace(datasetdir.rstrip('/'), new_datasetdir).rstrip('/')
        write_smpl(fname, nsys, erg, frcs, strs)

if __name__ == "__main__":

    main()
