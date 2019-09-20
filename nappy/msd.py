#!/usr/bin/env python
"""
Compute the mean square displacement (MSD) of atoms (see ids below)
from the sequential files in the arguments.
If several atom IDs are specified, an average MSD among those atoms are written.
Staggered measuring of MSD can be used for the statistical purpose.

Usage:
  msd.py [options] FILES [FILES...]

Options:
  -h, --help    Show this help message and exit.
  --id=ID       IDs of atoms whose paths are to be traced. 
                0 for all atoms and positive for atom-IDs.
                Several IDs can be specified separated by comma.
                [default: 0,]
  -m, --measure MEASURE
                Num of measuring lane. In case of 1, it is identical to non-staggered measuring. [default: 1]
  -s, --shift SHIFT
                Shift of each staggered lane. [default: 20]
  --spcs=SPCS   Species name whose MSD is to be computed. [default: None]
"""
from __future__ import print_function

import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem
from nappy.common import get_key

def anint(x):
    if x >= 0.5:
        return 1.0
    elif x < -0.5:
        return -1.0
    else:
        return 0.0

def get_ids(nsys,ids):
    atom_ids = []
    if 0 in ids:  # all the atoms
        atom_ids = [ i for i in range(len(nsys)) ]
        return atom_ids
    else:
        atom_ids = [ids]
    return atom_ids

    
def get_msd(files,ids0,nmeasure,nshift,spcs='None'):
    """
    Compute MSD specified species-ID from sequential structure FILES.
    
    Input:

      - files: list
            List of files used for the MSD calculation.
      - ids0: list
            List of atom-IDs (starting from 1) whose MSDs are to be computed.
      - nmeasure: int
            Number of staggered lanes to compute MSD for better statistics.
      - nshift: int
            Number of files to be skipped for each staggered lane.
      - spcs: string
            Species name. If it is None, no species is specified. [defalut: None]

    Output:
      - msd: Numpy array of dimension, (len(files),nmeasure,3).
    """
    if spcs != 'None':
        nsys = NAPSystem(fname=files[0])
        symbols = nsys.get_symbols()
        ids = [ i for i,s in enumerate(symbols) if s == spcs ]
    else:
        if 0 in ids0:
            nsys = NAPSystem(fname=files[0])
            ids = [ i for i in range(nsys.natm) ]
        else:
            ids = [ i-1 for i in ids0 ]

    p0= np.zeros((nmeasure,len(ids),3))
    pp= np.zeros((len(ids),3))
    msd= np.zeros((len(files),nmeasure,3))
    npbc= np.zeros((len(ids),3))
    hmat= np.zeros((3,3))
    for ifile in range(len(files)):
        fname= files[ifile]
        nsys= NAPSystem(fname=fname)
        hmat = nsys.get_hmat()
        for ia,idi in enumerate(ids):
            # #...human-readable ID to computer-oriented ID
            # i= idi - 1
            pi= nsys.poss[idi]
            if ifile == 0:
                pp[ia,:]= pi[:]
            else:
                #...correct periodic motion
                dev= pi -pp[ia]
                if dev[0] > 0.5:
                    npbc[ia,0] += -1.0
                elif dev[0] < -0.5:
                    npbc[ia,0] += 1.0
                if dev[1] > 0.5:
                    npbc[ia,1] += -1.0
                elif dev[1] < -0.5:
                    npbc[ia,1] += 1.0
                if dev[2] > 0.5:
                    npbc[ia,2] += -1.0
                elif dev[2] < -0.5:
                    npbc[ia,2] += 1.0
                # print npbc
                #...store current position
                pp[ia,:]= pi[:]

            for nm in range(nmeasure):
                if ifile == nm*nshift:
                    p0[nm,ia,0]= pi[0] +npbc[ia,0]
                    p0[nm,ia,1]= pi[1] +npbc[ia,1]
                    p0[nm,ia,2]= pi[2] +npbc[ia,2]
                if nm*nshift < ifile:
                    #...normalized to absolute
                    dev[0]= pi[0] +npbc[ia,0] -p0[nm,ia,0] 
                    dev[1]= pi[1] +npbc[ia,1] -p0[nm,ia,1] 
                    dev[2]= pi[2] +npbc[ia,2] -p0[nm,ia,2] 
                    dev= np.dot(hmat,dev)
                    msd[ifile-nm*nshift,nm,0] += dev[0]**2 
                    msd[ifile-nm*nshift,nm,1] += dev[1]**2 
                    msd[ifile-nm*nshift,nm,2] += dev[2]**2
                        

    for ifile in range(len(files)):
        for nm in range(nmeasure):
            msd[ifile-nm*nshift,nm,0] /= len(ids)
            msd[ifile-nm*nshift,nm,1] /= len(ids)
            msd[ifile-nm*nshift,nm,2] /= len(ids)
                    
    
    return msd

# def get_key(v):
#     import re
#     prefix, index = re.match(r'([a-z]+)_(\d+)', v).groups()
#     return prefix, -int(index)

if __name__ == "__main__":

    args = docopt(__doc__)

    files = args['FILES']
    ids = args['--id']
    spcs = args['--spcs']
    if spcs == 'None':
        ids = [ int(i) for i in ids.split(',') ]
    nmeasure = int(args['--measure'])
    nshift = int(args['--shift'])
    
    #...compute sampling time-window from nmeasure and nshift
    ntwindow= len(files) -(nmeasure-1)*nshift
    if ntwindow <= 0:
        err=' [Error] ntwindow <= 0 !!!\n'\
              + '  Chech the parameters nmeasure and nshift, and input files.'
        raise ValueError(err)

    files.sort(key=get_key,reverse=True)
    # for i in range(len(files)):
    #     print i,files[i]
    
    msd = get_msd(files,ids,nmeasure,nshift,spcs)

    #...make output data files
    outfname='out.msd'
    with open(outfname,'w') as f:
        f.write('#   data_ID,      msd_total,      msd_x,          msd_y,          msd_z\n')
        for ifile in range(len(files)-(nmeasure-1)*nshift):
            if ifile == 0:
                f.write(' {:10d}'.format(ifile)
                        +' {:15.7f} {:15.7f}'.format(0.0,0.0)
                        +' {:15.7f} {:15.7f}\n'.format(0.0,0.0))
            else:
                dev2= np.zeros((3,))
                for nm in range(nmeasure):
                    dev2 += msd[ifile,nm]
                dev2 /= nmeasure
                f.write(' {:10d}'.format(ifile)
                        +' {:15.7f}'.format((dev2[0]+dev2[1]+dev2[2]))
                        +' {:15.7f}'.format(dev2[0])
                        +' {:15.7f}'.format(dev2[1])
                        +' {:15.7f}'.format(dev2[2])
                        +' \n')
    print('Wrote a file: {0:s}'.format(outfname))
