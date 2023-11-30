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
                Several IDs can be specified separated by comma.
                If negative value is included in the list, all atoms are taken into account.
                [default: -1,]
  -m, --measure MEASURE
                Num of measuring lane. In case of 1, it is identical to non-staggered measuring. [default: 1]
  -s, --shift SHIFT
                Shift of each staggered lane. [default: -1]
  --dt DT       Time interval (fs) between sequential files. [default: -1.0]
  -o FILENAME   Output filename. [default: out.msd]
  --xyz         Decompose MSD to x,y,z-direction. [default: False]
  --com         Write MSD of center of motion (COM), which is not available with --xyz. [default: False]
  --specorder SPECORDER
                Species order of the given system, separated by comma. [default: None]
"""
import sys
import numpy as np
from docopt import docopt
import copy
from datetime import datetime

from nappy.napsys import NAPSystem
from nappy.common import get_key
from nappy.io import read
from nappy.util import gen_header

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
        atom_ids = [ i for i in range(nsys.num_atoms()) ]
        return atom_ids
    else:
        atom_ids = [ids]
    return atom_ids

    
def get_msd(files, ids0, nmeasure, nshift, specorder=None):
    """
    Compute MSD of specified species-ID from sequential structure FILES.
    
    Parameters
    ----------
    files: list
         List of files used for the MSD calculation.
    ids0: list
         List of atom-IDs (starting from 1) whose MSDs are to be computed.
    nmeasure: int
         Number of staggered lanes to compute MSD for better statistics.
    nshift: int
         Number of files to be skipped for each staggered lane.
    specorder: list
         Order of species.

    Returns
    -------
    msd : Numpy array of dimension, (len(files),nmeasure,nspc,3).
    specorder: list
    """
    from tqdm import tqdm

    nsys = read(fname=files[0],specorder=specorder)
    if specorder is None:
        specorder = copy.copy(nsys.specorder)
    
    nspc = len(specorder)
    if ids0 is not None:
        ids = [ i-1 for i in ids0 ]
        sids = nsys.atoms.sid
        naps = [ 0 for i in range(len(specorder)) ]
        for i in ids:
            sid = sids[i]
            naps[sid-1] += 1
    else:
        ids = [ i for i in range(nsys.num_atoms()) ]
        naps = nsys.natm_per_species()

    symbols = nsys.get_symbols()
    p0= np.zeros((nmeasure,len(ids),3))
    pp= np.zeros((len(ids),3))
    com = np.zeros((len(files),nspc,3))
    com0 = np.zeros((nmeasure,nspc,3))
    # msd= np.zeros((len(files),nmeasure,nspc,3))
    msd= np.zeros((len(files)-(nmeasure-1)*nshift+1, nmeasure, nspc, 3))
    msdcom= np.zeros((len(files)-(nmeasure-1)*nshift+1, nmeasure, nspc, 3))
    npbc= np.zeros((len(ids),3))
    hmat= np.zeros((3,3))
    
    for ifile in tqdm(range(len(files))):
        fname= files[ifile]
        # sys.stdout.write('\r{0:5d}/{1:d}: {2:s}'.format(ifile+1,len(files),fname),)
        # sys.stdout.flush()
        if ifile != 0:
            nsys = read(fname=fname,specorder=specorder)
        poss = nsys.get_scaled_positions()
        sids = nsys.atoms.sid
        
        hmat = nsys.get_hmat()
        for ia,idi in enumerate(ids):
            # #...human-readable ID to computer-oriented ID
            # i= idi - 1
            pi= poss[idi]
            sid = sids[idi] -1
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

            com[ifile,sid,:] += pi[:] +npbc[ia,:]
            for nm in range(nmeasure):
                if ifile == nm*nshift:
                    com0[nm,sid,:] += pi[:] +npbc[ia,:]

            for nm in range(nmeasure):
                if ifile == nm*nshift:
                    p0[nm,ia,0]= pi[0] +npbc[ia,0]
                    p0[nm,ia,1]= pi[1] +npbc[ia,1]
                    p0[nm,ia,2]= pi[2] +npbc[ia,2]
                if nm*nshift < ifile <= (nm+1)*nshift:
                    #...normalized to absolute
                    dev[0]= pi[0] +npbc[ia,0] -p0[nm,ia,0]
                    dev[1]= pi[1] +npbc[ia,1] -p0[nm,ia,1]
                    dev[2]= pi[2] +npbc[ia,2] -p0[nm,ia,2]
                    dev= np.dot(hmat,dev)
                    msd[ifile-nm*nshift,nm,sid,0] += dev[0]**2
                    msd[ifile-nm*nshift,nm,sid,1] += dev[1]**2
                    msd[ifile-nm*nshift,nm,sid,2] += dev[2]**2
                        
    for ifile in range(len(files)):
        for sid in range(nspc):
            com[ifile,sid,:] = np.dot(hmat,com[ifile,sid,:])
    for nm in range(nmeasure):
        for sid in range(nspc):
            com0[nm,sid,:] = np.dot(hmat,com0[nm,sid,:])

    for ifile in range(len(files)):
        for sid in range(nspc):
            for nm in range(nmeasure):
                if nm*nshift < ifile <= (nm+1)*nshift:
                    dev = com[ifile,sid,:] -com0[nm,sid,:]
                    msdcom[ifile-nm*nshift,nm,sid,:] = dev[:]**2 /naps[sid]

    for ifile in range(len(files)):
        for nm in range(nmeasure):
            if nm*nshift < ifile <= (nm+1)*nshift:
                #...NOTE: The code below could cause true_divide error,
                #...  since when atoms are specified via --ids,
                #...  any of naps elements could be zero...
                msd[ifile-nm*nshift,nm,:,0] /= naps[:]
                msd[ifile-nm*nshift,nm,:,1] /= naps[:]
                msd[ifile-nm*nshift,nm,:,2] /= naps[:]

    print('')
    return msd,msdcom,specorder

def main(files=[], dt=1.0, nmeasure=1, nshift=-1, ids=None,
         specorder=[], xyz=False, lcom=False):

    if nmeasure < 2:
        nmeasure = 1
        nshift = len(files)
        print(' Since nmeasure < 2, nmeasure and nshift are set to ',nmeasure,nshift)
    elif nshift < 0:
        nshift = int(len(files)/nmeasure)
        print(' Since nshift is not given, set nshift by len(files)/nmeasure = ',nshift)
    
    #...compute sampling time-window from nmeasure and nshift
    ntwindow= len(files) -(nmeasure-1)*nshift
    if ntwindow <= 0:
        err=' [Error] ntwindow <= 0 !!!\n'\
            + '  Chech the parameters nmeasure and nshift, and input files.'
        raise ValueError(err)

    files.sort(key=get_key,reverse=True)
    
    msd,msdcom,specorder = get_msd(files,ids,nmeasure,nshift,specorder)

    #...make output data files
    with open(outfname,'w') as f:
        f.write(gen_header(sys.argv))
        if dt > 0.0:
            f.write(f'# dt: {dt:0.1f}    ! Time interval in fs\n')
        if xyz:
            f.write('#   data_ID,')
            for spc in specorder:
                f.write(f'  msd_{spc:<2s} (x,y,z),                   ')
            f.write('\n')
            for idat in range(len(files)-(nmeasure-1)*nshift):
                if idat == 0:
                    f.write(f' {idat:10d}')
                    for isp in range(len(specorder)):
                        f.write(' {0:11.3e} {0:11.3e} {0:11.3e}'.format(0.0,))
                    f.write('\n')
                else:
                    f.write(f' {idat:10d}')
                    for isp in range(len(specorder)):
                        dev2= np.zeros(3)
                        for nm in range(nmeasure):
                            dev2[:] += msd[idat,nm,isp,:]
                        dev2 /= nmeasure
                        f.write(' {:11.3e}'.format(dev2[0])
                                +' {:11.3e}'.format(dev2[1])
                                +' {:11.3e}'.format(dev2[2]) )
                    f.write('\n')
        else:  # not xyz
            f.write('#   data_ID,')
            for spc in specorder:
                f.write(f'   msd_{spc:<2s},   ')
            if lcom:
                for spc in specorder:
                    f.write(f'   msdcom_{spc:<2s},')
            f.write('\n')
            for idat in range(len(files)-(nmeasure-1)*nshift):
                if idat == 0:
                    f.write(f' {idat:10d}')
                    for isp in range(len(specorder)):
                        f.write(' {:12.3e}'.format(0.))
                    if lcom:
                        for isp in range(len(specorder)):
                            f.write(' {:12.3e}'.format(0.))
                    f.write('\n')
                else:
                    f.write(f' {idat:10d}')
                    for isp in range(len(specorder)):
                        dev2 = np.zeros(3)
                        for nm in range(nmeasure):
                            dev2[:] += msd[idat,nm,isp,:]
                        dev2 /= nmeasure
                        f.write(' {:12.3e}'.format(dev2[0]+dev2[1]+dev2[2]))
                    if lcom:
                        for isp in range(len(specorder)):
                            dev2 = np.zeros(3)
                            for nm in range(nmeasure):
                                dev2[:] += msdcom[idat,nm,isp,:]
                            dev2 /= nmeasure
                            f.write(' {:12.3e}'.format(dev2[0]+dev2[1]+dev2[2]))
                    f.write('\n')
    print(' Wrote a file: {0:s}'.format(outfname))
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    files = args['FILES']
    if len(files) == 0:
        raise ValueError('No file is given.')
    ids = [ int(i) for i in args['--id'].split(',') if i != '' ]
    for i in ids:
        if i < 0:
            ids = None
            break
    nmeasure = int(args['--measure'])
    nshift = int(args['--shift'])
    dt = float(args['--dt'])
    if dt < 0.0:
        raise ValueError('Current version of msd.py requires --dt option to be set by the user.\n'
                         +'See the help with -h option.')
    outfname= args['-o']
    specorder = args['--specorder']
    if specorder == 'None' or specorder is None:
        specorder = None
    else:
        specorder = [ s for s in specorder.split(',')]
    xyz = args['--xyz']
    lcom = args['--com']

    main(files=files,
         dt=dt,
         nmeasure=nmeasure,
         nshift=nshift,
         ids=ids,
         specorder=specorder,
         xyz=xyz,
         lcom=lcom)
