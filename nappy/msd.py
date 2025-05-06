#!/usr/bin/env python
"""
Compute the mean square displacement (MSD) of atoms
from the sequential structre structures that could come from several files or one file.

Usage:
  msd.py [options] FILES [FILES...]

Options:
  -h, --help    Show this help message and exit.
  -i, --iid IID
                Num of independent sampling. If iid > 1, more than 1 identical sampling 
                will be taken from the sequential data,
                and the duration time becomes shorter to 1/iid. [default: 1]
  --dt DT       Time interval (fs) between sequential structures. [default: -1.0]
  -o FILENAME   Output filename. [default: out.msd]
  --xyz         Decompose MSD to x,y,z-direction. [default: False]
  --com         Write MSD of center of motion (COM), which is not available with --xyz. [default: False]
  --err         Write error of each MSD. [default: False]
  --spcs SPC    Specify a species that is written out. [default: None]
"""
import sys
import numpy as np
import scipy.stats
from docopt import docopt
import copy
from datetime import datetime

import nappy
from nappy.napsys import NAPSystem
from nappy.util import gen_header

__author__ = "RYO KOBAYASHI"
__version__ = "250506"

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

    
def get_msd(structs, nmeasure, nshift, error=False):
    """
    Compute MSD from sequential structure FILES.
    
    Parameters
    ----------
    structs: list
         List of NAPSystem structures used for the MSD calculation.
    nmeasure: int
         Number of staggered lanes to compute MSD.
    nshift: int
         Number of frames to be skipped for each staggered lanes.
    error: logical
         Whether or not compute errors if nmeasure > 1.

    Returns
    -------
    msd : Numpy array of dimension (time-window,nspc,3).
    msdc : Numpy array of dimension (time-window,nspcs,3).
    err : Numpy array of dimension (time-window,nspc,3).
    errc : Numpy array of dimension (time-window,nspc,3).
    specorder: list
    """
    from tqdm import tqdm
    
    nsys0 = structs[0]
    specorder = copy.copy(nsys0.specorder)
    natm = len(nsys0)
    nspc = len(specorder)
    ids = [ i for i in range(nsys0.num_atoms()) ]
    naps = nsys0.natm_per_species()
    sids = nsys0.atoms.sid
    ispcs = np.array([ sid-1 for sid in sids ], dtype=int)

    symbols = nsys0.get_symbols()
    npbc= np.zeros((len(ids),3))
    hmat= np.zeros((3,3))
    unwrapped = np.zeros((len(structs),natm,3))
    pp= np.zeros((natm,3))

    for istr in tqdm(range(len(structs))):
        nsys= structs[istr]
        poss = nsys.get_scaled_positions()
        
        hmat = nsys.get_hmat()
        for ia,idi in enumerate(ids):
            # #...human-readable ID to computer-oriented ID
            # i= idi - 1
            pi= poss[idi]
            isp = ispcs[idi]
            if istr == 0:
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

            unwrapped[istr,ia,:] = np.dot(hmat,pi[:]+npbc[ia,:])

    nt = len(structs) -(nmeasure-1)*nshift
    p0= np.zeros((nmeasure,len(ids),3))
    pc0 = np.zeros((nmeasure,nspc,3))
    msd= np.zeros((nt, nmeasure, nspc, 3))
    msdc = np.zeros((nt,nmeasure, nspc, 3))

    print('   Processing data...')
    for istr in tqdm(range(len(structs))):

        for im in range(nmeasure):
            if istr == im*nshift:
                p0[im,:,:]= unwrapped[istr,:,:]
                for isp in range(nspc):
                    lspc = ispcs == isp
                    pc0[im,isp,:] = unwrapped[istr,lspc,:].mean(axis=0)
            if im*nshift <= istr < im*nshift +nt:
                #...normalized to absolute
                devs = unwrapped[istr,:,:] -p0[im,:,:]
                it = istr -im*nshift
                for isp in range(nspc):
                    lspc = ispcs == isp
                    msd[it,im,isp,:] += (devs[lspc,:]**2).sum(axis=0) /naps[isp]
                    pc = unwrapped[istr,lspc,:].mean(axis=0)
                    dc = pc -pc0[im,isp,:]
                    msdc[it,im,isp,:] += dc[:]**2 *naps[isp]

    print('')
    if error:
        return msd.mean(axis=1), scipy.stats.sem(msd, axis=1), \
            msdc.mean(axis=1), scipy.stats.sem(msdc, axis=1), specorder
    else:
        return msd.mean(axis=1), msdc.mean(axis=1), specorder

def main(files=[], dt=1.0, n_iid=1, xyz=False, lcom=False,
         error=False, spcs=None, outfname="out.msd"):

    # if nmeasure < 2:
    #     nmeasure = 1
    #     nshift = len(files)
    #     error = False
    #     print(' Since nmeasure < 2, nmeasure and nshift are set to ',nmeasure,nshift)
    # elif nshift < 0:
    #     nshift = int(len(files)/nmeasure)
    #     print(' Since nshift is not given, set nshift by len(files)/nmeasure = ',nshift)
    
    # #...compute sampling time-window from nmeasure and nshift
    # ntwindow= len(files) -(nmeasure-1)*nshift
    # if ntwindow <= 0:
    #     errmsg=' [Error] ntwindow <= 0 !!!\n'\
    #         + '  Chech the parameters nmeasure and nshift, and input files.'
    #     raise ValueError(errmsg)

    #files.sort(key=get_key,reverse=True)
    nsyss = []
    try:
        print(' Reading data files...')
        for f in files:
            nsyss_tmp = nappy.io.read(f)
            if type(nsyss_tmp) is list:
                nsyss.extend(nsyss_tmp)
            else:  # it must be NAPSystem object
                nsyss.append(nsyss_tmp)
    except:
        raise

    nshift = int(len(nsyss)/n_iid)
    ntwindow = len(nsyss) -(n_iid -1)*nshift

    print(' Computing MSD...')
    if error:
        msd,err,msdc,errc,specorder = get_msd(nsyss,n_iid,nshift,error=error)
    else:
        msd,msdc,specorder = get_msd(nsyss,n_iid,nshift,error=error)

    if spcs is None:
        spcs = specorder
    elif type(spcs) != list and type(spcs) == str:
        spcs = [spcs]
    else:
        raise TypeError('spcs is in invalid form: ', spcs)

    #...make output data files
    with open(outfname,'w') as f:
        f.write(gen_header(sys.argv))
        if dt > 0.0:
            f.write(f'# dt: {dt/1000:0.3f}    ! Time interval in ps\n')
        if xyz:
            icol = 1
            f.write(f'#   {icol:d}:data_ID,')
            icol += 1
            f.write(f'  {icol:d}:time(ps),')
            icol += 1
            for spc in specorder:
                if spc not in spcs: continue
                f.write(f'  {icol:d}-{icol+2:d}:msd_{spc:<2s} (x,y,z),             ')
                icol += 3
                if error:
                    f.write(f'  {icol:d}-{icol+2:d}:err_{spc:<2s} (x,y,z),                ')
                    icol += 3
            if lcom:
                for spc in specorder:
                    if spc not in  spcs: continue
                    f.write(f'   {icol:d}-{icol+2:d}:msdc_{spc:<2s} (x,y,z),             ')
                    icol += 3
                    if error:
                        f.write(f'   {icol:d}-{icol+2:d}:errc_{spc:<2s} (x,y,z),            ')
                        icol += 3
            f.write('\n')
            for idat in range(len(msd)):
                f.write(f' {idat:10d}  {idat*dt/1000:11.4e}')
                for isp,spc in enumerate(specorder):
                    if spc not in spcs: continue
                    f.write(' {0:11.3e} {1:11.3e} {2:11.3e}'.format(*msd[idat,isp,:]))
                    if error:
                        f.write(' {0:11.3e} {1:11.3e} {2:11.3e}'.format(*err[idat,isp,:]))
                if lcom:
                    for isp,spc in enumerate(specorder):
                        if spc not in spcs: continue
                        f.write(' {0:11.3e} {1:11.3e} {2:11.3e}'.format(*msdc[idat,isp,:]))
                        if error:
                            f.write(' {0:11.3e} {1:11.3e} {2:11.3e}'.format(*errc[idat,isp,:]))
                f.write('\n')
        else:  # not xyz
            icol = 1
            f.write(f'# {icol:d}:data_ID,')
            icol += 1
            f.write(f'  {icol:d}:time(ps),')
            for spc in specorder:
                if spc not in spcs: continue
                icol += 1
                f.write(f'   {icol:d}:msd_{spc:<2s}, ')
                if error:
                    icol += 1
                    f.write(f'   {icol:d}:err_{spc:<2s}, ')
            if lcom:
                for spc in specorder:
                    if spc not in spcs: continue
                    icol += 1
                    f.write(f'   {icol:d}:msdc_{spc:<2s},')
                    if error:
                        icol += 1
                        f.write(f'   {icol:d}:errc_{spc:<2s},')
            f.write('\n')
            for idat in range(len(msd)):
                f.write(f' {idat:10d}  {idat*dt/1000:11.4e}')
                for isp,spc in enumerate(specorder):
                    if spc not in spcs: continue
                    f.write(f' {msd[idat,isp].sum():12.3e}')
                    if error:
                        f.write(f' {err[idat,isp].sum():12.3e}')
                if lcom:
                    for isp,spc in enumerate(specorder):
                        if spc not in spcs: continue
                        f.write(f' {msdc[idat,isp].sum():12.3e}')
                        if error:
                            f.write(f' {errc[idat,isp].sum():12.3e}')
                f.write('\n')
    print(' Wrote a file: {0:s}'.format(outfname))
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    files = args['FILES']
    if len(files) == 0:
        raise ValueError('No file is given.')
    n_iid = int(args['--iid'])
    dt = float(args['--dt'])
    if dt < 0.0:
        raise ValueError('Current version of msd.py requires --dt option to be set by the user.\n'
                         +'See the help with -h option.')
    outfname= args['-o']
    xyz = args['--xyz']
    lcom = args['--com']
    error = args['--err']
    spcs = args['--spcs']
    if spcs == 'None':
        spcs = None

    main(files=files,
         dt=dt,
         n_iid=n_iid,
         xyz=xyz,
         lcom=lcom,
         error= error,
         spcs= spcs,
         outfname=outfname)
