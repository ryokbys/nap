#!/usr/bin/env python
"""
Calculate the velocity auto correlation function from pmd/dump files.

Staggered measuring of auto correlation for the statistical purpose is available with -m and -s options.

Velocity data must be stored at 1st-3rd auxiliary data of each atomic data in pmd/dump files.

Usage:
  vel_auto_corr.py [options] INFILE [INFILE...]

Options:
  -h, --help   Show this help message and exit.
  -m, --measure MEASURE
               Num of measuring lane. In case of 1, it is identical to non-staggered measuring. [default: 1]
  -s, --shift SHIFT
               Shift of each staggered lane. [default: 20]
  -t, --time-interval TIME_INTERVAL
               Time interval between successive files in fs. [default: 1.0]
  --sigma SIGMA
               Sigma for Gaussian smoothing by integer. [default: 0]
"""
import os,sys,glob,time,copy
from datetime import datetime
import numpy as np
from docopt import docopt

from nappy.io import read
from nappy.common import get_key

from pwtools.signal import pad_zeros, welch
from scipy.fftpack import fft
from scipy.ndimage import gaussian_filter

__author__ = "Ryo KOBAYASHI"
__version__ = "230428"

def main(args):

    #...num of measuring lane, in case of 1, it is identical to non-staggered measuring
    nmeasure = int(args['--measure'])
    
    #...shift of each staggered lane
    nshift = int(args['--shift'])

    #...Ggaussian sigma
    sgm = int(args['--sigma'])
    
    #...time interval
    dt = float(args['--time-interval'])
    
    print(' nmeasure=',nmeasure)
    print(' nshift=',nshift)
    print(' dt=',dt,'fs')
    print(' len(args)=',len(args))
    
    #...parse arguments
    infiles = args['INFILE']
    infiles.sort(key=get_key,reverse=True)
    
    #...compute sampling time-window from nmeasure and nshift
    ntw= len(infiles) -(nmeasure-1)*nshift
    tmax = ntw *dt
    tdamp = tmax/np.sqrt(3.0)
    print(' ntw =',ntw)
    print(' tmax     =',tmax,'fs')
    if ntw <= 0:
        print(' [Error] ntw <= 0 !!!')
        print('  Chech the parameters nmeasure and nshift, and input files.')
        sys.exit()

    #...set global values
    infile= infiles[0]

    nsys0 = read(fname=infile)
    natm= len(nsys0)
    nspcs = len(nsys0.specorder)
    n_per_spcs = nsys0.natm_per_species()
    sids0 = nsys0.atoms['sid']
    print(' num of all atoms = ',natm)
    print(' num of atoms per species = ',n_per_spcs)
    
    print(' accumurating data',end='')
    actmp= np.zeros((nmeasure,ntw,3))
    acorr= np.zeros((ntw,3))
    v0= np.zeros((nmeasure,natm,3))
    v2= np.zeros((ntw,nmeasure,natm))
    for ifile in range(len(infiles)):
        print('.',end='')
        sys.stdout.flush()
        infile= infiles[ifile]

        nsys = read(fname=infile)
        hmat= nsys.get_hmat()
        vels = nsys.get_velocities(real=True)
    
        for im in range(nmeasure):
            if ifile == im*nshift:
                for ia in range(natm):
                    v0[im,ia,:] = vels[ia,:]
            if ifile >= im*nshift and ifile-im*nshift < ntw:
                for ia in range(natm):
                    v2[ifile-im*nshift,im,ia] = np.dot(vels[ia,:],v0[im,ia,:])
    print('')
    
    #.....output auto correlation function
    acfname='dat.autocorr'
    vdenom = np.zeros(nspcs)
    for im in range(nmeasure):
        for ia in range(natm):
            ispc = sids0[ia] -1
            vdenom[ispc] += np.dot(v0[im,ia,:],v0[im,ia,:])

    acfile= open(acfname,'w')
    acfile.write('# Velocity auto correlation from vel_auto_corr.py ' +
                 'at {0:s}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    acfile.write('#      t [fs],   Cvv(t) for each species,    sum of species-Cvv(t)\n')
    ac= np.zeros((ntw,nspcs))
    acm= np.zeros(nspcs)
    for it in range(ntw):
        for im in range(nmeasure):
            acm[:] = 0.0
            for ia in range(natm):
                ispc = sids0[ia] -1
                acm[ispc] += v2[it,im,ia]
            for ispc in range(nspcs):
                ac[it,ispc] += acm[ispc]
        ac[it,:] = ac[it,:] /vdenom[:]
        acfile.write(' {0:11.3e}'.format(it*dt) )
        sumac = 0.0
        for ispc in range(nspcs):
            acfile.write(' {0:11.3e}'.format(ac[it,ispc]))
            sumac += ac[it,ispc]
        acfile.write(f' {sumac:11.3e}\n')
    acfile.close()

    #...Power spectrum
    mwmax = int(ntw/2)
    dw = 2.0*np.pi /tmax
    ps0 = np.zeros((mwmax,nspcs))
    for mw in range(mwmax):
        w = dw *mw
        for it in range(ntw):
            t = it*dt
            decay = np.exp(-(t/tdamp)**2)
            coswt = np.cos(w*t)
            if it == 0:
                ps0[mw,:] += 1.0 *ac[it,:] *decay *coswt *dt
            else:
                ps0[mw,:] += 2.0 *ac[it,:] *decay *coswt *dt

    #...Normalize in THz freq unit
    pst = np.zeros(nspcs)
    for mw in range(mwmax):
        freq = float(mw) /(tmax/1000)
        pst[:] += ps0[mw,:]*freq
    for mw in range(mwmax):
        ps0[mw,:] /= pst[:]
        #...Negative values to zero
        for ispc in range(nspcs):
            ps0[mw,ispc] = max(ps0[mw,ispc],0.0)

    ps = copy.deepcopy(ps0)
    if sgm > 0:
        for ispc in range(nspcs):
            ps[:,ispc] = gaussian_filter(ps0[:,ispc], sigma=[sgm],)
    
    with open('dat.power','w') as f:
        f.write('# Power spectrum of velocity auto correlation from vel_auto_corr.py ' +
                'at {0:s}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        f.write('#      f [THz],   I(f) of each species,   sum of species-I(f)\n')
        for mw in range(mwmax):
            freq = float(mw) /(tmax/1000)
            f.write(f' {freq:11.3e}')
            sumps = 0.0
            for ispc in range(nspcs):
                f.write(' {0:11.3e}'.format(ps[mw,ispc]))
                sumps += ps[mw,ispc]
            f.write(f' {sumps:11.3e}\n')


    # pad = lambda x: pad_zeros(x, nadd=len(x) -1)
    # w = welch(ntw)
    # t = np.array([ dt*it/1000 for it in range(ntw) ])  # [ps]
    # freqs = np.fft.fftfreq(2*ntw-1,dt/1000)[:ntw]
    # ps0= np.zeros((ntw,nspcs))
    # for ispc in range(nspcs):
    #     ps0[:,ispc] = (abs(fft(pad(ac[:,ispc])))**2)[:ntw]/tmax

    # ps = copy.deepcopy(ps0)
    # if sgm > 0:
    #     for ispc in range(nspcs):
    #         ps[:,ispc] = gaussian_filter(ps0[:,ispc], sigma=[sgm],)
    
    # with open('dat.power','w') as f:
    #     f.write('# Power spectrum of velocity auto correlation from vel_auto_corr.py ' +
    #             'at {0:s}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    #     f.write('#      f [THz],   I(f) of each species,    sum of species-I(f)\n')
    #     for it in range(ntw):
    #         f.write(f' {freqs[it]:11.3e}' )
    #         sumps = 0.0
    #         for ispc in range(nspcs):
    #             f.write(' {0:11.3e}'.format(ps[it,ispc]))
    #             sumps += ps[it,ispc]
    #         f.write(f' {sumps:11.3e}\n')

    print(' Wrote dat.autocorr and dat.power.')
    
if __name__ == "__main__":

    args= docopt(__doc__)

    main(args)
    print(' vel_auto_corr.py finished correctly ;)')

