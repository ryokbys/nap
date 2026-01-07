#!/usr/bin/env python
"""
Calculate the velocity auto correlation function (VACF) from pmd/dump files.

Staggered measuring of auto correlation for the statistical purpose is available with -m and -s options.

Velocity data must be stored at 1st-3rd auxiliary data of each atomic data in pmd/dump files.

Usage:
  vacf.py [options] INFILE [INFILE...]

Options:
  -h, --help   Show this help message and exit.
  -m, --measure MEASURE
               Num of measuring lane. In case of 1, it is identical to non-staggered measuring. [default: 1]
  -s, --shift SHIFT
               Shift of each staggered lane. [default: 20]
  -t, --time-interval TIME_INTERVAL
               Time interval between successive files in fs. [default: 1.0]
  --normalize
               Whether or not normalize the VACF, <v(t).v(t0)> or <v(t).v(t0)>/<v(t0).v(t0)>. [default: False]
  --sigma SIGMA
               Sigma for Gaussian smoothing by integer. [default: 0]
"""
import os,sys,glob,time,copy
from datetime import datetime
import numpy as np
from docopt import docopt

from nappy.io import read, parse_filename
from nappy.common import get_key

# from pwtools.signal import pad_zeros, welch
from scipy.fftpack import fft
from scipy.ndimage import gaussian_filter

__author__ = "Ryo KOBAYASHI"
__version__ = "260107"

def get_nsyss(infiles):
    """
    Create list of nsyss from list of files.
    """
    nsyss = []
    
    if len(infiles) > 1:
        infiles.sort(key=get_key,reverse=True)

    for ifile in range(len(infiles)):
        infile = infiles[ifile]
        print(f' Read {infile}...')
        sys.stdout.flush()
        nsystmp = read(fname=infile)
        if type(nsystmp) == list:
            nsyss.extend(nsystmp)
        else:
            nsyss.append(nsystmp)
            
    return nsyss


def pad(x, nadd=None):
    """
    与えられた信号データ x の右に len(x)-1 の長さの０で埋めた行列を加えて返す．
    """
    if nadd is None:
        nadd = len(x) -1
    add_arr = np.zeros(nadd)
    return np.concatenate([x, add_arr])


def power_spectrum_old(ntw,nspcs,ac,dt,tmax,tdamp):
    """
    vacfをpower spectrumに変換する関数（以前のやり方）．
    """
    #...Compute power spectrum
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

    return ps0, mwmax


def power_spectrum(ntw,nspcs,ac,dt,tmax,tdamp):
    """
    Convert VACF to power spectrum.
    """
    dtsec = dt * 1.0e-15  # fs to sec
    freq = np.fft.fftfreq(2*ntw-1, dtsec) / 1e+12  # Hz to THz
    ps0 = np.zeros((len(freq),nspcs))
    acd = np.zeros_like(ac, dtype=float)
    for it in range(ntw):
        t = it * dt
        decay = np.exp(-(t/tdamp)**2)
        acd[it,:] = ac[it,:] * decay
        
    for ispc in range(nspcs):
        ps0[:,ispc] = (abs(fft(pad(acd[:,ispc])))**2)
    return ps0, freq


def main(args):

    #...num of measuring lane, in case of 1, it is identical to non-staggered measuring
    nmeasure = int(args['--measure'])
    
    #...shift of each staggered lane
    nshift = int(args['--shift'])

    #...Ggaussian sigma
    sgm = int(args['--sigma'])
    
    #...time interval
    dt = float(args['--time-interval'])

    #...Normalize or not
    normalize = args['--normalize']
    
    print(' nmeasure =',nmeasure)
    print(' nshift =',nshift)
    print(' normalize =',normalize)
    print(' len(args) =',len(args))

    #...parse arguments
    infiles = args['INFILE']
    if not parse_filename(infiles[0]) in ('pmd','extxyz'):
        raise ValueError('This file format is not available.')
        
    nsyss = get_nsyss(infiles)  # nsyss is a list
    #infiles.sort(key=get_key,reverse=True)
    
    #...compute sampling time-window from nmeasure and nshift
    ntw= len(nsyss) -(nmeasure-1)*nshift
    tmax = (ntw-1) *dt
    tdamp = tmax/np.sqrt(3.0)
    print(' ntw =',ntw)
    print(f' dt, tmax = {dt:.1f}, {tmax:.1f} fs')
    fmax = 1.0/(dt*1e-15)/1e+12  # THz
    df = 1.0/(tmax*1e-15)/1e+12  # THz
    print(f' df, fmax = {df:.3f}, {fmax:.1f} THz')
    if ntw <= 0:
        print(' [Error] ntw <= 0 !!!')
        print('  Chech the parameters nmeasure and nshift, and input files.')
        sys.exit()

    #...set global values
    # infile= infiles[0]

    #nsys0 = read(fname=infile)
    nsys0 = nsyss[0]
    natm= len(nsys0)
    nspcs = len(nsys0.specorder)
    n_per_spcs = nsys0.natm_per_species()
    sids0 = nsys0.atoms['sid']
    print(' num of all atoms = ',natm)
    print(' num of atoms per species = ',end='')
    for ispc in range(nspcs):
        spi = nsys0.specorder[ispc]
        nspi= n_per_spcs[ispc]
        print(f' {spi}:{nspi},',end='')
    print('')
    
    print(' accumulating data',end='')
    actmp= np.zeros((nmeasure,ntw,3))
    acorr= np.zeros((ntw,3))
    v0= np.zeros((nmeasure,natm,3))
    v2= np.zeros((ntw,nmeasure,natm))
    vc0= np.zeros((nmeasure,nspcs,3))
    vij2= np.zeros((ntw,nmeasure,natm))
    vc2= np.zeros((ntw,nmeasure,nspcs))
    vci= np.zeros((nspcs,3))
    #for ifile in range(len(infiles)):
    for isys in range(len(nsyss)):
        # print('.',end='')
        sys.stdout.flush()
        #infile= infiles[ifile]
        #nsys = read(fname=infile)
        nsys = nsyss[isys]
        hmat= nsys.get_hmat()
        vels = nsys.get_velocities(real=True)
        
        for im in range(nmeasure):
            if isys == im*nshift:
                vci[:,:] = 0.0
                for ia in range(natm):
                    v0[im,ia,:] = vels[ia,:]
                    isp = sids0[ia] -1
                    vci[isp,:] += vels[ia,:]
                for isp in range(nspcs):
                    vc0[im,isp,:] = vci[isp,:]/n_per_spcs[isp]
            if isys >= im*nshift and isys-im*nshift < ntw:
                vci[:,:] = 0.0
                for ia in range(natm):
                    v2[isys-im*nshift,im,ia] = np.dot(vels[ia,:],v0[im,ia,:])
                    isp = sids0[ia] -1
                    vci[isp,:] += vels[ia,:]
                    # for ja in range(natm):
                    #     if ja <= ia: continue
                    #     vij2[isys-im*nshift,im,ia] += 2*np.dot(vels[ia,:],v0[im,ja,:])
                for isp in range(nspcs):
                    vci[isp,:] /= n_per_spcs[isp]
                    vc2[isys-im*nshift,im,isp] = np.dot(vci[isp,:], vc0[im,isp,:])
                    
    print('')
    
    #.....output auto correlation function
    acfname='dat.vacf'
    vdenom = np.zeros(nspcs)
    for im in range(nmeasure):
        for ia in range(natm):
            ispc = sids0[ia] -1
            vdenom[ispc] += np.dot(v0[im,ia,:],v0[im,ia,:])
    vdenom_tot = vdenom.sum()

    acfile= open(acfname,'w')
    acfile.write('# Velocity auto correlation from vacf.py ' +
                 'at {0:s}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    i = 1
    acfile.write(f'#  {i:d}:t [fs],   ')
    i += 1
    for spc in nsys0.specorder:
        acfile.write(f'{i:d}:Cvv_{spc:<2s},   ')
        i += 1
    acfile.write(f'{i:d}:Cvv_tot,  ')
    i += 1
    for spc in nsys0.specorder:
        acfile.write(f'{i:d}:COM Cvv_{spc:s}, ')
        i += 1
    acfile.write(f'{i:d}:COM Cvv_tot')
    acfile.write('\n')
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
        if normalize:
            ac[it,:] = ac[it,:] /vdenom_tot
        else:
            ac[it,:] /= nmeasure *n_per_spcs[ispc]
        acfile.write(' {0:11.3e}'.format(it*dt) )
        sumac = 0.0
        for ispc in range(nspcs):
            acfile.write(' {0:11.3e}'.format(ac[it,ispc]))
            sumac += ac[it,ispc]/nspcs
        acfile.write(f' {sumac:11.3e}')
        #...charge (center of mass) velocity
        acm[:] = 0.0
        for im in range(nmeasure):
            acm[:] += vc2[it,im,:]
        for ispc in range(nspcs):
            acm[ispc] = acm[ispc] /nmeasure *n_per_spcs[ispc]
        sumacc = 0.0
        for ispc in range(nspcs):
            acfile.write(' {0:11.3e}'.format(acm[ispc]))
            sumacc += acm[ispc]/nspcs
        acfile.write(f' {sumacc:11.3e}')
        acfile.write('\n')
    acfile.close()

    #ps0, mwmax = power_spectrum_old(ntw,nspcs,ac,dt,tmax,tdamp)
    ps0, freq = power_spectrum(ntw,nspcs,ac,dt,tmax,tdamp)

    #...Gaussian smearing
    ps = copy.deepcopy(ps0)
    if sgm > 0:
        for ispc in range(nspcs):
            ps[:,ispc] = gaussian_filter(ps0[:,ispc], sigma=[sgm],)

    #...Normalize
    if normalize:
        sumps = ps.sum()*df
        ps = ps/sumps
    
    with open('dat.ps_vacf','w') as f:
        f.write('# Power spectrum of VACF from vacf.py ' +
                'at {0:s}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        # f.write('#     f [THz],   I(f) of each species,   sum of species-I(f)\n')
        i = 1
        f.write(f'#     {i:d}:f [THz],   ')
        i += 1
        for spc in nsys0.specorder:
            f.write(f'{i:d}:I(f)_{spc:<2s},  ')
            i += 1
        f.write(f'{i:d}:I(f)_total\n')
        for ifreq in range(len(freq)//2):
            #freq = float(mw) /(tmax/1000)
            fi = freq[ifreq]
            f.write(f' {fi:11.3e}')
            sumps = 0.0
            for ispc in range(nspcs):
                f.write(' {0:11.3e}'.format(ps[ifreq,ispc]))
                sumps += ps[ifreq,ispc]
            f.write(f' {sumps:11.3e}\n')

    print(' --> dat.vacf')
    print(' --> dat.ps_vacf')
    return None

if __name__ == "__main__":

    args= docopt(__doc__)

    main(args)

