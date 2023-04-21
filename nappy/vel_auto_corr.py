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
  --sid SID    Species-ID to be selected for power spectrum calculation.
               If it is 0, all the atoms are to be treated. [default: 0]
  -s, --shift SHIFT
               Shift of each staggered lane. [default: 20]
  -t, --time-interval TIME_INTERVAL
               Time interval between successive files in fs. [default: 1.0]
"""
import os,sys,glob,time,math
import numpy as np
from docopt import docopt

from nappy.io import read
from nappy.common import get_key

def main(args):

    #...num of measuring lane, in case of 1, it is identical to non-staggered measuring
    nmeasure = int(args['--measure'])
    
    #...shift of each staggered lane
    nshift = int(args['--shift'])
    
    #...species-ID
    sid = int(args['--sid'])
    if sid <= 0:
        sid = 0
    
    #...time interval
    dt = float(args['--time-interval'])
    
    print(' nmeasure=',nmeasure)
    print(' nshift=',nshift)
    print(' sid=',sid)
    print(' dt=',dt,'fs')
    print(' len(args)=',len(args))
    
    #...parse arguments
    infiles = args['INFILE']
    infiles.sort(key=get_key,reverse=True)
    
    #...compute sampling time-window from nmeasure and nshift
    ntwindow= len(infiles) -(nmeasure-1)*nshift
    tmax = ntwindow *dt
    print(' ntwindow =',ntwindow)
    print(' tmax     =',tmax,'fs')
    if ntwindow <= 0:
        print(' [Error] ntwindow <= 0 !!!')
        print('  Chech the parameters nmeasure and nshift, and input files.')
        sys.exit()

    #...set global values
    infile= infiles[0]

    nsys0 = read(fname=infile)
    natm= len(nsys0)
    psid= np.zeros((natm,),dtype=int)
    nas= 0
    sids = nsys0.atoms['sid']
    for ia in range(natm):
        if sid == 0:
            psid[ia] = 1
            nas += 1
        elif sids[ia] == sid:
            psid[ia] = 1
            nas += 1
    print(' num of all atoms = ',natm)
    print(' num of atoms to be considered = ',nas)
    
    print(' accumurating data',end='')
    actmp= np.zeros((nmeasure,ntwindow,3))
    acorr= np.zeros((ntwindow,3))
    v0= np.zeros((nmeasure,natm,3))
    v2= np.zeros((ntwindow,nmeasure,natm,3))
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
                    if psid[ia] == 0:
                        continue
                    v0[im,ia,:] = vels[ia,:]
            if ifile >= im*nshift and ifile-im*nshift < ntwindow:
                for ia in range(natm):
                    if psid[ia] == 0:
                        continue
                    v2[ifile-im*nshift,im,ia,:] = vels[ia,:]*v0[im,ia,:]
    print('')
    
    #.....output auto correlation function
    print(' writing dat.autocorr...')
    acfname='dat.autocorr'
    vdenom = 0.0
    for im in range(nmeasure):
        for ia in range(natm):
            if psid[ia] == 0:
                continue
            vdenom += v0[im,ia,0]**2 +v0[im,ia,1]**2 +v0[im,ia,2]**2

    acfile= open(acfname,'w')
    acfile.write('#      t [fs],   Cvv(t)x,    Cvv(t)y,    Cvv(t)z,   Cvv(t)\n')
    ac= np.zeros((ntwindow,3))
    acm= np.zeros(3)
    for it in range(ntwindow):
        for im in range(nmeasure):
            acm[:] = 0.0
            for ia in range(natm):
                if psid[ia] == 0:
                    continue
                acm[0] += v2[it,im,ia,0]
                acm[1] += v2[it,im,ia,1]
                acm[2] += v2[it,im,ia,2]
            ac[it,0] += acm[0]
            ac[it,1] += acm[1]
            ac[it,2] += acm[2]
        ac[it,0] = ac[it,0] /vdenom
        ac[it,1] = ac[it,1] /vdenom
        ac[it,2] = ac[it,2] /vdenom
        acfile.write(' {0:15.7e}'.format(it*dt) \
                     +' {0:15.7e} {1:15.7e} {2:15.7e}'.format(ac[it,0],
                                                              ac[it,1],
                                                              ac[it,2]) \
                     +' {0:15.7e}\n'.format(ac[it,0]+ac[it,1]+ac[it,2]))
    acfile.close()
    
    
if __name__ == "__main__":

    args= docopt(__doc__)

    main(args)
    print(' vel_auto_corr.py finished correctly ;)')

