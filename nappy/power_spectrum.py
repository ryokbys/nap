#!/bin/env python
#
u"""
Calculate the velocity auto correlation function and power spectrum
from akr files.

Staggered measuring of auto correlation for the statistical purpose.

Velocity data must be stored at 1st-3rd auxiliary data of each atomic
data in akr files.
"""

import os,sys,glob,time,math
import numpy as np
import optparse

from pmdsys import pmdsys

################################################## Functions ###########


############################################### Main routine ###########

usage= '%prog [options] akr0001 akr0002 akr0003 ...'

parser= optparse.OptionParser(usage=usage)
parser.add_option("--id",dest="id",type="int",default=1,
                  help="id of an atom whose path is to be traced.")
parser.add_option("-f","--file",dest="id_file",type="string",default=None,
                  help="file name which contains the list of IDs of atoms whose paths are to be traced.")
parser.add_option("-m","--measure",dest="measure",type="int",default=1,
                  help="num of measuring lane. In case of 1, it is identical to non-staggered measuring.")
parser.add_option("--sid",dest="sid",type="int",default=0,
                  help="species-ID to be selected for power spectrum calculation.")
parser.add_option("-s","--shift",dest="shift",type="int",default=20,
                  help="shift of each staggered lane.")
parser.add_option("-t","--time-interval",dest="tval",type="float",
                  default=1.0, 
                  help="time interval between successive akr files in [fs].")
parser.add_option("--relax",dest="trlx",type="float",
                  default=1e+10, 
                  help="relaxation time [fs] of decaying factor.")
(options,args)= parser.parse_args()

#...file name of the ID list
idfile= options.id_file
#...atom IDs whose trajectories are tracked.
if idfile is None:
    id= options.id
#...num of measuring lane, in case of 1, it is identical to non-staggered measuring
nmeasure= options.measure
#...shift of each staggered lane
nshift= options.shift
#...species-ID
sid= options.sid
#...time interval
tval= options.tval
trlx= options.trlx

print ' idfile=',idfile
print ' id=',id
print ' nmeasure=',nmeasure
print ' nshift=',nshift
print ' sid=',sid
print ' tval=',tval
print ' len(args)=',len(args)
# sys.exit()

if len(args) < 1:
    print ' [Error] number of arguments wrong.'
    print usage
    sys.exit()

#...parse arguments
infiles= []
for file in args:
    infiles.append(file)
infiles.sort()

#...compute sampling time-window from nmeasure and nshift
ntwindow= len(infiles) -(nmeasure-1)*nshift
print ' ntwindow=',ntwindow
if ntwindow <= 0:
    print ' [Error] ntwindow <= 0 !!!'
    print '  Chech the parameters nmeasure and nshift, and input files.'
    sys.exit()


#...set global values
infile= infiles[0]
system=pmdsys()
system.read_akr(infile)
natm= len(system.atoms)
psid= np.zeros((natm,),dtype=np.int8)
nas= 0
for ia in range(natm):
    if sid == 0:
        psid[ia] = 1
        nas += 1
    elif system.atoms[ia].sid == sid:
        psid[ia] = 1
        nas += 1
print ' num of all atoms = ',natm
print ' num of atoms to be considered = ',nas

print ' accumurating data',
actmp= np.zeros((nmeasure,ntwindow,3))
acorr= np.zeros((ntwindow,3))
v02= np.zeros((nmeasure,natm,3))
v2= np.zeros((ntwindow,nmeasure,natm,3))
hmat= np.zeros((3,3))
for ifile in range(len(infiles)):
    print '.',
    sys.stdout.flush()
    infile= infiles[ifile]
    system= pmdsys()
    system.read_akr(infile)
    hmat[0]= system.a1 *system.alc
    hmat[1]= system.a2 *system.alc
    hmat[2]= system.a3 *system.alc
    atoms= system.atoms

    for im in range(nmeasure):
        if ifile == im*nshift:
            for ia in range(natm):
                if psid[ia] == 0:
                    continue
                vi=atoms[ia].vel
                v02[im,ia,0] = vi[0]
                v02[im,ia,1] = vi[1]
                v02[im,ia,2] = vi[2]
        if ifile >= im*nshift and ifile-im*nshift < ntwindow:
            for ia in range(natm):
                if psid[ia] == 0:
                    continue
                vi= atoms[ia].vel
                v2[ifile-im*nshift,im,ia,0] = vi[0]*v02[im,ia,0]
                v2[ifile-im*nshift,im,ia,1] = vi[1]*v02[im,ia,1]
                v2[ifile-im*nshift,im,ia,2] = vi[2]*v02[im,ia,2]
print ' '

#.....output auto correlation function
print ' writing dat.autocorr...'
acfname='dat.autocorr'
acfile= open(acfname,'w')
acfile.write('#      t [fs],   Cvv(t)x,    Cvv(t)y,    Cvv(t)z,   Cvv(t)\n')
ac= np.zeros((len(infiles)-(nmeasure-1)*nshift,3))
for it in range(len(infiles)-(nmeasure-1)*nshift):
    for im in range(nmeasure):
        acm= np.zeros(3)
        for ia in range(natm):
            if psid[ia] == 0:
                continue
            vdenom= v02[im,ia,0]**2 \
                    +v02[im,ia,1]**2 \
                    +v02[im,ia,2]**2
            acm[0] += v2[it,im,ia,0]/vdenom
            acm[1] += v2[it,im,ia,1]/vdenom
            acm[2] += v2[it,im,ia,2]/vdenom
        ac[it,0] += acm[0]/nas
        ac[it,1] += acm[1]/nas
        ac[it,2] += acm[2]/nas
    decay= np.exp(-it*tval/trlx)
    ac[it,0] = ac[it,0]/nmeasure *decay
    ac[it,1] = ac[it,1]/nmeasure *decay
    ac[it,2] = ac[it,2]/nmeasure *decay
    acfile.write(' {0:10.1f}'.format(it*tval) \
                 +' {0:15.7f}'.format(ac[it,0]) \
                 +' {0:15.7f}'.format(ac[it,1]) \
                 +' {0:15.7f}'.format(ac[it,2]) \
                 +' {0:15.7f}\n'.format(ac[it,0]+ac[it,1]+ac[it,2]))
acfile.close()

#.....output power spectrum
print ' writing dat.power...'
psfname='dat.power'
psfile= open(psfname,'w')
mmax= ntwindow/2
psfile.write('#     w[cm^-1],      Ix(w),     Iy(w),     Iz(w),     I(w)\n')
speed= 299792458.0*100 /1.0e+15 # cm/fs
for mw in range(mmax):
    ps = np.zeros(3)
    for it in range(ntwindow):
        coswt= np.cos(2.0*np.pi*mw*it/ntwindow)
        ps[:] += ac[it,:] *coswt
    ps[:] = ps[:] *tval /np.pi
    psfile.write(' {0:15.2f}'.format(float(mw)/(ntwindow*tval)/speed) \
                 +' {0:15.7f}'.format(ps[0]) \
                 +' {0:15.7f}'.format(ps[1]) \
                 +' {0:15.7f}'.format(ps[2]) \
                 +' {0:15.7f}\n'.format(ps[0]+ps[1]+ps[2]) )
psfile.close()
print ' power_spectrum.py finished correctly ;)'
