#!/bin/env python
#
u"""
Get mean square displacements (MSDs) of given atoms (see ids below)
from the akr files in the arguments.

This utility assumes that the cell is fixed during the simulation.

Staggered measuring of MSD for for the statistical purpose.
"""
from __future__ import print_function

import os,sys,glob,time
import numpy as np
import optparse

from nappy.napsys import NAPSystem

usage= '%prog [options] akr0001 akr0002 akr0003 ...'

parser= optparse.OptionParser(usage=usage)
parser.add_option("--id",dest="id",type="int",default=1,
                  help="id of an atom whose path is to be traced.")
parser.add_option("-m","--measure",dest="measure",type="int",default=1,
                  help="num of measuring lane. In case of 1, it is identical to non-staggered measuring.")
parser.add_option("-s","--shift",dest="shift",type="int",default=20,
                  help="shift of each staggered lane.")
(options,args)= parser.parse_args()

#...atom IDs whose trajectories are tracked.
#ids=(55,)
id= options.id
#...num of measuring lane, in case of 1, it is identical to non-staggered measuring
#nmeasure= 10
nmeasure= options.measure
#...shift of each staggered lane
#nshift= 50
nshift= options.shift

# print id
# print nmeasure
# print nshift
# print args
# exit()

def anint(x):
    if x >= 0.5:
        return 1.0
    elif x < -0.5:
        return -1.0
    else:
        return 0.0

if len(args) < 1:
    print(' [Error] number of arguments wrong.')
    print(usage)
    sys.exit()

#...parse arguments
infiles= []
#print os.getcwd()
for file in args:
    infiles.append(file)

#...compute sampling time-window from nmeasure and nshift
ntwindow= len(infiles) -(nmeasure-1)*nshift
if ntwindow <= 0:
    print(' [Error] ntwindow <= 0 !!!')
    print('  Chech the parameters nmeasure and nshift, and input files.')
    sys.exit()


#...make output data files
outfname='dat.msd-{0}'.format(id)
outfile= open(outfname,'w')

p0= np.zeros((nmeasure,3))
pp= np.zeros(3)
msd= np.zeros((len(infiles),nmeasure,3))
npbc= np.zeros((3,),dtype=int)
hmat= np.zeros((3,3))
for ifile in range(len(infiles)):
    file= infiles[ifile]
    system= NAPSystem()
    system.read_akr(file)
    hmat[0]= system.a1 *system.alc
    hmat[1]= system.a2 *system.alc
    hmat[2]= system.a3 *system.alc
    #...human-readable ID to computer-oriented ID
    i= id-1
    ai= system.atoms[i]
    pi= ai.pos
    if ifile == 0:
        pp= pi
    else:
        #...correct periodic motion
        dev= pi -pp
        if dev[0] > 0.5:
            npbc[0] += -1
        elif dev[0] < -0.5:
            npbc[0] += 1
        if dev[1] > 0.5:
            npbc[1] += -1
        elif dev[1] < -0.5:
            npbc[1] += 1
        if dev[2] > 0.5:
            npbc[2] += -1
        elif dev[2] < -0.5:
            npbc[2] += 1
        # print npbc
        #...store current position
        pp= pi
                        
    for nm in range(nmeasure):
        if ifile == nm*nshift:
            p0[nm,0]= pi[0] +float(npbc[0])
            p0[nm,1]= pi[1] +float(npbc[1])
            p0[nm,2]= pi[2] +float(npbc[2])
        if nm*nshift < ifile:
            #...normalized to absolute
            dev[0]= pi[0] -p0[nm,0] +float(npbc[0])
            dev[1]= pi[1] -p0[nm,1] +float(npbc[1])
            dev[2]= pi[2] -p0[nm,2] +float(npbc[2])
            dev= np.dot(hmat.T,dev)
            msd[ifile-nm*nshift,nm,0] = dev[0]**2 
            msd[ifile-nm*nshift,nm,1] = dev[1]**2 
            msd[ifile-nm*nshift,nm,2] = dev[2]**2 
                

for ifile in range(len(infiles)-(nmeasure-1)*nshift):
    if ifile == 0:
        outfile.write(' {:10d}'.format(ifile)
                      +' {:15.7f} {:15.7f}'.format(0.0,0.0)
                      +' {:15.7f} {:15.7f}\n'.format(0.0,0.0))
    else:
        dev= np.zeros((3,))
        for nm in range(nmeasure):
            dev += msd[ifile,nm]
        dev /= nmeasure
        outfile.write(' {:10d}'.format(ifile)
                      +' {:15.7f}'.format(dev[0])
                      +' {:15.7f}'.format(dev[1])
                      +' {:15.7f}'.format(dev[2])
                      +' {:15.7f}'.format((dev[0]+dev[1]+dev[2])/3)
                      +' \n')
