#!/usr/local/bin/python
#
u"""
Get mean square displacements (MSDs) of given atoms (see ids below)
from the akr files in the arguments.

This utility assumes that the cell is fixed during the simulation.

Staggered measuring of MSD for for the statistical purpose.
"""

import os,sys,glob,time
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__))
                +'/../../nappy')
from AtomSystem import AtomSystem

#...atom IDs whose trajectories are tracked.
ids=(7,164,225,)
#...num of measuring lane, in case of 1, it is identical to non-staggered measuring
nmeasure= 10
#...shift of each staggered lane
nshift= 20

def anint(x):
    if x >= 0.5:
        return 1.0
    elif x < -0.5:
        return -1.0
    else:
        return 0.0

if len(sys.argv) < 2:
    print ' [Error] number of arguments wrong.'
    print '  Usage: ./get_msd_staggered.py akr0001 akr0002 akr0003 ...'
    sys.exit()

#...parse arguments
infiles= []
#print os.getcwd()
for i in range(len(sys.argv)):
    if i==0: continue
    if sys.argv[i].find('*') >= 0 or sys.argv[i].find('?') >= 0:
        files= glob.glob(sys.argv[i])
        for f in files:
            infiles.append(f)
    infiles.append(sys.argv[i])

#...compute sampling time-window from nmeasure and nshift
ntwindow= len(infiles) -(nmeasure-1)*nshift
if ntwindow <= 0:
    print ' [Error] ntwindow <= 0 !!!'
    print '  Chech the parameters nmeasure and nshift, and input files.'
    sys.exit()


#...make output data files
outfiles= []
for l in ids:
    file= open('dat.msd-{}'.format(l),'w')
    outfiles.append(file)

p0= np.zeros((nmeasure,len(ids),3))
pp= np.zeros((len(ids),3))
msd= np.zeros((len(infiles),nmeasure,len(ids),3))
npbc= np.zeros((len(ids),3),dtype=int)
hmat= np.zeros((3,3))
for ifile in range(len(infiles)):
    file= infiles[ifile]
    system= AtomSystem()
    system.read_akr(file)
    hmat[0]= system.a1 *system.alc
    hmat[1]= system.a2 *system.alc
    hmat[2]= system.a3 *system.alc
    for l in range(len(ids)):
        #...human-readable ID to computer-oriented ID
        i= ids[l]-1
        ai= system.atoms[i]
        pi= ai.pos
        if ifile == 0:
            pp[l]= pi
        else:
            #...correct periodic motion
            dev= pi -pp[l]
            if dev[0] > 0.5:
                npbc[l,0] += -1
            elif dev[0] < -0.5:
                npbc[l,0] += 1
            if dev[1] > 0.5:
                npbc[l,1] += -1
            elif dev[1] < -0.5:
                npbc[l,1] += 1
            if dev[2] > 0.5:
                npbc[l,2] += -1
            elif dev[2] < -0.5:
                npbc[l,2] += 1
            # print npbc
            #...store current position
            pp[l]= pi
                        
        for nm in range(nmeasure):
            if ifile == nm*nshift:
                p0[nm,l,0]= pi[0] +float(npbc[l,0])
                p0[nm,l,1]= pi[1] +float(npbc[l,1])
                p0[nm,l,2]= pi[2] +float(npbc[l,2])
            if nm*nshift < ifile:
                #...normalized to absolute
                dev[0]= pi[0] -p0[nm,l,0] +float(npbc[l,0])
                dev[1]= pi[1] -p0[nm,l,1] +float(npbc[l,1])
                dev[2]= pi[2] -p0[nm,l,2] +float(npbc[l,2])
                dev= np.dot(hmat.T,dev)
                msd[ifile-nm*nshift,nm,l,0] = dev[0]**2 
                msd[ifile-nm*nshift,nm,l,1] = dev[1]**2 
                msd[ifile-nm*nshift,nm,l,2] = dev[2]**2 
                

for ifile in range(len(infiles)-(nmeasure-1)*nshift):
    for l in range(len(ids)):
        if ifile == 0:
            outfiles[l].write(' {:10d}'.format(ifile)
                              +' {:15.7f} {:15.7f}'.format(0.0,0.0)
                              +' {:15.7f} {:15.7f}\n'.format(0.0,0.0))
        else:
            dev= np.zeros((3,))
            for nm in range(nmeasure):
                dev += msd[ifile,nm,l]
            dev /= nmeasure
            outfiles[l].write(' {:10d}'.format(ifile)
                              +' {:15.7f}'.format(dev[0])
                              +' {:15.7f}'.format(dev[1])
                              +' {:15.7f}'.format(dev[2])
                              +' {:15.7f}'.format((dev[0]+dev[1]+dev[2])/3)
                              +' \n')
