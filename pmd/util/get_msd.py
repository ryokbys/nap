#!/usr/local/bin/python
#
u"""
Get mean square displacements (MSDs) of given atoms (see ids below)
from the akr files in the arguments.

This utility assumes that the cell is fixed during the simulation.
"""

import os,sys,glob
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__))
                +'/../../nappy')
from AtomSystem import AtomSystem

#...atom IDs whose trajectories are tracked.
ids=(7,164,225,)

def anint(x):
    if x >= 0.5:
        return 1.0
    elif x < -0.5:
        return -1.0
    else:
        return 0.0

if len(sys.argv) < 2:
    print ' [Error] number of arguments wrong.'
    print '  Usage: ./get_msd.py akr0001 akr0002 akr0003 ...'
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

#...make output data files
outfiles= []
for l in ids:
    file= open('dat.msd-{}'.format(l),'w')
    outfiles.append(file)

p0= np.zeros((len(ids),3))
pp= np.zeros((len(ids),3))
msd= np.zeros((len(ids),3))
npbc= np.zeros((len(ids),3),dtype=int)
pit= np.zeros((3,))
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
            p0[l]= pi
            outfiles[l].write(' {:10d}'.format(ifile)
                              +' {:15.7f} {:15.7f}'.format(0.0,0.0)
                              +' {:15.7f} {:15.7f}\n'.format(0.0,0.0))
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
            # pit[0]= pi[0] -anint(dev[0])
            # pit[1]= pi[1] -anint(dev[1])
            # pit[2]= pi[2] -anint(dev[2])
            #...normalized to absolute
            #dev= pit -pp[l]
            dev[0]= pi[0] -p0[l,0] +float(npbc[l,0])
            dev[1]= pi[1] -p0[l,1] +float(npbc[l,1])
            dev[2]= pi[2] -p0[l,2] +float(npbc[l,2])
            dev= np.dot(hmat,dev)
            # msd[l,0] += dev[0]**2 
            # msd[l,1] += dev[1]**2 
            # msd[l,2] += dev[2]**2 
            msd[l,0] = dev[0]**2 
            msd[l,1] = dev[1]**2 
            msd[l,2] = dev[2]**2 
            
            outfiles[l].write(' {:10d}'.format(ifile)
                              +' {:15.7f}'.format(msd[l,0])
                              +' {:15.7f}'.format(msd[l,1])
                              +' {:15.7f}'.format(msd[l,2])
                              +' {:15.7f}'.format((msd[l,0]
                                                   +msd[l,1]
                                                   +msd[l,2])/3)
                              +' \n')
            #...store current position
            pp[l]= pi
            
