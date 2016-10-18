#!/usr/bin/env python

import numpy as np

def make_pair_list(atoms,rcut=3.0,maxnn=50):
    h = atoms.get_cell()
    rc2 = rcut*rcut
    hi = np.linalg.inv(h)
    lcx = int(1.0/np.sqrt(hi[0,0]**2 +hi[1,0]**2 +hi[2,0]**2)/rcut)
    lcy = int(1.0/np.sqrt(hi[0,1]**2 +hi[1,1]**2 +hi[2,1]**2)/rcut)
    lcz = int(1.0/np.sqrt(hi[0,2]**2 +hi[1,2]**2 +hi[2,2]**2)/rcut)
    if lcx == 0: lcx= 1
    if lcy == 0: lcy= 1
    if lcz == 0: lcz= 1
    lcyz= lcy*lcz
    lcxyz= lcx*lcy*lcz
    print "lcx,y,z,xyz = ",lcx,lcy,lcz,lcxyz
    rcx= 1.0/lcx
    rcy= 1.0/lcy
    rcz= 1.0/lcz
    rcxi= 1.0/rcx
    rcyi= 1.0/rcy
    rczi= 1.0/rcz
    lscl= np.zeros((len(atoms),),dtype=int)
    lshd= np.zeros((lcxyz,),dtype=int)
    lscl[:]= -1
    lshd[:]= -1

    spos = atoms.get_scaled_positions()
    for i in range(len(spos)):
        spi= spos[i]
        #...assign a vector cell index
        mx= int(spi[0]*rcxi)
        my= int(spi[1]*rcyi)
        mz= int(spi[2]*rczi)
        m= mx*lcyz +my*lcz +mz
        # print i,pi,mx,my,mz,m
        lscl[i]= lshd[m]
        lshd[m]= i
    
    # make a pair list
    nlspr= np.zeros((len(spos)),dtype=int)
    lspr= np.zeros((len(spos),maxnn),dtype=int)
    lspr[:]= -1

    for ia in range(len(atoms)):
        pi= spos[ia]
        mx= int(pi[0]*rcxi)
        my= int(pi[1]*rcyi)
        mz= int(pi[2]*rczi)
        m= mx*lcyz +my*lcz +mz
        #print 'ia,pi,mx,my,mz,m=',ia,pi[0:3],mx,my,mz,m
        for kuz in range(-1,2):
            m1z= mz +kuz
            if m1z < 0: m1z += lcz
            if m1z >= lcz: m1z -= lcz
            for kuy in range(-1,2):
                m1y= my +kuy
                if m1y < 0: m1y += lcy
                if m1y >= lcy: m1y -= lcy
                for kux in range(-1,2):
                    m1x= mx +kux
                    if m1x < 0: m1x += lcx
                    if m1x >= lcx: m1x -= lcx
                    m1= m1x*lcyz +m1y*lcz +m1z
                    ja= lshd[m1]
                    if ja== -1: continue
                    scan_j_in_cell(ia,pi,ja,lscl,h,rc2,spos,nlspr,lspr,maxnn)
    return nlspr,lspr


def scan_j_in_cell(ia,pi,ja,lscl,h,rc2,spos,nlspr,lspr,maxnn):
    if ja == ia: ja = lscl[ja]
    if ja == -1: return 0
    
    if not ja in lspr[ia]:
        pj= spos[ja]
        xij= pj-pi
        xij= xij -np.round(xij)
        rij= np.dot(h,xij)
        rij2= rij[0]**2 +rij[1]**2 +rij[2]**2
        #print " ia,ja = ",ia,ja, np.sqrt(rij2),np.sqrt(rc2)
        if rij2 < rc2:
            n= nlspr[ia]
            lspr[ia,n]= ja
            nlspr[ia] += 1
            if nlspr[ia] >= maxnn:
                print nlspr[ia]
                print lspr[ia]
                raise RuntimeError(' [Error] nlspr[{0}] >= _maxnn !!!'.format(ia))

    ja= lscl[ja]
    scan_j_in_cell(ia,pi,ja,lscl,h,rc2,spos,nlspr,lspr,maxnn)
