import sys,os,time,optparse
import numpy as np
import math,copy
from math import cos,sin,sqrt
from random import random

sys.path.append('../')
from pmdsys import PMDSystem
from atom import Atom
import cell_maker


#...constants
rcut= 2.8553
dmin= 1.5


class Grain(object):
    
    def __init__(self,point,angle):
        self.point= point
        self.angle= angle
        self.rmat= self.get_rotation_matrix(self.angle)
        
    def get_rotation_matrix(self,angle):
        a= angle[0]
        b= angle[1]
        c= angle[2]
        rmat= np.zeros((3,3))
        #...Eular angle?
        # rmat[0,0]= cos(a)*cos(b)*cos(c) -sin(a)*sin(c)
        # rmat[0,1]=-cos(a)*cos(b)*sin(c) -sin(a)*cos(c)
        # rmat[0,2]= cos(a)*sin(b)
        # rmat[1,0]= sin(a)*cos(b)*cos(c) +cos(a)*sin(c)
        # rmat[1,1]=-sin(a)*cos(b)*sin(c) +cos(a)*sin(c)
        # rmat[1,2]= sin(a)*sin(b)
        # rmat[2,0]= -sin(b)*cos(c)
        # rmat[2,1]= sin(b)*sin(c)
        # rmat[2,2]= cos(b)
        #...yaw-pitch-rolling rotation
        rmx= np.zeros((3,3))
        rmy= np.zeros((3,3))
        rmz= np.zeros((3,3))
        rmx[0,:]= [ 1.0, 0.0,     0.0 ]
        rmx[1,:]= [ 0.0, cos(a), -sin(a) ]
        rmx[2,:]= [ 0.0, sin(a),  cos(a) ]
        rmy[0,:]= [ cos(b), 0.0, sin(b) ]
        rmy[1,:]= [ 0.0,    1.0,  0.0 ]
        rmy[2,:]= [-sin(b), 0.0, cos(b) ]
        rmz[0,:]= [ cos(c), -sin(c), 0.0 ]
        rmz[1,:]= [ sin(c),  cos(c), 0.0 ]
        rmz[2,:]= [ 0.0,     0.0,    1.0 ]
        rmat= np.dot(rmx,rmy)
        rmat= np.dot(rmat,rmz)
        return rmat

def anint(x):
    if x >= 0.5:
        x = 1.0
    elif x < -0.5:
        x =-1.0
    else:
        x = 0.0
    return x

def distance(p1,p2):
    a0= p1[0]-p2[0] -anint(p1[0]-p2[0])
    a1= p1[1]-p2[1] -anint(p1[1]-p2[1])
    a2= p1[2]-p2[2] -anint(p1[2]-p2[2])
    return sqrt( a0**2 +a1**2 +a2**2 )

def pbc(x):
    if x >= 1.0:
        x -= 1.0
    elif x < 0.0:
        x += 1.0
    return x

def shift_vector():
    sv= np.zeros((27,3))
    n=0
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                sv[n,0]= float(i)
                sv[n,1]= float(j)
                sv[n,2]= float(k)
                n += 1
    return sv

def make_polycrystal(grns,uc,n1,n2,n3):
    u"""
    This routine is not that universal.
    Each grain has to have neighboring grains within a supercell,
    otherwise there will be some unexpecting grain boundries.
    In order to do so, the system should be large enough and
    the number of grains should be large enough.
    """
    sv= shift_vector()
    system= AtomSystem()
    system.set_lattice(uc.alc,uc.a1*n1,uc.a2*n2,uc.a3*n3)
    hmat= np.zeros((3,3))
    hmat[0]= system.a1 *system.alc
    hmat[1]= system.a2 *system.alc
    hmat[2]= system.a3 *system.alc
    hmati= np.linalg.inv(hmat)
    for ig in range(len(grns)):
        grain= grns[ig]
        rmat= grain.rmat
        pi= grain.point
        api= np.dot(hmat,pi)
        for ix in range(-n1/2-1,n1/2+2):
            print 'ix=',ix
            for iy in range(-n2/2-1,n2/2+2):
                for iz in range(-n3/2-1,n3/2+2):
                    for m in range(len(uc.atoms)):
                        rt= np.zeros((3,))
                        rt[0]= (uc.atoms[m].pos[0]+ix)/n1
                        rt[1]= (uc.atoms[m].pos[1]+iy)/n2
                        rt[2]= (uc.atoms[m].pos[2]+iz)/n3
                        #...rt to absolute position
                        art= np.dot(hmat,rt)
                        #...rotate
                        ari= np.dot(rmat,art)
                        ari[0]= ari[0]+api[0]
                        ari[1]= ari[1]+api[1]
                        ari[2]= ari[2]+api[2]
                        # #...reduce into the supercell
                        # ri= np.dot(hmati,ari)
                        # ri[0]= pbc(ri[0])
                        # ri[1]= pbc(ri[1])
                        # ri[2]= pbc(ri[2])
                        # ari= np.dot(hmat,ri)
                        #...check distance from all the grain points
                        di= distance(ari,api)
                        isOutside= False
                        for jg in range(len(grns)):
                            gj= grns[jg]
                            for isv in range(27):
                                pj= gj.point
                                if jg == ig and isv == 13:
                                    continue
                                svi= sv[isv]
                                pj= pj +svi
                                apj= np.dot(hmat,pj)
                                dj= distance(ari,apj)
                                if dj < di:
                                    isOutside= True
                                    break
                            if isOutside:
                                break
                        if isOutside:
                            break
                        #...here ri is inside this grain, register it
                        atom= Atom()
                        ri= np.dot(hmati,ari)
                        ri[0]= pbc(ri[0])
                        ri[1]= pbc(ri[1])
                        ri[2]= pbc(ri[2])
                        atom.set_pos(ri[0],ri[1],ri[2])
                        atom.set_sid(uc.atoms[m].sid)
                        system.add_atom(atom)
    #...remove too-close atoms at the grain boundaries
    system.make_pair_list(rcut)
    ls_remove= []
    dmin2= dmin**2
    xij= np.zeros((3,))
    for ia in range(system.num_atoms()):
        ai= system.atoms[ia]
        pi= ai.pos
        nlst= system.nlspr[ia]
        lst= system.lspr[ia]
        for j in range(nlst):
            ja= lst[j]
            aj= system.atoms[ja]
            pj= aj.pos
            xij[0]= pj[0]-pi[0] -anint(pj[0]-pi[0])
            xij[1]= pj[1]-pi[1] -anint(pj[1]-pi[1])
            xij[2]= pj[2]-pi[2] -anint(pj[2]-pi[2])
            xij= np.dot(hmat,xij)
            d2= xij[0]**2 +xij[1]**2 +xij[2]**2
            if d2 < dmin2:
                ls_remove.append(ia)
                ls_remove.append(ja)
    #...one of two will survive
    #print ls_remove
    count= [ ls_remove.count(ls_remove[i]) for i in range(len(ls_remove))]
    #print count
    for i in range(0,len(ls_remove),2):
        if count[i] > count[i+1]:
            ls_remove[i+1]= -1
        elif count[i] < count[i+1]:
            ls_remove[i]= -1
        else:
            n= int(random()*2.0) # 0 or 1
            ls_remove[i+n]= -1
    ls_remove.sort()
    for ia in range(len(ls_remove)-1,-1,-1):
        n= ls_remove[ia]
        if ia != len(ls_remove)-1:
            if n == nprev: continue
        system.atoms.pop(n)
        nprev= n
    return system

if __name__ == '__main__':
    usage= '%prog [options] nx ny nz'
    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-n",dest="ng",type="int",default=3,
                      help="number of grains in the system.")
    (options,args)= parser.parse_args()

    print '{0:=^72}'.format(' make_polycrystal_py ')
    t0= time.time()

    n1= int(args[0])
    n2= int(args[1])
    n3= int(args[2])
    ngrain= options.ng
    
    print "number of grains= {0:3d}".format(ngrain)
    print "number of unit cells= {0:3d} {1:3d} {2:3d}".format(n1,n2,n3)

#     n1= 10
#     n2= 10
#     n3= 10
#     ngrain= 3
    grains= []
    for i in range(ngrain):
        pi= np.zeros((3,))
        ai= np.zeros((3,))
        pi[0]= random()
        #pi[1]= random()
        pi[1]= 0.5
        pi[2]= random()
        # ai[0]= random()*math.pi*2 -math.pi
        # ai[1]= random()*(math.pi/2) -math.pi/2
        # ai[2]= random()*math.pi*2 -math.pi
        ai[0]= 0.0
        ai[1]= random()*math.pi*2 -math.pi
        ai[2]= 0.0
        print 'point=',pi
        print 'angle=',ai
        gi= Grain(pi,ai)
        grains.append(gi)
    uc= cell_maker.make_bcc()
    #uc.alc= 3.204
    uc.alc= 2.8553
    # uc.write_pmd('uc0000')
    system= make_polycrystal(grains,uc,n1,n2,n3)
    system.write_pmd('pmd00000')
    system.write_akr('akr0000')

    print '{0:=^72}'.format(' finished correctly ')
    print '   Elapsed time = {0:12.2f}'.format(time.time()-t0)
