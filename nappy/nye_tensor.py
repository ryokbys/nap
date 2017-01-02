#!/usr/bin/env python
"""
Compute nye tensor around dislocation core.
The INPUT should be a LAMMPS dump file.

Usage:
  nye_tensor.py dumpconf
  nye_tensor.py calc INPUT [CONFIG]
  nye_tensor.py plot [options] [CONFIG]

Options:
  -h,--help   Show this message and exit.
  --component COMPONENT
              Specifies a component to be plot. [default: 33]
  -x X        Origin of x. [default: 0.0]
  -y Y        Origin of y. [default: 0.0]
  --width W   Half width in unit of Burgers vector. [default: 10]
  --height H  Half height in unit of Burgers vector. [default: 8]
"""
from __future__ import print_function

import os,sys,math,json
from docopt import docopt
import numpy as np
import numba

from ase.io import read

__author__ = "RYO KOBAYASHI"
__version__ = "160930"

########################################################################
#...parameters to be set to match the system
########################################################################
a0 = 4.045
c_o_a = 1.0
phi_min = 27.
iwindow = 1
xmin = 150.  # in Angstrom
xmax = 200.
ymin = 150.
ymax = 200.
nei_perf = np.zeros((12,3),dtype=float)
nei_perf[:,:] = [
    [ .5,  .5,  0.],
    [ .5, -.5,  0.],
    [-.5,  .5,  0.],
    [-.5, -.5,  0.],
    [ .5,  0.,  .5],
    [ .5,  0., -.5],
    [-.5,  0.,  .5],
    [-.5,  0., -.5],
    [ 0.,  .5,  .5],
    [ 0., -.5,  .5],
    [ 0.,  .5, -.5],
    [ 0., -.5, -.5]
]
#...rotation matrix from reference to deformed system
#...to be normalized in the script
block = np.array([[ 1.,-2., 1.],
                   [-1.,-1.,-1.],
                   [ 1., 0.,-1.]])

#...max number of nearest neighbors
maxnn = 20

def nye_tensor(atoms):
    """
    Compute Nye tensor around z-axis of given structure.
    """
    global block,maxnn,nei_perf,a0,\
        xmin,xmax,ymin,ymax,iwindow,phi_min

    phi_min = math.cos(phi_min*math.pi/180)
    pos = atoms.get_positions()
    cot = [ [i,pos[i,0],pos[i,1],pos[i,2]] for i in range(len(atoms)) ]
    co1 = np.array(cot)
    cell = atoms.get_cell()
    #...assuming that cell[2,2] is z-period
    period = cell[2,2]
    ntot = len(atoms)
    print(' num of atoms = ',ntot)

    xmin0 = np.min(co1[:,1])
    xmax0 = np.max(co1[:,1])
    ymin0 = np.min(co1[:,2])
    ymax0 = np.max(co1[:,2])
    zmin0 = np.min(co1[:,3])
    zmax0 = np.max(co1[:,3])
    zmin = zmin0 - 0.05
    zmax = zmax0 + 0.05

    for mm in range(3):
        block[mm,:] /= np.linalg.norm(block[mm,:])
    block = block.T
    print(' block after modification:')
    for mm in range(3):
        print(' {0:8.2f} {1:8.2f} {2:8.2f}'.format(block[mm,0],
                                                   block[mm,1],
                                                   block[mm,2]))
        
    nei_perf *= a0
    b = nei_perf[0,:]
    r1 = np.sqrt(np.dot(b,b))
    rc = (r1 +a0)/2
    rc2 = rc*rc
    print(' r1 = ',r1)
    print(' rc,rc2 = ',rc,rc2)

    #...Rotate the lattice vectors to the block system (what is block system?)
    # print(' nei_perf before rotation:')
    # for i in range(len(nei_perf)):
    #     print('   i,nei_perf =',i,nei_perf[i,:])
    for i in range(len(nei_perf)):
        nei_perf[i,:] = np.dot(nei_perf[i,:],block)
    # print(' nei_perf after rotation:')
    # for i in range(len(nei_perf)):
    #     print('   i,nei_perf =',i,nei_perf[i,:])

    #...select atoms within the window
    if iwindow > 0:
        tmp = []
        for i in range(ntot):
            if (xmin < co1[i,1] < xmax) and (ymin < co1[i,2] < ymax) and (zmin < co1[i,3] < zmax):
                tmp.append(co1[i,:])
        nfreew = len(tmp)
        co = np.array(tmp)

    else:
        xmin = xmin0
        xmax = xmax0
        ymin = ymin0
        ymax = ymax0
        zmin = zmax0
        zmax = zmax0
        nfreew = ntot
        co = copy.deepcopy(co1)
    print(' nfreew = ',nfreew)

    print(' >>> make_neighbors...')
    lsnn,lspr,lssft = make_neighbors_fast(co,co1,maxnn,
                                          rc,period)
    #check_neighbors(lsnn,lspr,lssft,co,co1,nei_perf)
    print(' >>> get_lattice_correspondance...')
    emat = get_lattice_correspondance(co,co1,
                                      lsnn,lspr,lssft,
                                      nei_perf,maxnn,
                                      period,r1,phi_min)
    print(' >>> compute_nye...')
    nye,inner = compute_nye(co,co1,lsnn,lspr,lssft,
                            emat,rc,period,
                            xmin,xmax,ymin,ymax,zmin,zmax)

    #...write to file
    f = open('dat.nye','w')
    for jj in range(len(inner)):
        i = inner[jj]
        f.write(' {0:8.3f} {1:8.3f} {2:8.3f}'.format(co[i,1],co[i,2],co[i,3]))
        for mm in range(3):
            f.write(' {0:8.4f} {1:8.4f} {2:8.4f}'.format(nye[jj,mm,0],
                                                        nye[jj,mm,1],
                                                        nye[jj,mm,2]))
        f.write('\n')
    f.close()
    print(' Nye tensor was written to file: dat.nye')
    return None

#@numba.jit
def make_neighbors(co,co1,maxnn,rc2,period):
    """
    Make neighbors of atoms in the window.
    Note that the neighbors can be outside the window.
    """
    #...indentify neighbors within the window
    nfreew = len(co)
    nnmax = 0
    nnmin= 20
    ntemp = np.zeros(maxnn,dtype=int)
    nshift= np.zeros(maxnn,dtype=int)
    lsnn = np.zeros(nfreew,dtype=int)
    lspr = np.zeros((nfreew,maxnn),dtype=int)
    lssft= np.zeros((nfreew,maxnn),dtype=int)
    for i in range(nfreew):
        ri = co[i,1:4]
        nn = 0
        for j in range(len(co1)):
            rj = co1[j,1:4]
            for k in (-1,0,1):
                rij = rj[:] -ri[:]
                rij[2] += period*k
                #r2 = np.dot(rij,rij)
                r2 = rij[0]*rij[0] +rij[1]*rij[1] +rij[2]*rij[2]
                if 0.001 < r2 < rc2:
                    ntemp[nn] = j
                    nshift[nn]= k
                    nn += 1
                    if nn >= maxnn:
                        raise RuntimeError('nn >= maxnn, nn =',nn)

        lsnn[i] = nn
        lspr[i,0:nn] = ntemp[0:nn]
        lssft[i,0:nn] = nshift[0:nn]
        if nn > nnmax:
            nnmax = nn
        if nn < nnmin:
            nnmin = nn
    print(' nnmax = ',nnmax)
    print(' nnmin = ',nnmin)
    return lsnn,lspr,lssft


def make_neighbors_fast(co,co1,maxnn,rc,period):
    """
    Make neighbors of atoms in the window using cell list method, 
    which is much faster than Brute-force method.
    Assuming that the PERIOD, which is z-lengh, is short enough
    not to use cell list method.
    Note that the neighbors can be outside the window.
    """
    tiny = 1.e-3
    xmin1 = np.min(co1[:,1])-tiny
    xmax1 = np.max(co1[:,1])+tiny
    ymin1 = np.min(co1[:,2])-tiny
    ymax1 = np.max(co1[:,2])+tiny
    zmin1 = np.min(co1[:,3])-tiny
    zmax1 = np.max(co1[:,3])+tiny
    wx = xmax1 - xmin1
    wy = ymax1 - ymin1
    print(' xmin1,xmax1,wx =',xmin1,xmax1,wx)
    print(' ymin1,ymax1,wy =',ymin1,ymax1,wy)
    lcx = int(wx/rc)
    lcy = int(wy/rc)
    lcxy= lcx*lcy
    print(' lcx,lcy,lcxy = ',lcx,lcy,lcxy)
    rcx = wx /lcx
    rcy = wy /lcy
    rcxi= 1.0 /rcx
    rcyi= 1.0 /rcy
    lscl= np.zeros(len(co1),dtype=int)
    lshd= np.zeros(lcxy,dtype=int)
    lscl[:] = -1
    lshd[:] = -1

    for i in range(len(co1)):
        #...assign a vector cell index
        mx = int((co1[i,1]-xmin1)*rcxi)
        my = int((co1[i,2]-ymin1)*rcyi)
        m = mx*lcy + my
        try:
            lscl[i] = lshd[m]
        except:
            print(' i,mx,my,m =',i,mx,my,m)
            print(' co1 =',co1[i,:])
            raise
        lshd[m] = i

    natm = len(co)
    natm1= len(co1)
    lsnn = np.zeros(natm,dtype=int)
    lspr = np.zeros((natm,maxnn),dtype=int)
    lssft= np.zeros((natm,maxnn),dtype=int)
    lspr[:] = -1

    rc2 = rc*rc
    for ia in range(len(co)):
        ri = co[ia,1:4]
        mx = int((ri[0]-xmin1)*rcxi)
        my = int((ri[1]-ymin1)*rcyi)
        m = mx*lcy + my
        for kuy in (-1,0,1):
            m1y = my + kuy
            if m1y < 0 or m1y >= lcy: continue
            for kux in (-1,0,1):
                m1x = mx + kux
                if m1x < 0 or m1x >= lcx: continue
                m1 = m1x*lcy + m1y
                ja = lshd[m1]
                if ja == -1: continue
                scan_j_in_cell(ia,ri,ja,lscl,rc2,co1,
                               lsnn,lspr,lssft,maxnn,
                               period)
    
    print(' nnmax = ',np.max(lsnn))
    print(' nnmin = ',np.min(lsnn))
    return lsnn,lspr,lssft


def scan_j_in_cell(ia,ri,ja,lscl,rc2,co1,
                   lsnn,lspr,lssft,maxnn,period):
    if ja == ia: ja = lscl[ja]
    if ja == -1: return 0
    
    if not ja in lspr[ia]:
        rj= co1[ja,1:4]
        for k in (-1,0,1):
            rij= rj[:] - ri[:]
            rij[2] += period*k
            r2 = rij[0]*rij[0] +rij[1]*rij[1] +rij[2]*rij[2]
            if 0.001 < r2 < rc2:
                n = lsnn[ia]
                lspr[ia,n] = ja
                lssft[ia,n]= k
                lsnn[ia] += 1
                # if lsnn[ia] > 12:
                #     print(' lsnn[ia] > 12, ia =',ia)
                #     print('   lspr = ',lspr[ia,:])
                #     print('   lssft= ',lssft[ia,:])
                if lsnn[ia] >= maxnn:
                    raise RuntimeError(' lsnn >= maxnn,'+
                                       ' ia =',ia)

    ja= lscl[ja]
    scan_j_in_cell(ia,ri,ja,lscl,rc2,co1,
                   lsnn,lspr,lssft,maxnn,period)


def check_neighbors(lsnn,lspr,lssft,co,co1,nei_perf):
    #...check neighbor list and rotated nei_perf
    print(' >>> check_neighbors...')
    idx = lsnn.argmax()
    print(' idx =',idx)
    ri = co[idx,1:4]
    print('Neighbor vectors in deformed system:')
    for jj in range(lsnn[idx]):
        j = lspr[idx,jj]
        rj = co1[j,1:4]
        rij = rj[:] -ri[:]
        rij[2] += period*lssft[idx,jj]
        print('jj,j,rij=',jj,j,rij)
    
    for j in range(len(nei_perf)):
        print('j,nei_perf =',j,nei_perf[j,:])
    return None


# @numba.jit
def get_lattice_correspondance(co,co1,lsnn,lspr,lssft,nei_perf,
                               maxnn,period,r1,phi_min):
    """
    Compute lattice correspondance between deformed and reference lattices.
    Local strain is also computed in the same loop as well.

    Returns
    --------
    emat : ndarray
           Local correspondance strain matrices for each atom.
    
    """

    nfreew = len(co)
    iden = np.identity(3)
    isym = 0
    iunique = 0
    iiu = 0
    npr = 0
    rr = np.zeros(maxnn)
    best = np.zeros(maxnn)
    km = [ 0 for i in range(maxnn) ]
    igood = np.zeros((nfreew,maxnn))
    nsym = np.zeros(nfreew)
    best_fit = np.ones(nfreew)
    emat = np.zeros((nfreew,3,3))
    angle = np.zeros(len(nei_perf))
    #...loop for atoms in deformed system
    for i in range(nfreew):
        ri = co[i,1:4]
        b = np.zeros((3,3))
        a = np.zeros((3,3))
        isym1 = isym
        iunique1 = iunique
        #...loop for neighbors of atom-i in deformed system
        for jj in range(lsnn[i]):
            j = lspr[i,jj]
            rj = co1[j,1:4]
            rij = rj[:] - ri[:]
            rij[2] += period*lssft[i,jj]
            drij= np.linalg.norm(rij)
            rr[jj] = abs(drij-r1)
            erij = rij[:]/drij
            #...check angle of bond from reference structure
            for k in range(len(nei_perf)):
                rn = nei_perf[k,0:3]
                ern = rn[:]/r1
                angle[k] = np.dot(erij,ern)
            best[jj] = np.max(angle)
            kmax = [ k for k in range(len(angle)) if abs(angle[k]-best[jj]) < 1.e-5 ]
            
            if not kmax:
                print('i,jj,j,kmax=',i,jj,j,kmax)
                print('ri,rj=',ri,rj)
                print('angle =',angle)
                raise ValueError('kmax is [] which should not happen.')
            
            if len(kmax) > 1:
                if isym == 0:
                    print(' problems with symmetry: isym==0, len(kmax)>1')
                    print('     i,jj,j =',i,jj,j)
                isym = isym + 1
                igood[i,jj] = 0
            else:
                igood[i,jj] = 1
                if best[jj] < best_fit[jj]:
                    best_fit[jj] = best[jj]
            if igood[i,jj] == 1: # if the deformed vector unique
                km[jj] = kmax
                repeat = ismember(km[0:jj],km[jj])
                if np.sum(repeat) > 0:
                    try:
                        jjj = km[0:jj].index(km[jj])
                    except:
                        print('kmax = ',kmax)
                        print('km = ',km)
                        raise
                    dr1 = rr[jjj]
                    dr2 = rr[jj]
                    #...choose deformed vector with smaller stretch
                    if dr2 > dr1:
                        iunique = iunique + 1
                        igood[i,jj] = 0
                    elif dr2 < dr1:
                        igood[i,jjj] = 0
                        iunique = iunique + 1
                    else:
                        if iiu == 0:
                            print('Lattice correspondance is not unique.')
                            iiu = 1
                        iunique = iunique + 1
            if igood[i,jj] == 1 and best[jj] < phi_min:
                igood[i,jj] = 0
            if igood[i,jj] == 1:
                p = nei_perf[kmax,:]
                b = b + np.outer(rij,p)
                a = a + np.outer(rij,rij)
        if isym > isym1:
            npr = npr + 1
            nsym[npr] = i
        try:
            ee = np.dot(np.linalg.inv(a),b) # local correspondence strain
        except:
            print('i = ',i)
            print('a  =',a)
            raise
        strn = (iden - np.dot(ee,ee)) /2 # finite local strain
        small= ( (iden - ee) + (iden -ee.T) ) /2 # local small strain
        inv1 = np.trace(small) # 1st invariant of strain
        inv2 = (np.sum(small*small) -inv1*inv1) /2 # 2nd invariant of strain
        rot = ( (iden - ee) - (iden - ee).T) /2 # local spin
        vspin = [ -rot[1,2], rot[0,2], -rot[0,1]] # local angular velocity
        emat[i,:,:] = ee[:,:]
    
    if npr > 0:
        print(' correspondence is not unique.')
        print(' number of problematic atoms =',npr)
    else:
        print(' no problematic atoms.')

    return emat

@numba.jit()
def compute_nye(co,co1,lsnn,lspr,lssft,emat,rc,period,
                xmin,xmax,ymin,ymax,zmin,zmax):
    """
    Compute Nye tensor as curl of G.

    Returns
    -------
    nye : ndarray

    """
    rcc = 1.1*rc
    # inner = [ i for i in range(len(co)) \
    #     if xmin+rcc < co[i,1] < xmax-rcc and \
    #        ymin+rcc < co[i,2] < ymax-rcc and \
    #        zmin-rcc/1000 < co[i,3] < zmax+rcc/1000 ]
    n=0
    for i in range(len(co)):
        if xmin+rcc < co[i,1] < xmax-rcc and \
           ymin+rcc < co[i,2] < ymax-rcc and \
           zmin-rcc/1000 < co[i,3] < zmax+rcc/1000:
            n += 1
    inner = np.zeros(n,dtype=int)
    n= 0
    for i in range(len(co)):
        if xmin+rcc < co[i,1] < xmax-rcc and \
           ymin+rcc < co[i,2] < ymax-rcc and \
           zmin-rcc/1000 < co[i,3] < zmax+rcc/1000:
            inner[n] = i
            n += 1
    
    mi = len(inner)
    # co2 = np.zeros((mi,4))
    # for i in range(mi):
    #     co2[i,:] = co[inner[i],0:4]
    
    grd = np.zeros((mi,3,3,3))
    for i1 in range(3):
        for j1 in range(3):
            for jj in range(mi):
                n = inner[jj]
                r = co[n,1:4]
                a = np.zeros((3,3))
                c = np.zeros((3))
                e0 = emat[n,i1,j1]
                for j in range(lsnn[n]):
                    nb = lspr[n,j]
                    nnb = find_idx_in_co(nb,co)
                    rn = co1[nb,1:4]
                    dr = rn[:] -r[:]
                    dr[2] += period*lssft[n,j]
                    a = a + np.outer(dr,dr)
                    de = emat[nnb,i1,j1] - e0
                    c = c + de*dr
                aa = np.dot(c,np.linalg.inv(a))
                grd[jj,i1,j1,:] = aa[:]
    
    #...compute curl and get Nye tensor
    nye = np.zeros((mi,3,3))
    for jj in range(mi):
        for mm in range(3):
            nye[jj,0,mm] = -grd[jj,2,mm,1] + grd[jj,1,mm,2]
            nye[jj,1,mm] = -grd[jj,0,mm,2] + grd[jj,2,mm,0]
            nye[jj,2,mm] = -grd[jj,1,mm,0] + grd[jj,0,mm,1]

    return nye,inner


#@numba.jit
def ismember(a,b):
    """
    Mimics ismember function in Matlab.
    Returns an array of the same length as A in which elements are 0 or 1.
    If an element of A exists in the array B, the element of returning array is 1,
    otherwise it is 0.
    """
    arr = np.zeros((len(a)),dtype=int)
    for i in range(len(a)):
        if a[i] in b:
            arr[i] = 1
    return arr

@numba.jit
def find_idx_in_co(idx,co):
    """
    Find the index such that co[index,0]==idx.
    """
    for i in range(len(co)):
        if co[i,0] == idx:
            return i
    return None


def dumpconf():
    global a0,c_o_a,phi_min,iwindow,xmin,xmax,ymin,ymax, \
        nei_perf,block,maxnn

    confname = 'nye_config.json'
    dic = {
        "a0": a0,
        "c_o_a": c_o_a,
        "phi_min": phi_min,
        "iwindow": iwindow,
        "xmin": xmin,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax,
        "nei_perf": nei_perf.tolist(),
        "block": block.tolist(),
        "maxnn":maxnn,
    }
    with open(confname,'w') as f:
        json.dump(dic,f,indent=2,sort_keys=True)
    print(' config file '+confname+' was written.')
    return None

def loadconf(confname):
    global a0,c_o_a,phi_min,iwindow,xmin,xmax,ymin,ymax, \
        nei_perf,block,maxnn

    with open(confname,'r') as f:
        js = json.load(f)
    a0 = js['a0']
    c_o_a = js['c_o_a']
    phi_min = js['phi_min']
    iwindow = js['iwindow']
    xmin = js['xmin']
    xmax = js['xmax']
    ymin = js['ymin']
    ymax = js['ymax']
    nei_perf = np.array(js['nei_perf'])
    block = np.array(js['block'])
    maxnn = js['maxnn']
    return None
    

def read_dat_nye(fname='dat.nye',):
    with open(fname,'r') as f:
        lines = f.readlines()
        natm = len(lines)
        pos = np.zeros((natm,3),dtype=float)
        nye = np.zeros((natm,3,3),dtype=float)
        for il,l in enumerate(lines):
            ds = l.split()
            pos[il,:] = [ float(ds[i]) for i in range(3) ]
            k = 3
            for i in range(3):
                for j in range(3):
                    nye[il,i,j] = float(ds[k])
                    k += 1
    return pos,nye
    

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['INPUT']
    confname = args['CONFIG']
    print(args)

    exit

    if args['dumpconf']:
        dumpconf()
        
    elif args['calc']:
        if confname:
            loadconf(confname)
        #...output for checking
        print(' block:')
        for i in range(3):
            print(' {0:8.2f} {1:8.2f} {2:8.2f}'.format(block[i,0],
                                                       block[i,1],
                                                       block[i,2]))
        print(' nei_perf:')
        for i in range(len(nei_perf)):
            print(' {0:4d}:'.format(i)+
                  ' {0:8.2f} {1:8.2f} {2:8.2f}'.format(nei_perf[i,0],
                                                       nei_perf[i,1],
                                                       nei_perf[i,2]))
        atoms = read(fname,format='lammps-dump')
        nye_tensor(atoms)
        
    elif args['plot']:
        #...See http://matplotlib.org/examples/pylab_examples/griddata_demo.html
        #...for contour plot example.
        from matplotlib.mlab import griddata
        import matplotlib.pyplot as plt
        if confname:
            loadconf(confname)
        comp = args['--component']
        xo = float(args['-x'])
        yo = float(args['-y'])
        w = int(args['--width'])
        h = int(args['--height'])
        if len(comp) != 2:
            raise ValueError('len(comp) != 2')
        cmpnnt = [ int(comp[0])-1, int(comp[1])-1 ]
        pos,nye = read_dat_nye(fname='dat.nye')
        x = pos[:,0] -xo
        y = pos[:,1] -yo
        z = nye[:,cmpnnt[0],cmpnnt[1]]
        nei_perf *= a0
        b = nei_perf[0,:]
        b = np.sqrt(np.dot(b,b))
        print(' b = ',b)
        print(' xmin,xmax = ',x.min(),x.max())
        print(' ymin,ymax = ',y.min(),y.max())
        xi = np.linspace(x.min(),x.max(),100)
        yi = np.linspace(y.min(),y.max(),100)
        zi = griddata(x,y,z,xi,yi,interp='linear')
        cs = plt.contourf(xi,yi,zi, 15, cmap=plt.cm.rainbow,
                          vmax=abs(zi).max(),
                          vmin=-abs(zi).max())
        plt.colorbar(shrink=float(h)/w,aspect=20.*float(h)/w)
        plt.scatter(x,y,marker='o',s=40,zorder=10,
                    facecolors='white',edgecolors='black')
        print(' -w*b,w*b = ',-w*b,w*b)
        print(' -h*b,h*b = ',-h*b,h*b)
        plt.xlim(-w*b,w*b)
        plt.ylim(-h*b,h*b)
        xtics = [ i*b for i in range(-w,w+1,1) ]  
        xticss = [ '{0:d}b'.format(i) for i in range(-w,w+1,1) ]  
        ytics = [ i*b for i in range(-h,h+1,1) ]
        yticss = [ '{0:d}b'.format(i) for i in range(-h,h+1,1) ]
        plt.xticks( xtics, xticss, fontname='sans-serif', 
                   fontsize=16 )
        plt.yticks( ytics, yticss, fontname='sans-serif', 
                   fontsize=16 )
        plt.xlabel('[-1 1 0] -->', fontname='sans-serif', 
                   fontsize=16)
        plt.ylabel('[ 1 1 1] -->', fontname='sans-serif', 
                   fontsize=16)
        plt.title('Nye component '+comp)
        plt.axes().set_aspect(1.0)
        plt.axes().tick_params(axis='x',
                               which='major',
                               pad=10)
        plt.axes().tick_params(axis='y',
                               which='major',
                               pad=10)
        plt.show()
