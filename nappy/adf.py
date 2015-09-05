#!/usr/bin/env python
"""
Calculate the angular distribution function (ADF) from files.
Ensemble averaging about atoms in a file and about files are taken.

Usage:
    adf.py [options] ID0 ID1 ID2 INFILE [INFILE...]

ID0, ID1, and ID2 are species-IDs consisting bonds around the angle,
like ID1-ID0-ID2.

Options:
    -h, --help  Show this help message and exit.
    -w DEG      Width of the angular degree. [default: 1.0]
    -r RCUT     Cutoff radius of the bonding pair. [default: 3.0]
    -s FMT      Input file format. If is not *akr*, users must specify it. [default: akr]
    --gsmear=SIGMA
                Width of Gaussian smearing, zero means no smearing. [default: 0]
    -p          Plot a graph on the screen. [default: False]
"""

import os,sys
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from pmdsys import pmdsys
from gaussian_smear import gsmear


def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def ndang(ia,dang,rcut,asys,id1=0,id2=0):
    """
    Compute number of atoms in the every range of angle [0:180].
    """
    na= int(180.0/dang) +1
    hmat= (asys.alc *np.array([asys.a1,asys.a2,asys.a3])).transpose()
    nda= np.zeros(na,dtype=np.int)
    natm= asys.num_atoms()
    rcut2= rcut*rcut
    pi= asys.atoms[ia].pos
    for ji in range(asys.nlspr[ia]):
        ja= asys.lspr[ia,ji]
        if ja == ia:
            continue
        if id1 != 0 and asys.atoms[ja].sid != id1:
            continue
        pj= asys.atoms[ja].pos
        pij= pj-pi
        pij= pij -np.round(pij)
        vij= np.dot(hmat,pij)
        rij2= np.dot(vij,vij)
        if rij2 >= rcut2:
            continue
        rij= np.sqrt(rij2)
        for ki in range(asys.nlspr[ia]):
            ka= asys.lspr[ia,ki]
            if ka == ia or ka <= ja:
                continue
            if id2 != 0 and asys.atoms[ka].sid != id2:
                continue
            pk= asys.atoms[ka].pos
            pik= pk-pi
            pik= pik -np.round(pik)
            vik= np.dot(hmat,pik)
            rik2= np.dot(vik,vik)
            if rik2 >= rcut2:
                continue
            rik= np.sqrt(rik2)
            cs= np.dot(vij,vik)/rij/rik
            rad= np.arccos(cs)
            deg= rad/np.pi *180.0
            nda[int(deg/dang)] += 1
    return nda

def adf(asys,dang,rcut,id0=0,id1=0,id2=0):

    natm0= asys.num_atoms()

    #.....expand system if the cell size smaller than 2*rcut
    a1norm= norm(asys.a1*asys.alc)
    a2norm= norm(asys.a1*asys.alc)
    a3norm= norm(asys.a1*asys.alc)
    print ' norm a1,a2,a3=',a1norm,a2norm,a3norm
    vol0= a1norm*a2norm*a3norm
    if a1norm/(2.0*rcut) < 1.0 or a2norm/(2.0*rcut) < 1.0 \
       or a3norm/(2.0*rcut) < 1.0:
        n1= int(2.0*rcut/a1norm)+1
        n2= int(2.0*rcut/a1norm)+1
        n3= int(2.0*rcut/a1norm)+1
        print ' system to be expanded, n1,n2,n3=',n1,n2,n3
        asys.expand(n1,n2,n3)
    print ' a1=',asys.a1
    print ' a2=',asys.a2
    print ' a3=',asys.a3

    print ' making pair list...'
    asys.make_pair_list(rcut=rcut)

    natme= asys.num_atoms()
    na= int(180.0/dang)+1
    print " rcut,dang,na=",rcut,dang,na
    anda= np.zeros(na,dtype=np.float)
    angd= np.array([ dang*ia for ia in range(na) ])
    nsum= 0
    for ia in range(natm0):
        if id0==0 or asys.atoms[ia].sid==id0:
            nsum += 1
            nda= ndang(ia,dang,rcut,asys,id1,id2)
            for iang in range(na):
                anda[iang]= anda[iang] +nda[iang]
    anda /= nsum
    return angd,anda

def adf_file_average(infiles,ffmt='akr',dang=1.0,rcut=3.0,
                     id0=0,id1=0,id2=0):
    na= int(180.0/dang) +1
    df= np.zeros(na,dtype=float)
    aadf= np.zeros(na,dtype=float)
    for infname in infiles:
        if not os.path.exists(infname):
            print "[Error] File, {0}, does not exist !!!".format(infname)
            sys.exit()
        asys= pmdsys(fname=infname,ffmt=ffmt)
        print ' infname=',infname
        angd,df= adf(asys,dang,rcut,id0,id1,id2)
        aadf += df
    aadf /= len(infiles)
    return angd,aadf

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    id0= int(args['ID0'])
    id1= int(args['ID1'])
    id2= int(args['ID2'])
    dang= float(args['-w'])
    drad= np.pi *dang/180.0
    rcut= float(args['-r'])
    flag_plot= args['-p']
    sigma= int(args['--gsmear'])
    ffmt= args['-s']

    na= int(180.0/dang) +1
    angd,agr= adf_file_average(infiles,ffmt=ffmt,dang=dang,
                               rcut=rcut,id0=id0,id1=id1,id2=id2)

    if not sigma == 0:
        angd,agr= gsmear(angd,agr,sigma)

    if flag_plot:
        plt.plot(angd, agr, '-', linewidth=1)
        plt.show()
        
    outfile= open('out.adf','w')
    for i in range(na):
        outfile.write(' {0:10.4f} {1:15.7f}\n'.format(angd[i],agr[i]))
    outfile.close()

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * out.adf'
