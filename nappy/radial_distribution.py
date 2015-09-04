#!/bin/env python
"""
Calculate the radial distribution function from *akr* files.
Ensemble averaging about atoms in a file and about files are taken.

Usage:
    radial_distribution.py [options] INFILE [INFILE...]

Options:
    -h, --help  Show this help message and exit.
    -d DR       Width of the bin. [default: 0.1]
    -r RMAX     Cutoff radius of radial distribution. [default: 5.0]
    --src-sid=IDSRC
                Species-ID of source atoms to be selected.
                Zero means all the species will be selected. [default: 0]
    --dst-sid=IDDST
                Species-ID of source atoms to be selected.
                Zero means all the species will be selected. [default: 0]
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

def read_akr(fname):
    h=np.zeros((3,3),dtype=float)
    f=open(fname,'r')
    hunit= float(f.readline())
    h[0,0:3]= [ float(x) for x in f.readline().split() ]
    h[1,0:3]= [ float(x) for x in f.readline().split() ]
    h[2,0:3]= [ float(x) for x in f.readline().split() ]
    natm= int(f.readline().split()[0])
    isp= np.zeros((natm,),dtype=int)
    ra= np.zeros((natm,3),dtype=float)
    for ia in range(natm):
        data= f.readline().split()
        isp[ia]= int(data[0])
        ra[ia,0]= float(data[1])
        ra[ia,1]= float(data[2])
        ra[ia,2]= float(data[3])
    f.close()
    return hunit,h,natm,isp,ra

def compute_ndr(ia,dr,rmax,asys):
    """
    This routine can be only applied to cubic systems.
    """
    nr= int(rmax/dr) +1
    ndr= [0 for n in range(nr)]
    natm= asys.num_atoms()
    ri= asys.atoms[ia].pos
    for ja in range(natm):
        if ja == ia:
            continue
        r1= asys.atoms[ja].pos[0] -ri[0]
        r2= asys.atoms[ja].pos[1] -ri[1]
        r3= asys.atoms[ja].pos[2] -ri[2]
        #print ia,ja,ra[ia][0:3],ra[ja][0:3]
        r1= (r1 -round(r1))
        r2= (r2 -round(r2))
        r3= (r3 -round(r3))
        rx= asys.alc*(r1*asys.a1[0] +r2*asys.a2[0] +r3*asys.a3[0])
        ry= asys.alc*(r1*asys.a1[1] +r2*asys.a2[1] +r3*asys.a3[1])
        rz= asys.alc*(r1*asys.a1[2] +r2*asys.a2[2] +r3*asys.a3[2])
        rr2= rx*rx +ry*ry +rz*rz
        rr= np.sqrt(rr2)
        #print ia,ja,rx,ry,rz,rr,r
        if rr >= rmax+dr/2:
            continue
        rrdr= rr/dr
        ir= int(round(rrdr))
        #print "ia,ja,rr,rx,ry,rz,ir=",ia,ja,rr,rx,ry,rz,ir
        ndr[ir]= ndr[ir] +1
    return ndr

def gr_of_file(infname,dr,rmax):

    if not os.path.exists(infname):
        print "[Error] file does not exist !!!"
        exit()
    asys= pmdsys()
    asys.read_akr(infname)
    natm0= asys.num_atoms()

    #.....expand system if the cell size smaller than 2*rmax
    a1norm= norm(asys.a1*asys.alc)
    a2norm= norm(asys.a1*asys.alc)
    a3norm= norm(asys.a1*asys.alc)
    print ' norm a1,a2,a3=',a1norm,a2norm,a3norm
    vol0= a1norm*a2norm*a3norm
    rho= float(natm0)/vol0
    print " natm,vol,rho=",natm0,vol0,rho
    if a1norm/(2.0*rmax) < 1.0 or a2norm/(2.0*rmax) < 1.0 \
       or a3norm/(2.0*rmax) < 1.0:
        n1= int(2.0*rmax/a1norm)+1
        n2= int(2.0*rmax/a1norm)+1
        n3= int(2.0*rmax/a1norm)+1
        print ' system to be expanded, n1,n2,n3=',n1,n2,n3
        asys.expand(n1,n2,n3)
    print ' a1=',asys.a1
    print ' a2=',asys.a2
    print ' a3=',asys.a3

    natme= asys.num_atoms()
    nr= int(rmax/dr)+1
    print " rmax,dr,nr=",rmax,dr,nr
    nadr= [0.0 for n in range(nr)]
    rd= [ dr*ir for ir in range(nr) ]
    for ia in range(natm0):
        #print "ia=",ia
        ndr= compute_ndr(ia,dr,rmax,asys)
        for ir in range(nr):
            nadr[ir]= nadr[ir] +ndr[ir]
    for ir in range(1,nr):
        r= dr *ir
        nadr[ir]= float(nadr[ir])/(4.0*np.pi*rho*r*r*dr)/natm0
    return rd,nadr

def gr_file_average(infiles,dr,rmax):
    agr= np.zeros(nr,dtype=float)
    for infname in infiles:
        print ' infname=',infname
        rd,gr= gr_of_file(infname,dr,rmax)
        agr += gr
    agr /= len(infiles)
    return rd,agr

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    idsrc= int(args['--src-sid'])
    iddst= int(args['--dst-sid'])
    dr= float(args['-d'])
    rmax= float(args['-r'])
    flag_plot= args['-p']
    sigma= int(args['--gsmear'])

    nr= int(rmax/dr) +1
    rd,agr= gr_file_average(infiles,dr,rmax)

    if not sigma == 0:
        rd,agr= gsmear(rd,agr,sigma)

    if flag_plot:
        plt.plot(rd, agr, '-', linewidth=1)
        plt.show()
        
    outfile= open('out.radial_distribution','w')
    for i in range(nr):
        outfile.write(' {0:10.4f} {1:15.7f}\n'.format(rd[i],agr[i]))
    outfile.close()

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * out.radial_distribution'
