#!/usr/bin/env python
"""
Calculate the radial distribution function (RDF) from *akr* files.
Ensemble averaging about atoms in a file and about files are taken.

Usage:
    rdf.py [options] IDSRC IDDST INFILE [INFILE...]

IDSRC is the species-ID of atoms to be an origin of RDF,
IDDST is the species-ID of atoms to be searched around IDSRC.
Zero for IDSRC and IDDST means any species is taken into account.

Options:
    -h, --help  Show this help message and exit.
    -d DR       Width of the bin. [default: 0.1]
    -r RMAX     Cutoff radius of radial distribution. [default: 5.0]
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

def compute_ndr(ia,dr,rmax,asys,iddst=0):
    """
    Compute number of atoms in the every shell [r:r+dr] up to *rmax*.
    This routine can be only applied to cubic systems.
    """
    nr= int(rmax/dr) +1
    ndr= [0 for n in range(nr)]
    natm= asys.num_atoms()
    ri= asys.atoms[ia].pos
    for ja in range(natm):
        if ja == ia:
            continue
        if iddst != 0 and asys.atoms[ja].sid != iddst:
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

def gr_of_file(infname,dr,rmax,idsrc=0,iddst=0):

    if not os.path.exists(infname):
        print "[Error] file does not exist !!!"
        sys.exit()
    asys= pmdsys(fname=infname)
    natm0= asys.num_atoms()

    #.....expand system if the cell size smaller than 2*rmax
    a1norm= norm(asys.a1*asys.alc)
    a2norm= norm(asys.a1*asys.alc)
    a3norm= norm(asys.a1*asys.alc)
    print ' norm a1,a2,a3=',a1norm,a2norm,a3norm
    vol0= a1norm*a2norm*a3norm
    if iddst == 0:
        rho= float(natm0)/vol0
    else:
        num_species= asys.num_species()
        rho= float(num_species[iddst-1])/vol0
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
        if idsrc==0 or asys.atoms[ia].sid==idsrc:
            ndr= compute_ndr(ia,dr,rmax,asys,iddst)
            for ir in range(nr):
                nadr[ir]= nadr[ir] +ndr[ir]
    for ir in range(1,nr):
        r= dr *ir
        nadr[ir]= float(nadr[ir])/(4.0*np.pi*rho*r*r*dr)/natm0
    return rd,nadr

def gr_file_average(infiles,dr,rmax,idsrc=0,iddst=0):
    agr= np.zeros(nr,dtype=float)
    for infname in infiles:
        print ' infname=',infname
        rd,gr= gr_of_file(infname,dr,rmax,idsrc,iddst)
        agr += gr
    agr /= len(infiles)
    return rd,agr

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    idsrc= int(args['IDSRC'])
    iddst= int(args['IDDST'])
    dr= float(args['-d'])
    rmax= float(args['-r'])
    flag_plot= args['-p']
    sigma= int(args['--gsmear'])

    nr= int(rmax/dr) +1
    rd,agr= gr_file_average(infiles,dr,rmax,idsrc,iddst)

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
