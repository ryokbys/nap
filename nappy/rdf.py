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
    -s FMT      Input file format. If is not *akr*, users must specify it. [default: akr]
    --gsmear=SIGMA
                Width of Gaussian smearing, zero means no smearing. [default: 0]
    -o OUT      Output file name. [default: out.rdf]
    -p          Plot a graph on the screen. [default: False]
"""

import os,sys
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from pmdsys import PMDSystem
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
    ndr= np.zeros(nr,dtype=np.int)
    natm= asys.num_atoms()
    hmat= (asys.alc *np.array([asys.a1,asys.a2,asys.a3])).transpose()
    pi= asys.atoms[ia].pos
    for ja in range(natm):
        if ja == ia:
            continue
        if iddst != 0 and asys.atoms[ja].sid != iddst:
            continue
        pj= asys.atoms[ja].pos
        pij= pj -pi
        pij= pij -np.round(pij)
        vij= np.dot(hmat,pij)
        rij2= np.dot(vij,vij)
        rij= np.sqrt(rij2)
        #print ia,ja,rx,ry,rz,rr,r
        if rij >= rmax+dr/2:
            continue
        rrdr= rij/dr
        ir= int(round(rrdr))
        #print "ia,ja,rr,rx,ry,rz,ir=",ia,ja,rr,rx,ry,rz,ir
        ndr[ir]= ndr[ir] +1
        #print 'ia,ja,rij,ir,ndr[ir]=',ia,ja,rij,ir,ndr[ir]
    return ndr

def rdf(asys,dr,rmax,idsrc=0,iddst=0):

    natm0= asys.num_atoms()
    vol= asys.volume()
    ns= asys.num_species()
    nsrc= asys.num_atoms(sid=idsrc)
    ndst= asys.num_atoms(sid=iddst)
    rho= float(ndst)/vol
    print ' natm0,nsrc,ndst,vol,rho=',natm0,nsrc,ndst,vol,rho

    n1,n2,n3= asys.get_expansion_num(2.0*rmax)
    if not (n1==1 and n2==1 and n3==1):
        print ' system to be expanded, n1,n2,n3=',n1,n2,n3
        asys.expand(n1,n2,n3)
    # print ' a1=',asys.a1
    # print ' a2=',asys.a2
    # print ' a3=',asys.a3

    nr= int(rmax/dr)+1
    # print " rmax,dr,nr=",rmax,dr,nr
    nadr= [0.0 for n in range(nr)]
    rd= [ dr*ir for ir in range(nr) ]
    for ia in range(natm0):
        if idsrc==0 or asys.atoms[ia].sid==idsrc:
            ndr= compute_ndr(ia,dr,rmax,asys,iddst)
            for ir in range(nr):
                nadr[ir]= nadr[ir] +ndr[ir]
    #...normalize
    for ir in range(1,nr):
        r= dr *ir
        nadr[ir]= float(nadr[ir])/(4.0*np.pi*rho*r*r*dr)
    return rd,nadr,nsrc

def rdf_average(infiles,ffmt='akr',dr=0.1,rmax=3.0,
                idsrc=0,iddst=0):
    agr= np.zeros(nr,dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print "[Error] File, {0}, does not exist !!!".format(infname)
            sys.exit()
        asys= PMDSystem(fname=infname,ffmt=ffmt)
        print ' infname=',infname
        rd,gr,n= rdf(asys,dr,rmax,idsrc,iddst)
        nsum += n
        agr += gr
    # agr /= len(infiles)
    #print ' nsum=',nsum
    agr /= nsum
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
    ffmt= args['-s']
    ofname= args['-o']

    nr= int(rmax/dr) +1
    rd,agr= rdf_average(infiles,ffmt=ffmt,dr=dr,rmax=rmax,
                        idsrc=idsrc,iddst=iddst)

    if not sigma == 0:
        rd,agr= gsmear(rd,agr,sigma)

    if flag_plot:
        plt.plot(rd, agr, '-', linewidth=1)
        plt.show()
        
    outfile= open(ofname,'w')
    for i in range(nr):
        outfile.write(' {0:10.4f} {1:15.7f}\n'.format(rd[i],agr[i]))
    outfile.close()
