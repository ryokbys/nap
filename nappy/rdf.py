#!/usr/bin/env python
"""
Calculate the radial distribution function (RDF) from files.
Ensemble averaging about atoms in a file and about files are taken.

Usage:
    rdf.py [options] INFILE [INFILE...]

Options:
    -h, --help  Show this help message and exit.
    -d DR       Width of the bin. [default: 0.1]
    -r RMAX     Cutoff radius of radial distribution. [default: 5.0]
    -s FMT      Input file format. If is not *pmd*, users must specify it. [default: pmd]
    --gsmear=SIGMA
                Width of Gaussian smearing, zero means no smearing. [default: 0]
    -o OUT      Output file name. [default: out.rdf]
    --num-species=NSPCS
                Number of species in the system. [default: 1]
"""

import os,sys
import numpy as np
from docopt import docopt
from napsys import NAPSystem
from gaussian_smear import gsmear


def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def compute_ndr(ia,isid,dr,rmax,asys,nspcs):
    """
    Compute number of atoms in the every shell [r:r+dr] up to *rmax*.
    This routine can be only applied to cubic systems.
    """
    nr= int(rmax/dr) +1
    #ndr= np.zeros(nr,dtype=np.int)
    ndr = np.zeros((nspcs+1,nspcs+1,nr),dtype=np.int)
    natm= asys.num_atoms()
    hmat= (asys.alc *np.array([asys.a1,asys.a2,asys.a3])).transpose()
    pi= asys.atoms[ia].pos
    for ja in range(natm):
        if ja == ia:
            continue
        jsid = asys.atoms[ja].sid
        pj= asys.atoms[ja].pos
        pij= pj -pi
        pij= pij -np.round(pij)
        vij= np.dot(hmat,pij)
        rij2= np.dot(vij,vij)
        rij= np.sqrt(rij2)
        if rij >= rmax:
            continue
        rrdr= rij/dr
        ir = int(rrdr)
        ndr[0,0,ir]= ndr[0,0,ir] +1
        ndr[isid,jsid,ir] = ndr[isid,jsid,ir] +1
    return ndr

def rdf(asys,nspcs,dr,rmax):

    natm0= asys.num_atoms()
    vol= asys.volume()
    rho= float(natm0)/vol
    print ' natm0,vol,rho=',natm0,vol,rho

    n1,n2,n3= asys.get_expansion_num(2.0*rmax)
    if not (n1==1 and n2==1 and n3==1):
        print ' system to be repeated, n1,n2,n3=',n1,n2,n3
        asys.repeat(n1,n2,n3)
    # print ' a1=',asys.a1
    # print ' a2=',asys.a2
    # print ' a3=',asys.a3

    nr= int(rmax/dr)+1
    # print " rmax,dr,nr=",rmax,dr,nr
    nadr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    rd= [ dr*ir+dr/2 for ir in range(nr) ]
    for ia in range(natm0):
        isid = asys.atoms[ia].sid
        ndr= compute_ndr(ia,isid,dr,rmax,asys,nspcs)
        for ir in range(nr):
            nadr[:,:,ir] += ndr[:,:,ir]
    #nadr /= nsrc
    #print nadr
    #...normalize
    for ir in range(1,nr):
        r= dr *ir
        nadr[:,:,ir]= nadr[:,:,ir]/(4.0*np.pi*rho*r*r*dr)
    #print nadr
    return rd,nadr,natm0

def rdf_average(infiles,nspcs,nr,ffmt='akr',dr=0.1,rmax=3.0,):
    agr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print "[Error] File, {0}, does not exist !!!".format(infname)
            sys.exit()
        asys= NAPSystem(fname=infname,ffmt=ffmt)
        print ' infname=',infname
        rd,gr,n= rdf(asys,nspcs,dr,rmax)
        nsum += n
        agr += gr
    # agr /= len(infiles)
    #print ' nsum=',nsum
    #print agr
    agr /= nsum
    return rd,agr

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    dr= float(args['-d'])
    rmax= float(args['-r'])
    sigma= int(args['--gsmear'])
    ffmt= args['-s']
    ofname= args['-o']
    nspcs = int(args['--num-species'])

    nr= int(rmax/dr) +1
    rd,agr= rdf_average(infiles,nspcs,nr,ffmt=ffmt,dr=dr,rmax=rmax,)

    if not sigma == 0:
        print ' Gaussian smearing...'
        #...Smearing of total RDF
        agrt= gsmear(rd,agr[0,0],sigma)
        agr[0,0,:] = agrt[:]
        #...Smearing of inter-species RDF
        for isid in range(1,nspcs+1):
            for jsid in range(1,nspcs+1):
                agrt= gsmear(rd,agr[isid,jsid],sigma)
                agr[isid,jsid,:] = agrt[:]

    outfile= open(ofname,'w')
    outfile.write('# {0:10s} {1:15s}'.format('rd[i],','agr[0,0,i],'))
    for isid in range(1,nspcs+1):
        for jsid in range(1,nspcs+1):
            outfile.write(' {0:10s}'.format('agr[{0:d}-{1:d}]'.format(isid,jsid)))
    outfile.write('\n')
    for i in range(nr):
        outfile.write(' {0:10.4f} {1:15.7f}'.format(rd[i],agr[0,0,i]))
        for isid in range(1,nspcs+1):
            for jsid in range(1,nspcs+1):
                outfile.write(' {0:10.3f}'.format(agr[isid,jsid,i]))
        outfile.write('\n')
    outfile.close()
    print ' Check '+ofname+' with gnuplot, like'
    print ''
    print " > plot "+ofname+"  us 1:2  w l"
    print ''
