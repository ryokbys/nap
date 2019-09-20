#!/usr/bin/env python
"""
Calculate the radial distribution function (RDF) from files.
Statistical averaging about atoms in a file and over files are taken.

Usage:
  rdf.py [options] INFILE [INFILE...]

Options:
  -h, --help  Show this help message and exit.
  -d DR       Width of the bin. [default: 0.1]
  -r RMAX     Cutoff radius of radial distribution. [default: 5.0]
  --gsmear=SIGMA
              Width of Gaussian smearing, zero means no smearing. [default: 0]
  -o OUT      Output file name. [default: out.rdf]
  --specorder=SPECORDER
              Order of species separated by comma, like, --specorder=W,H. [default: None]
  --skip=NSKIP 
              Skip first NSKIP steps from the statistics. [default: 0]
  --no-average
              Not to take average over files.
  --no-normalize
              Not to normalize by the density.
  --plot      Plot figures. [default: False]
"""
from __future__ import print_function

import os,sys
import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem
from nappy.gaussian_smear import gsmear
from nappy.common import get_key

def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def compute_ndr(ia,isid,dr,rmax,nsys,nspcs):
    """
    Compute number of atoms in the every shell [r:r+dr] up to *rmax*.
    This routine can be only applied to cubic systems.
    """
    nr= int(rmax/dr) +1
    #ndr= np.zeros(nr,dtype=np.int)
    ndr = np.zeros((nspcs+1,nspcs+1,nr),dtype=np.int)
    natm= nsys.natm
    hmat= (nsys.alc *np.array([nsys.a1,nsys.a2,nsys.a3])).transpose()
    pi= nsys.poss[ia]
    for ja in range(natm):
        if ja == ia:
            continue
        jsid = nsys.sids[ja]
        pj= nsys.poss[ja]
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

def rdf(nsys,nspcs,dr,rmax,normalize=True):

    natm0= nsys.natm
    vol= nsys.volume()
    rho= float(natm0)/vol

    n1,n2,n3= nsys.get_expansion_num(2.0*rmax)
    if not (n1==1 and n2==1 and n3==1):
        print(' system to be repeated, n1,n2,n3=',n1,n2,n3)
        nsys.repeat(n1,n2,n3)

    nr= int(rmax/dr)+1
    nadr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    rd= [ dr*ir+dr/2 for ir in range(nr) ]
    for ia in range(natm0):
        isid = nsys.sids[ia]
        ndr= compute_ndr(ia,isid,dr,rmax,nsys,nspcs)
        for ir in range(nr):
            nadr[:,:,ir] += ndr[:,:,ir]
    #nadr /= nsrc

    #...normalize
    if normalize:
        for ir in range(1,nr):
            r= dr *ir
            nadr[:,:,ir]= nadr[:,:,ir]/(4.0*np.pi*rho*r*r*dr)

    return rd,nadr,natm0

def rdf_average(infiles,nr,specorder,dr=0.1,rmax=3.0,average=True,
                normalize=True):
    nspcs = len(specorder)
    agr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print("[Error] File, {0}, does not exist !!!".format(infname))
            sys.exit()
        nsys= NAPSystem(fname=infname,specorder=specorder)
        print(' File =',infname)
        rd,gr,n= rdf(nsys,nspcs,dr,rmax,normalize=normalize)
        nsum += n
        agr += gr
    # agr /= len(infiles)
    if average:
        agr /= nsum
    return rd,agr

def plot_figures(specorder,rd,agr):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(context='talk',style='ticks')

    nspcs = len(specorder)
    plt.figure(figsize=(8,6))
    x = rd
    y = agr[0,0,:]
    plt.plot(x,y,'r-',label='Total RDF')
    plt.xlabel('Distance (A)')
    plt.ylabel('RDF')
    plt.legend()
    plt.savefig("graph_rdf_total.png", format='png', dpi=300, bbox_inches='tight')

    if nspcs == 1:
        return
    plt.clf()
    fig, axes = plt.subplots(nspcs,nspcs,figsize=(15,10),sharex=True)
    for i in range(nspcs):
        isp = i + 1
        for j in range(nspcs):
            jsp = j + 1
            if j < i:
                axes[i,j].axis('off')
                continue
            ax = axes[i,j]
            y = agr[isp,jsp,:]
            ax.plot(x,y,'r-')
            ax.set_title('{0:d}-{1:d}'.format(isp,jsp))
            if isp==jsp:
                ax.set_xlabel('Distance (A)')
            if isp==1 and jsp==1:
                ax.set_ylabel('RDF')
    plt.savefig("graph_rdfs.png", format='png', dpi=300, bbox_inches='tight')
    return

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    dr= float(args['-d'])
    rmax= float(args['-r'])
    sigma= int(args['--gsmear'])
    ofname= args['-o']
    specorder = [ x for x in args['--specorder'].split(',') ]
    if specorder == ['None']:
        specorder = []
    no_average = args['--no-average']
    no_normalize = args['--no-normalize']
    average = not no_average
    normalize = not no_normalize
    plot = args['--plot']
    nskip = int(args['--skip'])

    nspcs = len(specorder)
    if nspcs < 1:
        raise ValueError('--specorder must be set.')

    infiles.sort(key=get_key,reverse=True)
    del infiles[:nskip]

    nr= int(rmax/dr) +1
    rd,agr= rdf_average(infiles,nr,specorder,dr=dr,rmax=rmax,
                        average=average,normalize=normalize)

    if not sigma == 0:
        print(' Gaussian smearing...')
        #...Smearing of total RDF
        agrt= gsmear(rd,agr[0,0],sigma)
        agr[0,0,:] = agrt[:]
        #...Smearing of inter-species RDF
        for isid in range(1,nspcs+1):
            for jsid in range(1,nspcs+1):
                agrt= gsmear(rd,agr[isid,jsid],sigma)
                agr[isid,jsid,:] = agrt[:]

    outfile= open(ofname,'w')
    outfile.write('# 1:{0:10s}  2:all-all,  '.format('rd[i],'))
    n = 2
    for isid in range(1,nspcs+1):
        si = specorder[isid-1]
        for jsid in range(isid,nspcs+1):
            sj = specorder[jsid-1]
            n += 1
            #outfile.write(' {0:d}:{1:10s}'.format(n,'agr[{0:d}-{1:d}]'.format(isid,jsid)))
            outfile.write('  {0:d}:{1:s}-{2:s},   '.format(n,si,sj))
    outfile.write('\n')
    for i in range(nr):
        outfile.write(' {0:10.4f} {1:13.5e}'.format(rd[i],agr[0,0,i]))
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
                outfile.write(' {0:12.4e}'.format(agr[isid,jsid,i]))
        outfile.write('\n')
    outfile.close()

    if plot:
        plot_figures(rd,agr)
        print('')
        print(' RDF graphes are plotted.')
        if nspcs == 1:
            print(' Check graph_rdf_total.png')
        else:
            print(' Check graph_rdf_total.png and graph_rdfs.png')
    else:
        print(' Check {0:s} with gnuplot, like'.format(ofname))
        print('')
        print(" > plot '{0:s}' us 1:2  w l".format(ofname))
        print('')
