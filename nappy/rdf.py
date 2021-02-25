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
  --out4fp    Flag to write out in general fp.py format. [default: Fault]
  --pairs PAIRS
              Pairs to be extracted, available only if out4fp is specified.
              hyphen-connected, comma separated, e.g.) Li-O,P-O [default: None]
  --skip=NSKIP
              Skip first NSKIP steps from the statistics. [default: 0]
  --no-pairwise
              Not to take averaging by pairwise.
  --plot      Plot figures. [default: False]
  --SQ        Calc and output S(Q) converted from RDF to out.sq
  -q QMAX     Cutoff wavenumber. [default: 25.0]
"""
from __future__ import print_function

import os,sys
import numpy as np
from docopt import docopt

#from nappy.napsys import NAPSystem
import nappy
from nappy.gaussian_smear import gsmear
from nappy.common import get_key

__author__ = "Ryo KOBAYASHI"
__version__ = "200505"

def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def rdf_of_atom(ia,nsys,rmax=5.0,dr=0.1,sigma=0):
    """
    Compute RDF of specified atom.
    """
    #...Radial points
    nr = int(rmax/dr) +1
    rd = np.array([dr*ir+dr/2 for ir in range(nr)],)
    
    nspcs = len(nsys.specorder)
    ndri = np.zeros((nspcs+1,nr),dtype=float)
    spos = nsys.atoms.pos
    sids = nsys.atoms.sid
    natm = nsys.num_atoms()
    pi = spos[ia]
    hmat = nsys.get_hmat()
    #...Compute ndr of atom-ia
    for ja in nsys.neighbors_of(ia,rcut=rmax):
        pj = spos[ja]
        jsid = sids[ja]
        pij = pj -pi
        pij = pij -np.round(pij)
        vij = np.dot(hmat,pij)
        rij2 = np.dot(vij,vij)
        rij = np.sqrt(rij2)
        ir = int(rij/dr)
        ndri[0,ir] += 1.0
        ndri[jsid,ir] += 1.0

    #...Normalize to get raw RDF(ia)
    #.....Total RDF(ia)
    tmp = 4.0 *np.pi *(natm-1) *dr
    for ir in range(1,nr):
        r = dr *ir
        ndri[0,ir] /= tmp*r*r
    #.....Species-decomposed RDF(ia)
    natms = [ float(natm) ]
    for isp in range(1,nspcs+1):
        natms.append(float(nsys.num_atoms(isp)))
    vol = nsys.get_volume()
    isid = sids[ia]
    tmp0 = 4.0 *np.pi *dr /vol
    for jsid in range(1,nspcs+1):
        nj = natms[jsid]
        if jsid == isid:
            tmp = tmp0 *(nj-1)
        else:
            tmp = tmp0 *nj
        for ir in range(1,nr):
            r = dr *ir
            ndri[jsid,ir] /= tmp*r*r

    rdfi = np.zeros(ndri.shape)
    #...Gaussian smearing
    if not sigma == 0:
        #...Total 
        rdfi[0,:] = gsmear(rd,ndri[0,:],sigma)
        #...Species-decomposed
        for jsid in range(1,nspcs+1):
            rdfi[jsid,:] = gsmear(rd,ndri[jsid,:],sigma)

    return rd,rdfi

def rdf_desc_of(ia,nsys,rmax=5.0,dr=0.1):
    """
    RDF descriptor of a given atom-ia.
    RDF descriptor has the information of 1st and 2nd peak positions
    and height of each species.
    For example, in case of Mg-Si-O 3-component system, an atoms has the following data:
       R1(Mg), H1(Mg), R2(Mg), H2(Mg), R1(Si), H1(Si), R2(Si), H2(Si), R1(O), H1(O), R2(O), H2(O)
    Thus the dimension of RDF descriptor of N-component system is 4N and 
    if the peaks are not found, the descriptor values are set 0.0.
    """
    
    rd,rdfi = rdf_of_atom(ia,nsys,rmax=rmax,dr=dr,sigma=2)
    nspcs = len(nsys.specorder)
    nr = len(rd)

    rdf_desci = np.zeros(4*nspcs,dtype=float)
    for jsp in range(1,nspcs+1):
        signs = np.zeros(len(rd),dtype=int)
        rdfij = rdfi[jsp,:]
        for ir in range(1,nr-1):
            diff = (rdfij[ir+1]-rdfij[ir-1])
            if diff < 0.0:
                signs[ir] = -1
            else:
                signs[ir] = +1
        found_1st = False
        for ir in range(1,nr-1):
            if signs[ir] *signs[ir-1] < 0:
                if not found_1st:
                    rdf_desci[(jsp-1)*4 +0] = rd[ir]
                    rdf_desci[(jsp-1)*4 +1] = rdfij[ir]
                    found_1st = True
                elif signs[ir-1] > 0:  # 2nd peak
                    rdf_desci[(jsp-1)*4 +2] = rd[ir]
                    rdf_desci[(jsp-1)*4 +3] = rdfij[ir]
                    break
    return rdf_desci

def rdf_desc(nsys,rmax=5.0,dr=0.1,progress=False):
    """
    Compute RDF descriptor of the given system.
    """
    import sys
    natm = nsys.num_atoms()
    nspcs = len(nsys.specorder)
    desc = np.zeros((natm,4*nspcs),dtype=float)
    interval = max(int(natm/20),1)
    for ia in range(natm):
        if ia % interval == 0 and progress:
            print('ia/natm = {0:8d}/{1:8d}'.format(ia,natm))
            sys.stdout.flush()
        desc[ia,:] = rdf_desc_of(ia,nsys,rmax=rmax,dr=dr)

    return desc

def read_rdf(fname='out.rdf'):
    """
    Read RDF data from a file.
    The format is the same as that of write_normal function.
    """
    with open(fname,'r') as f:
        lines = f.readlines()
    #...Count num of data
    nr = 0
    for line in lines:
        if line[0] != '#':
            nr += 1
    nd = len(line[-1].split()) -1
    rdfs = np.zeros((nr,nd),dtype=float)
    rs = np.zeros(nr,dtype=float)
    ir = 0
    for il,line in enumerate(lines):
        if line[0] == '#': continue
        dat = line.split()
        rs[ir] = float(data[0])
        rdfs[ir,:] = [ float(x) for x in dat[1:]]
    return rs,rdfs
    
def rdf(nsys0,nspcs,dr,rmax,pairwise=False):
    import copy

    natm0= nsys0.num_atoms()
    vol= nsys0.get_volume()
    natms = [ float(natm0) ]
    for ispcs in range(1,nspcs+1):
        natms.append(float(nsys0.num_atoms(ispcs)))

    nsys = copy.deepcopy(nsys0)
    n1,n2,n3= nsys.get_expansion_num(2.0*rmax)
    if not (n1==1 and n2==1 and n3==1):
        print(' Extend system by {0:d}x{1:d}x{2:d}'.format(n1,n2,n3))
        nsys.repeat(n1,n2,n3)

    r2max = rmax*rmax
    nr= int(rmax/dr)+1
    rd= [ dr*ir+dr/2 for ir in range(nr) ]
    hmat = nsys.get_hmat()
    # Since an access to pandas DataFrame is much slower than that to numpy array,
    # use numpy arrays in the most time consuming part.
    poss = np.array(nsys.atoms.pos,)
    sids = np.array(nsys.atoms.sid,)
    natm = len(nsys.atoms)
    nadr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    ndr = np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsys.make_pair_list(rcut=rmax,distance=True)
    for ia in range(natm0):
        isid = sids[ia]
        pi = poss[ia]
        ndr[:,:] = 0.0
        for ja,dij in nsys.neighbors_of(ia,distance=True):
            jsid = sids[ja]
            rrdr= dij/dr
            ir = int(rrdr)
            ndr[0,0,ir] += 1.0
            ndr[isid,jsid,ir] += 1.0
        for ir in range(nr):
            nadr[:,:,ir] += ndr[:,:,ir]

    #...normalize
    if pairwise:
        tmp = 4.0 *np.pi *natms[0]*(natms[0]-1)/vol *dr
        for ir in range(1,nr):
            r= dr *(ir-0.5)
            nadr[0,0,ir] /= tmp*r*r
        for isid in range(1,nspcs+1):
            ni = natms[isid]
            for jsid in range(isid,nspcs+1):
                nj = natms[jsid]
                tmp = 4.0*np.pi*dr  /vol
                if isid == jsid:
                    tmp *= ni*(ni-1)
                else:
                    tmp *= ni*nj
                for ir in range(1,nr):
                    r= dr *(ir-0.5)
                    nadr[isid,jsid,ir] /= tmp*r*r
    else:
        tmp = 4.0 *np.pi *natms[0]*natms[0]/vol *dr
        for ir in range(1,nr):
            r= dr *(ir -0.5)
            nadr[:,:,ir] /= tmp*r*r

    return rd,nadr

def gr_to_SQ(rs,gr,rho,rcut=5.0,qcut=25.0):
    """
    Convert RDF to S(Q).
    """
    nbins = len(rs)
    sq = np.zeros(nbins)
    qs = np.zeros(nbins)
    dq = qcut /nbins
    dr = rcut /nbins
    for ib in range(nbins):
        q = (ib+0.5)*dq
        tmp = 0.0
        qs[ib] = q
        for jb in range(1,nbins):
            r = (jb +0.5)*dr
            jbm = jb -1
            rm = (jbm +0.5)*dr
            tmp1 = (gr[jbm]-1.0)*np.sin(q*rm) /(q*rm) *rm*rm
            tmp2 = (gr[jb] -1.0)*np.sin(q*r) /(q*r) *r*r
            tmp = tmp +0.5*dr *(tmp1+tmp2)
        sq[ib] = 1.0 +4.0*np.pi*rho*tmp
    return qs,sq
    
def rdf_average(infiles,nr,specorder,dr=0.1,rmax=3.0,pairwise=False):
    nspcs = len(specorder)
    agr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print("[Error] File, {0}, does not exist !!!".format(infname))
            sys.exit()
        nsys = nappy.io.read(fname=infname,specorder=specorder)
        print(' File =',infname)
        rd,gr= rdf(nsys,nspcs,dr,rmax,pairwise=pairwise)
        agr += gr
        nsum += 1
    agr /= nsum
    return rd,agr

def write_normal(fname,specorder,nspcs,rd,agr,nr):
    """
    Write out RDF data in normal RDF format.
    """
    outfile= open(fname,'w')
    outfile.write('# 1:{0:10s}  2:all-all,  '.format('rd[i],'))
    n = 2
    for isid in range(1,nspcs+1):
        si = specorder[isid-1]
        for jsid in range(isid,nspcs+1):
        # for jsid in range(1,nspcs+1):
            sj = specorder[jsid-1]
            n += 1
            outfile.write('  {0:d}:{1:s}-{2:s},   '.format(n,si,sj))
    outfile.write('\n')
    for i in range(nr-1):
        outfile.write(' {0:10.4f} {1:13.5e}'.format(rd[i],agr[0,0,i]))
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
            # for jsid in range(1,nspcs+1):
                outfile.write(' {0:12.4e}'.format(agr[isid,jsid,i]))
        outfile.write('\n')
    outfile.close()
    return None

def write_out4fp(fname,specorder,nspcs,agr,nr,rmax,pairs,nperline=6):
    """
    Write out RDF data in general fp.py format.

    Parameters
    ----------
    nperline : int
           Number of data in a line. [default: 6]
    """
    ndat = (nr-1) *len(pairs)
    data = np.zeros(ndat)
    n = 0
    for pair in pairs:
        isid,jsid = pair
        for i in range(nr-1):
            data[n] = agr[isid,jsid,i]
            n += 1

    with open(fname,'w') as f:
        f.write('# RDF for pairs: ')
        for pair in pairs:
            si = specorder[pair[0]-1]
            sj = specorder[pair[1]-1]
            f.write(' {0:s}-{1:s},'.format(si,sj))
        f.write('\n')
        f.write('# rmax, nr = {0:.3f}, {1:d}\n'.format(rmax,nr))
        f.write('#\n')
        #...Num of data, weight for the data
        f.write(' {0:6d}  {1:7.3f}\n'.format(ndat, 1.0))
        j0 = 0
        while True:
            f.write('  '.join('{0:12.4e}'.format(data[j]) for j in range(j0,j0+nperline) if j < ndat))
            f.write('\n')
            j0 += nperline
            if j0 >= ndat:
                break

    return None

def write_sq_normal(fname,qs,sq):
    """
    Write S(Q) data in normal, gnuplot-readable format.
    """
    nd = len(qs)
    with open(fname,'w') as f:
        f.write('# S(Q) computed in rdf.py\n')
        f.write('# Q,             S(Q)\n')
        for i in range(nd):
            f.write(' {0:10.4f}  {1:10.5f}\n'.format(qs[i],sq[i]))
    return None
        

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
    out4fp = args['--out4fp']
    if out4fp:
        pairwise = True
        pairs0 = args['--pairs'].split(',')
        pairs = []
        for pair in pairs0:
            spi,spj = pair.split('-')
            isid = specorder.index(spi)+1
            jsid = specorder.index(spj)+1
            if jsid < isid:
                itmp = jsid
                jsid = isid
                isid = itmp
            pairs.append((isid,jsid))
    else:
        no_pairwise = args['--no-pairwise']
        pairwise = not no_pairwise
    plot = args['--plot']
    nskip = int(args['--skip'])
    sq = args['--SQ']
    if sq:
        qmax = float(args['-q'])

    nspcs = len(specorder)
    if nspcs < 1:
        raise ValueError('--specorder must be set.')

    if len(infiles) > 1:
        infiles.sort(key=get_key,reverse=True)
    del infiles[:nskip]

    nr= int(rmax/dr) +1
    rd,agr= rdf_average(infiles,nr,specorder,dr=dr,rmax=rmax,
                        pairwise=pairwise)

    if not sigma == 0:
        print(' Gaussian smearing...')
        #...Smearing of total RDF
        agrt= gsmear(rd,agr[0,0],sigma)
        agr[0,0,:] = agrt[:]
        #...Smearing of inter-species RDF
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
                agrt= gsmear(rd,agr[isid,jsid],sigma)
                agr[isid,jsid,:] = agrt[:]

    if out4fp:
        write_out4fp(ofname,specorder,nspcs,agr,nr,rmax,pairs)
    else:
        write_normal(ofname,specorder,nspcs,rd,agr,nr,)

    if sq:
        nsys = nappy.io.read(infiles[0])
        rho = float(nsys.num_atoms()) /nsys.get_volume()
        qs,sq = gr_to_SQ(rd,agr[0,0,:],rho,rcut=rmax,qcut=qmax)
        write_sq_normal('out.sq',qs,sq)

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
