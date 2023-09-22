#!/usr/bin/env python
"""
Calculate the radial distribution function (RDF) from files.
It takes statistical averaging over atoms in a file and over files.

Usage:
  {0:s} [options] INFILE [INFILE...]

Options:
  -h, --help  Show this help message and exit.
  -d DR       Width of the bin. [default: 0.1]
  -r,--rmax RMAX
              Cutoff radius of radial distribution. [default: 5.0]
  --rmin RMIN 
              Minimum radius to be considered. [default: 0.0]
  --gsmear=SIGMA
              Width of Gaussian smearing, zero means no smearing. [default: 0]
  --nnmax=NNMAX
              Max num of neighbors when counting neighbors. [default: 100]
  -o OUT      Name of file to be written in addition to out.rdf if specified. [default: None]
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
  --qmax QMAX    Cutoff wavenumber. [default: 20.0]
  --qmin QMIN    Shortest wavenumber. [default: 0.7]
  --scatter-length LENGTHS
              Scattering lengths of corresponding species. [default: None]
  --fortran   Try using fortran function for computing RDF.
"""

import os,sys
from datetime import datetime
import numpy as np
import copy
from docopt import docopt

import nappy
from nappy.gaussian_smear import gsmear
from nappy.common import get_key

__author__ = "Ryo KOBAYASHI"
__version__ = "230107"

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
    nr = int(rmax/dr) #+1
    rd = np.array([dr*ir+dr/2 for ir in range(nr)],)
    
    nspcs = len(nsys.specorder)
    ndri = np.zeros((nspcs+1,nr),dtype=float)
    spos = nsys.get_scaled_positions()
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
        for ir in range(1,nr):
            diff = (rdfij[ir+1]-rdfij[ir-1])
            if diff < 0.0:
                signs[ir] = -1
            else:
                signs[ir] = +1
        found_1st = False
        for ir in range(1,nr):
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
    
def rdf(nsys0,nspcs,dr,nr,rmax0,pairwise=False,rmin=0.0,
        nnmax=100,fortran=False,mask=None):
    """
    Compute RDF of the givin system.
    Number of bins is determined from DR, RMAX0, and RMIN as, int((RMAX0-RMIN)/DR)+1.
    Position of a bin is defined as the point of center of the bin, r_i = DR*(i+0.5).
    """
    natm0= nsys0.num_atoms()
    #...NOTE: mask is not available now
    # if mask:
    #     if type(mask) != list:
    #         raise TypeError('mask is not list.')
    #     if len(mask) != natm0:
    #         raise ValueError('len(mask) != len(nsys0)')
        
    rmax = rmin +dr*nr  # Use corrected rmax to cover regions of NR bins
    r2max = rmax*rmax
    
    nsys = copy.deepcopy(nsys0)
    n1,n2,n3= nsys.get_expansion_num(3.0*rmax)
    if not (n1==1 and n2==1 and n3==1):
        print(' Extend system by {0:d}x{1:d}x{2:d}'.format(n1,n2,n3))
        nsys.repeat(n1,n2,n3)

    natm = len(nsys)
    natms = [ float(natm) ]
    for ispcs in range(1,nspcs+1):
        natms.append(float(nsys.num_atoms(ispcs)))
    vol= nsys.get_volume()
        
    hmat = nsys.get_hmat()
    # Since an access to pandas DataFrame is much slower than that to numpy array,
    # use numpy arrays in the most time consuming part.
    poss = nsys.get_scaled_positions()
    sids = np.array(nsys.atoms.sid,)
    natm = len(nsys)
    if fortran:
        try:
            # import nappy.pmd.mods as pmods
            import nappy.pmd.pmd_wrapper as pw
            hmati = nsys.get_hmat_inv()
            tags = nsys.get_tags()
            iprint = 0
            l1st = True
            # lspr = pmods.pairlist.mk_lspr_sngl(natm,nnmax,tags,
            #                                    poss.T,rmax,hmat,hmati,
            #                                    iprint,l1st)
            # rd,rdfs= pmods.distfunc.calc_rdf(tags,hmat,poss.T,rmax,rmin,
            #                                  lspr,iprint,l1st,pairwise,
            #                                  nspcs,nr)
            rd,rdfs= pw.wrap_calc_rdf(poss.T,tags,hmat,hmati,rmax,rmin,l1st,
                                      pairwise,nr,nspcs)
            return rd, rdfs.T
        except Exception as e:
            print(' Since failed to use the fortran routines, use python instead...')
            print(e)
            pass
    
    rd= np.array([ rmin +dr*(ir+0.5) for ir in range(nr) ])
    nadr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    ndr = np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsys.make_pair_list(rcut=rmax,nnmax=nnmax)
    hmat = nsys.get_hmat()
    for ia in range(natm):
        isid = sids[ia]
        # ndr[:,:] = 0.0
        pi = poss[ia]
        for ja in nsys.neighbors_of(ia,):
            jsid = sids[ja]
            pj = poss[ja]
            pij = pj -pi -np.round(pj-pi)
            vij = np.dot(hmat,pij)
            dij2 = np.dot(vij,vij)
            dij = np.sqrt(dij2)
            rrdr= (dij-rmin)/dr
            if rrdr < 0.0:
                continue
            ir = min(int(rrdr),nr-1)
            if ir < 0:
                print('Something is wrong: ir<0 ')
                print(ia,ja,dij,rrdr,ir)
                raise
            # ndr[0,0,ir] += 1.0
            # ndr[isid,jsid,ir] += 1.0
            nadr[0,0,ir] += 1.0
            nadr[isid,jsid,ir] += 1.0
        # for ir in range(nr):
        #     nadr[:,:,ir] += ndr[:,:,ir]

    #...normalize
    if pairwise:
        #...Total
        tmp = 4.0 *np.pi *natms[0]*(natms[0]-1)/vol *dr
        for ir in range(nr):
            # r= rmin +dr*(ir+0.5)
            r = rd[ir]
            nadr[0,0,ir] /= tmp*r*r
        #...Pairwise
        for isid in range(1,nspcs+1):
            ni = natms[isid]
            if ni == 0: continue
            for jsid in range(isid,nspcs+1):
                nj = natms[jsid]
                if nj == 0: continue
                tmp = 4.0*np.pi*dr  /vol
                if isid == jsid:
                    if ni == 1: continue
                    tmp *= ni*(ni-1)
                else:
                    tmp *= ni*nj
                for ir in range(nr):
                    # r= dr *(ir-0.5)
                    r = rd[ir]
                    nadr[isid,jsid,ir] /= tmp*r*r
    else:
        tmp = 4.0 *np.pi *natms[0]*(natms[0]-1)/vol *dr
        for ir in range(nr):
            # r= dr *(ir -0.5)
            r = rd[ir]
            nadr[:,:,ir] /= tmp*r*r

    return rd,nadr

def rdf_average(infiles,specorder,dr=0.1,rmin=0.0,rmax=3.0,pairwise=False,nnmax=100,fortran=False):
    nspcs = len(specorder)
    tiny = 1.0e-8
    nr = int((rmax-rmin+tiny)/dr) #+1 , no need to add 1
    agr= np.zeros((nspcs+1,nspcs+1,nr),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print("[Error] File, {0}, does not exist !!!".format(infname))
            sys.exit()
        nsys = nappy.io.read(fname=infname,specorder=specorder)
        print(' File =',infname)
        rd,gr= rdf(nsys,nspcs,dr,nr,rmax,rmin=rmin,
                   pairwise=pairwise,nnmax=nnmax,fortran=fortran)
        if rd.shape[-1] != nr:
            raise ValueError('The shape of radius data is wrong.')
        agr += gr
        nsum += 1
    agr /= nsum
    return rd,agr

def gr_to_SQ(rs,gr,rho,qmin=0.7,qmax=20.0,nq=100):
    """
    Convert g(r) to S(Q).
    """
    nr = len(rs)
    rmin = min(rs)
    rmax = max(rs)
    sq = np.zeros(nq)
    qs = np.zeros(nq)
    dq = (qmax-qmin) /nq
    dr = (rmax-rmin) /nr
    for iq in range(nq):
        q = iq*dq +qmin
        tmp = 0.0
        qs[iq] = q
        if q < 1.0e-15:
            raise ValueError('qmin should not be 0.')
        for jr in range(1,nr):
            #r = jr*dr +rmin
            jrm = jr -1
            r = rs[jr]
            #rm = jrm*dr +rmin
            rm = rs[jrm]
            if abs(rm) < 1.0e-15:
                continue
            tmp1 = (gr[jrm]-1.0)*np.sin(q*rm) /(q*rm) *rm*rm
            tmp2 = (gr[jr] -1.0)*np.sin(q*r) /(q*r) *r*r
            tmp += 0.5 *(tmp1+tmp2) *dr
        sq[iq] = 1.0 +4.0*np.pi*rho*tmp
    return qs,sq

def SQ_to_gr(qs,sq,rho,rmin=0.0,rmax=5.0,nr=100):
    """
    Convert S(Q) to g(r).
    """
    nbins = len(qs)
    qmin = min(qs)
    qmax = max(qs)
    rs = np.zeros(nr)
    gr = np.zeros(nr)
    dq = (qmax-qmin) /nbins
    dr = (rmax-rmin) /nr
    for ir in range(nr):
        r = ir*dr +rmin
        tmp = 0.0
        rs[ir] = r
        if r < 1.0e-15:
            gr[ir] = 0.0
            continue
        for jb in range(1,nbins):
            #q = jb*dq +qmin
            jbm = jb -1
            #qm = jbm*dq +qmin
            q = qs[jb]
            qm = qs[jbm]
            if qm < 1.0e-15:
                continue
            tmp1 = (sq[jbm]-1.0) *np.sin(qm*r) *qm
            tmp2 = (sq[jb] -1.0) *np.sin(q*r) *q
            tmp += 0.5*(tmp1+tmp2)*dq
        gr[ir] = 1.0 +1.0/(2.0*np.pi**2 *r *rho) *tmp
    return rs, gr

def write_rdf_normal(fname,specorder,nspcs,rd,agr,nr):
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
    for i in range(nr):
        outfile.write(' {0:10.4f} {1:13.5e}'.format(rd[i],agr[0,0,i]))
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
            # for jsid in range(1,nspcs+1):
                outfile.write(' {0:12.4e}'.format(agr[isid,jsid,i]))
        outfile.write('\n')
    outfile.close()
    return None

def write_rdf_out4fp(fname,specorder,nspcs,agr,nr,rmax,pairs=None,rmin=0.0,nperline=6):
    """
    Write out RDF data in general fp.py format.

    Parameters
    ----------
    nperline : int
           Number of data in a line. [default: 6]
    """
    if pairs != None:
        if type(pairs) not in (list,tuple):
            raise TypeError('pairs must be list or tuple.')
        ndat = nr *len(pairs)
        data = np.zeros(ndat)
        n = 0
        for pair in pairs:
            isid,jsid = pair
            for i in range(nr):
                data[n] = agr[isid,jsid,i]
                n += 1
    else:
        print('Since pairs are not specified, use total RDF instead.')
        ndat = nr
        data = np.zeros(ndat)
        n = 0
        for i in range(nr):
            data[n] = agr[0,0,i]
            n += 1
    
    with open(fname,'w') as f:
        cmd = ' '.join(s for s in sys.argv)
        f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        f.write('#  {0:s}\n'.format(cmd))
        if pairs != None:
            f.write('# RDF for pairs: ')
            for pair in pairs:
                isid,jsid = pair
                if isid == 0:
                    si = 'All'
                else:
                    si = specorder[pair[0]-1]
                if jsid == 0:
                    sj = 'All'
                else:
                    sj = specorder[pair[1]-1]
                f.write(' {0:s}-{1:s},'.format(si,sj))
            f.write('\n')
        f.write('# rmin,rmax, nr = {0:.3f}, {1:.3f}, {2:d}\n'.format(rmin,rmax,nr))
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
        

def write_sq_out4fp(fname,qs,sq,nperline=6):
    """
    Write S(Q) data in normal, gnuplot-readable format.
    """
    nd = len(qs)
    with open(fname,'w') as f:
        f.write('# S(Q) computed in rdf.py\n')
        f.write('# Output S(Q) in out4fp format.\n')
        f.write('#\n')
        #...Num of data, weight for the data
        f.write(' {0:6d}  {1:7.3f}\n'.format(nd, 1.0))
        j0 = 0
        while True:
            f.write('  '.join('{0:12.4e}'.format(sq[j]) for j in range(j0,j0+nperline) if j < nd))
            f.write('\n')
            j0 += nperline
            if j0 >= nd:
                break
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
            # ax.set_title('{0:d}-{1:d}'.format(isp,jsp))
            ax.text(0.05, 0.8, '{0:s}-{1:s}'.format(specorder[i],specorder[j]),
                    transform=ax.transAxes, ha="left")
            if isp==jsp:
                ax.set_xlabel('Distance (A)')
            if isp==1 and jsp==1:
                ax.set_ylabel('RDF')
    plt.savefig("graph_rdfs.png", format='png', dpi=300, bbox_inches='tight')
    return

def nbplot(nsys,dr=0.1,rmin=0.0,rmax=5.0,nnmax=200,pairs=None,sigma=0):
    """
    Plot RDFs of given nsys on the jupyter notebook.

    pairs option should be in the form of list of species pair list, e.g., (('Si','Si'),('Si','O')).
    """
    if not 'JPY_PARENT_PID' in os.environ:
        raise Exception('This routine must be called on jupyter notebook.')

    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
        sns.set(context='talk',style='ticks')
    except:
        pass

    nspcs = len(nsys.specorder)
    try:
        nr = int((rmax-rmin)/dr) #+1
        rd,gr= rdf(nsys,nspcs,dr,nr,rmax,rmin=rmin,nnmax=nnmax)
    except:
        raise Exception('rdf(..) failed.')

    if sigma > 0:
        #...Smearing of total RDF
        grt= gsmear(rd,gr[0,0,:],sigma)
        gr[0,0,:] = grt[:]
        #...Smearing of inter-species RDF
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
                grt= gsmear(rd,gr[isid,jsid,:],sigma)
                gr[isid,jsid,:] = grt[:]

    plt.figure(figsize=(8,6))
    x = rd
    if not pairs:
        for i in range(nspcs):
            isp = i +1
            spi = nsys.specorder[i]
            for j in range(nspcs):
                jsp =  j +1
                spj = nsys.specorder[j]
                if j<i:
                    continue
                y = gr[isp,jsp,:]
                plt.plot(x,y,label='{0:s}-{1:s}'.format(spi,spj))
    else:  # only specified pairs are plotted in addition to total RDF
        for p in pairs:
            spi,spj = p
            try:
                isp = nsys.specorder.index(spi) +1
                jsp = nsys.specorder.index(spj) +1
            except:
                raise ValueError('No such species or pairs.')
            y = gr[isp,jsp,:]
            plt.plot(x,y,label='{0:s}-{1:s}'.format(spi,spj))

    #...Total RDF
    y = gr[0,0,:]
    plt.plot(x,y,'r--',label='Total RDF')

    plt.xlabel('Distance ($\mathrm{\AA}$)')
    plt.ylabel('RDF')
    plt.legend(bbox_to_anchor=(1.05,1))
    plt.show()
    return None

def main():

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])), version=__version__)
    
    infiles= args['INFILE']
    dr= float(args['-d'])
    rmax= float(args['--rmax'])
    rmin = float(args['--rmin'])
    sigma= int(args['--gsmear'])
    nnmax = int(args['--nnmax'])
    ofname= args['-o']

    if nnmax < int(rmax**3):
        newnnmax = int(rmax**3)
        print(' nnmax is updated from {0:d} to {1:d} according to rmax.'.format(nnmax,newnnmax))
        nnmax = newnnmax

    if ofname == 'None':
        ofname = None
    specorder = [ x for x in args['--specorder'].split(',') ]
    if specorder == ['None']:
        specorder = []
    plot = args['--plot']
    nskip = int(args['--skip'])
    SQ = args['--SQ']
    if SQ:
        qmax = float(args['--qmax'])
        qmin = float(args['--qmin'])
        lscatter = [ float(x) for x in args['--scatter-length'].split(',') ]
        if len(lscatter) != len(specorder):
            raise ValueError('--scatter-length is not set correctly.')
    out4fp = args['--out4fp']
    fortran = args['--fortran']
    if out4fp and ofname is None:
        raise ValueError("Output file name must be specified with option -o.")
    if out4fp and not SQ:
        pairwise = True
        pairs0 = args['--pairs'].split(',')
        pairs = []
        for pair in pairs0:
            spi,spj = pair.split('-')
            try:
                isid = specorder.index(spi)+1
            except:
                isid = 0
            try:
                jsid = specorder.index(spj)+1
            except:
                jsid = 0
            if jsid < isid:
                itmp = jsid
                jsid = isid
                isid = itmp
            pairs.append((isid,jsid))
    else:
        no_pairwise = args['--no-pairwise']
        pairwise = not no_pairwise
        pairs = None

    nspcs = len(specorder)
    if nspcs < 1:
        raise ValueError('--specorder must be set.')

    if len(infiles) > 1:
        infiles.sort(key=get_key,reverse=True)
    del infiles[:nskip]
    if len(infiles) < 1:
        raise ValueError('No input files to be processed.')
    print(' Number of files to be processed: ',len(infiles))

    tiny = 1.0e-8
    nr= int((rmax-rmin+tiny)/dr) #+1
    rd,agr= rdf_average(infiles,specorder,dr=dr,rmin=rmin,rmax=rmax,
                        pairwise=pairwise,nnmax=nnmax,fortran=fortran)

    if not sigma == 0:
        #print(' Gaussian smearing...')
        #...Smearing of total RDF
        agrt= gsmear(rd,agr[0,0,:],sigma)
        agr[0,0,:] = agrt[:]
        #...Smearing of inter-species RDF
        for isid in range(1,nspcs+1):
            for jsid in range(isid,nspcs+1):
                agrt= gsmear(rd,agr[isid,jsid,:],sigma)
                agr[isid,jsid,:] = agrt[:]

    if SQ:
        nsys = nappy.io.read(infiles[0])
        rho = float(nsys.num_atoms()) /nsys.get_volume()
        if nspcs > 1:
            #...Redfine total RDF as weighted sum of g_{ij}(r) in case of multiple species
            natms = [ float(nsys.num_atoms()) ]
            cs = [ 1.0 ] 
            for ispcs in range(1,nspcs+1):
                natms.append(float(nsys.num_atoms(ispcs)))
                cs.append(natms[ispcs]/natms[0])
            bmean = 0.0
            for isid in range(1,nspcs+1):
                bi = lscatter[isid-1]
                ci = cs[isid]
                bmean += ci*bi
            agr[0,0,:] = 0.0
            for isid in range(1,nspcs+1):
                bi = lscatter[isid-1]
                ci = cs[isid]
                for jsid in range(isid,nspcs+1):
                    bj = lscatter[jsid-1]
                    cj = cs[jsid]
                    wij = ci*cj*bi*bj/bmean
                    if isid == jsid:
                        agr[0,0,:] += agr[isid,jsid,:] *wij
                    else:
                        agr[0,0,:] += 2.0*agr[isid,jsid,:] *wij
        qs,sqs = gr_to_SQ(rd,agr[0,0,:],rho,qmin=0.7,qmax=qmax,nq=200)

    #...Regardless ofname, write out.rdf in normal format
    write_rdf_normal('out.rdf',specorder,nspcs,rd,agr,nr,)
    if SQ:
        write_sq_out4fp('out.sq',qs,sqs)
    #...Format of output (named by ofname) depends on out4fp
    if ofname is not None:
        if out4fp:
            write_rdf_out4fp(ofname,specorder,nspcs,agr,nr,rmax,pairs=pairs,rmin=rmin)
        else:
            write_rdf_normal(ofname,specorder,nspcs,rd,agr,nr,)

    if plot:
        plot_figures(specorder,rd,agr)
        print('')
        print(' RDF graphes are plotted.')
        if nspcs == 1:
            print(' Check graph_rdf_total.png')
        else:
            print(' Check graph_rdf_total.png and graph_rdfs.png')
    else:
        print(' Check out.rdf with gnuplot, like')
        print("   gnuplot> plot 'out.rdf' us 1:2  w l")
        print('')
        if ofname is not None:
            print(" In addition to out.rdf, {0:s} is also written.".format(ofname))
            print('')

    return None

if __name__ == "__main__":

    main()
