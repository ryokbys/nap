#!/usr/bin/env python
"""
Calculate the angular distribution function (ADF) from files.
Take an average over atoms in a file or files.

Usage:
  {0:s} [options] INFILE [INFILE...]

Options:
  -h, --help  Show this help message and exit.
  -w DEG      Width of the angular degree. [default: 1.0]
  -r RCUT     Cutoff radius of the bonding pair. [default: 3.0]
  --nnmax NNMAX
              Max num of neighbors when counting neighbors. [default: 100]
  --gsmear=SIGMA
              Width of Gaussian smearing, zero means no smearing. [default: 0]
  --specorder=SPECORDER
              Order of species separated by comma, like, --specorder=Si,O. [default: None]
  --triplets=TRIPLETS
              Triplets whose angles are to be computed. Three species should be specified connected by hyphen,
              and separated by comma, in which first one is the center of the three.
              e.g.) P-O-O,Li-O-O. [default: None]
  -o OUT      Output file name [default: None]
  --out4fp    Flag to write out in general fp.py format. [default: Fault]
  --skip=NSKIP
              Skip first NSKIP steps from the statistics. [default: 0]
  --plot      Plot figures. [default: False]
  --fortran   Try using fortran routine for ADF calculation.
  --format FORMAT
              Input file format. [default: None]
"""

import os,sys
import numpy as np
from docopt import docopt
from datetime import datetime

from nappy.gaussian_smear import gsmear
from nappy.common import get_key
from nappy.io import read

__author__ = "Ryo KOBAYASHI"
__version__ = "230107"

def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def adf_atom(ia,dang,rcut,nsys,poss,lspr,symbols,sj,sk):
    """
    Compute number of atoms in the every range of angle [0:180].
    """
    na= int(180.0/dang) +1
    hmat= nsys.get_hmat()
    nda= np.zeros(na,dtype=int)
    natm= nsys.num_atoms()
    rcut2= rcut*rcut
    # pi= nsys.get_atom_attr(ia,'pos')
    # lspri = nsys.get_atom_attr(ia,'lspr')
    pi = poss[ia]
    lspri = lspr[ia]
    for ji in range(len(lspri)):
        ja= lspri[ji]
        if ja == ia:
            continue
        sji = symbols[ja]
        if sji not in (sj,sk):
            continue
        # pj= nsys.get_atom_attr(ja,'pos')
        pj = poss[ja]
        pij= pj-pi
        pij= pij -np.round(pij)
        vij= np.dot(hmat,pij)
        rij2= np.dot(vij,vij)
        if rij2 >= rcut2:
            continue
        rij= np.sqrt(rij2)
        for ki in range(ji+1,len(lspri)):
            ka= lspri[ki]
            # if ka == ia or ka <= ja:
            #     continue
            ski = symbols[ka]
            if set((sji,ski)) != set((sj,sk)):
                continue
            # pk= nsys.get_atom_attr(ka,'pos')
            pk = poss[ka]
            pik= pk-pi
            pik= pik -np.round(pik)
            vik= np.dot(hmat,pik)
            rik2= np.dot(vik,vik)
            if rik2 >= rcut2:
                continue
            rik= np.sqrt(rik2)
            cs= np.dot(vij,vik)/rij/rik
            if cs <= -1.0:
                rad= np.pi
            else:
                rad= np.arccos(cs)
            deg= rad/np.pi *180.0
            nda[int(deg/dang)] += 1
    return nda

def adf(nsys,dang,rcut,triplets,fortran=False,nnmax=100):

    natm0= nsys.num_atoms()

    n1,n2,n3= nsys.get_expansion_num(2.0*rcut)
    if not (n1==1 and n2==1 and n3==1):
        print(' system to be repeated, n1,n2,n3=',n1,n2,n3)
        nsys.repeat(n1,n2,n3)
    nsys.assign_pbc()

    na= int(180.0/dang)
    poss = nsys.get_scaled_positions()

    if fortran:
        try:
            # import nappy.pmd.mods as pmods
            import nappy.pmd.pmd_wrapper as pw
            natm = len(nsys)
            tags = nsys.get_tags()
            hmat = nsys.get_hmat()
            hmati = nsys.get_hmat_inv()
            iprint = 0
            l1st = True
            itriples = np.zeros((len(triplets),3),dtype=int)
            specorder = nsys.specorder
            for it,t in enumerate(triplets):
                si,sj,sk = t
                isp = specorder.index(si)+1
                jsp = specorder.index(sj)+1
                ksp = specorder.index(sk)+1
                itriples[it,0] = isp
                itriples[it,1] = jsp
                itriples[it,2] = ksp
            # lspr = pmods.pairlist.mk_lspr_sngl(natm,nnmax,tags,poss.T,
            #                                    rcut,hmat,hmati,
            #                                    iprint,l1st)
            # angd,adfs = pmods.distfunc.calc_adf(tags,hmat,poss.T,rcut,
            #                                     lspr,itriples.T,
            #                                     dang,na)
            angd,adfs = pw.wrap_calc_adf(poss.T,tags,hmat,hmati,rcut,
                                         itriples.T,na,l1st)
            return angd, adfs.T
        except Exception as e:
            print(' Since failed to use the fortran routines, use python instead')
            print(e)
            pass

    nsys.make_pair_list(rcut=rcut,nnmax=nnmax)
    anda= np.zeros((len(triplets),na),dtype=float)
    angd= np.array([ dang/2 +dang*ia for ia in range(na) ])
    symbols = nsys.get_symbols()
    lspr = nsys.atoms.neighbors
    for it,t in enumerate(triplets):
        si,sj,sk = t
        for ia in range(natm0):
            if symbols[ia] != si:
                continue
            adfa= adf_atom(ia,dang,rcut,nsys,poss,lspr,symbols,sj,sk)
            for iang in range(na):
                anda[it,iang]= anda[it,iang] +adfa[iang]
    return angd,anda

def adf_average(infiles,dang=1.0,rcut=3.0,triplets=[],
                specorder=None,fortran=False,nnmax=100,
                format=None):
    na= int(180.0/dang)
    aadf= np.zeros((len(triplets),na),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print("[Error] File, {0}, does not exist !!!".format(infname))
            sys.exit()
        #nsys= NAPSystem(fname=infname,specorder=specorder)
        print(' File = ',infname)
        nsys = read(fname=infname,specorder=specorder,format=format)
        angd,df= adf(nsys,dang,rcut,triplets,fortran=fortran,nnmax=nnmax)
        #...NOTE that df is not averaged over the atoms in nsys
        aadf += df
        nsum += 1
    if nsum != 0:
        aadf /= nsum
    return angd,aadf

def write_normal(fname,triplets,na,angd,agr):
    """
    Write out ADF data in normal ADF format.
    """
    outfile= open(fname,'w')
    outfile.write('# 1:theta[i], ')
    for it,t in enumerate(triplets):
        outfile.write(' {0:d}:{1:s}-{2:s}-{3:s},'.format(it+2,*t))
    outfile.write('\n')
    for i in range(na):
        outfile.write(' {0:10.4f}'.format(angd[i]))
        for it,t in enumerate(triplets):
            outfile.write(' {0:11.3e}'.format(agr[it,i]))
        outfile.write('\n')
    outfile.close()
    return None

def write_out4fp(fname,triplets,na,angd,agr,rcut,nperline=6):
    """
    Write out ADF data in general fp.py format.

    Parameters
    ----------
    nperline : int
           Number of data in a line. [default: 6]
    """
    ndat = na*len(triplets)
    data = np.zeros(ndat)
    n = 0
    for it,tri in enumerate(triplets):
        for i in range(na):
            data[n] = agr[it,i]
            n += 1

    with open(fname,'w') as f:
        # cmd = ' '.join(s for s in sys.argv)
        cmd = 'adf.py'
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        cwd = os.getcwd()
        f.write(f'# Output from {cmd:s}\n')
        f.write(f'#   at {cwd:s}\n')
        f.write(f'#   at {now:s}\n')
        f.write('#\n')
        #...Num of data, weight for the data
        f.write(' {0:6d}  {1:7.3f}\n'.format(ndat,1.0))
        j0 = 0
        while True:
            f.write('  '.join('{0:12.4e}'.format(data[j]) for j in range(j0,j0+nperline) if j < ndat))
            f.write('\n')
            j0 += nperline
            if j0 >= ndat:
                break

    return None


def plot_figures(angd,agr,triplets):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(context='talk',style='ticks')
    for i,t in enumerate(triplets):
        plt.plot(angd,agr[i],legend='{0:s}-{1:s}-{2:s}'.format(*t))
    plt.xlabel('Angle (degree)')
    plt.ylabel('ADF')
    plt.savefig("graph_adf.png", format='png', dpi=300, bbox_inches='tight')
    return None


def read_adf(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    nd = len(lines) - 1
    angls = np.zeros(nd)
    trios = []
    adfs = {}
    for il, line in enumerate(lines):
        data = line.split()
        if il == 0:
            for d in data[2:]:
                dd = d.strip(',').split(':')[-1]
                trio = tuple(dd.split('-'))
                trios.append(trio)
                adfs[trio] = np.zeros(nd)
            continue
        angls[il-1] = float(data[0])
        for i, d in enumerate(data[1:]):
            trio = trios[i]
            adfs[trio][il-1] = float(d)
    return angls, adfs


def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),version=__version__)

    infiles= args['INFILE']
    triplets = args['--triplets']
    specorder = [ x for x in args['--specorder'].split(',') ]
    if specorder == ['None']:
        specorder = []
    if triplets == 'None':
        raise ValueError('Triplets must be specified.')
    triplets = [ t.split('-') for t in triplets.split(',') ]
    if len(triplets) == 0:
        raise ValueError('There must be at least one triplet.')
    out4fp = args['--out4fp']
    dang= float(args['-w'])
    drad= np.pi *dang/180.0
    rcut= float(args['-r'])
    nnmax= int(args['--nnmax'])
    sigma= int(args['--gsmear'])
    fortran = args['--fortran']
    fmt = args['--format']
    if fmt == 'None':
        fmt = None
    ofname= args['-o']
    if ofname == 'None':
        ofname = None
    flag_plot= args['--plot']
    nskip = int(args['--skip'])

    if out4fp and ofname is None:
        raise ValueError("Output file name must be specified with option -o.")

    if nskip > len(infiles):
        raise ValueError('NSKIP must be less than num of files given: ',len(infiles))
    if len(infiles) > 1:
        infiles.sort(key=get_key,reverse=True)
    if nskip > 0 and len(infiles) > nskip:
        del infiles[:nskip]
    print(' Number of files to be processed: ',len(infiles))

    na= int(180.0/dang)
    angd,agr= adf_average(infiles,dang=dang,
                          rcut=rcut,triplets=triplets,
                          specorder=specorder,fortran=fortran,nnmax=nnmax,
                          format=fmt)

    if not sigma == 0:
        print(' Gaussian smearing...')
        for it,t in enumerate(triplets):
            agr[it,:] = gsmear(angd,agr[it,:],sigma)

    if flag_plot:
        plot_figures(angd,agr,triplets)
        print('')
        print(' RDF graphes are plotted.')
        print(' Check graph_adf.png')

    #...Regardless ofname, write out.adf in normal format
    write_normal('out.adf',triplets,na,angd,agr)
    #...Format of output (named by ofname) depends on out4fp
    if ofname is not None:
        if out4fp:
            write_out4fp(ofname,triplets,na,angd,agr,rcut)
        else:
            write_normal(ofname,triplets,na,angd,agr)

    print(' --> out.adf')
    if ofname is not None:
        print(' --> {0:s}'.format(ofname))
    return None

if __name__ == "__main__":

    main()
