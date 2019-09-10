#!/usr/bin/env python
"""
Calculate the angular distribution function (ADF) from files.
Take an average over atoms in a file or files.

Usage:
    adf.py [options] INFILE [INFILE...]

Options:
    -h, --help  Show this help message and exit.
    -w DEG      Width of the angular degree. [default: 1.0]
    -r RCUT     Cutoff radius of the bonding pair. [default: 3.0]
    --gsmear=SIGMA
                Width of Gaussian smearing, zero means no smearing. [default: 0]
    --triplets=TRIPLETS
                Triplets whose angles are to be computed. Three species should be specified connected by hyphen,
                and separated by comma, e.g.) P-O-O,Li-O-O. [default: None]
    -o OUT      Output file name [default: out.adf]
    --no-average
                Not to take average over files.
    --plot      Plot figures. [default: False]
"""
from __future__ import print_function

import os,sys
import numpy as np
from docopt import docopt

from nappy.napsys import NAPSystem
from nappy.gaussian_smear import gsmear


def norm(vector):
    norm= 0.0
    for e in vector:
        norm += e*e
    return np.sqrt(norm)

def adf_atom(ia,dang,rcut,nsys,symbols,sj,sk):
    """
    Compute number of atoms in the every range of angle [0:180].
    """
    na= int(180.0/dang) +1
    hmat= nsys.get_hmat()
    nda= np.zeros(na,dtype=np.int)
    natm= nsys.natm
    rcut2= rcut*rcut
    pi= nsys.poss[ia]
    for ji in range(nsys.nlspr[ia]):
        ja= nsys.lspr[ia,ji]
        if ja == ia:
            continue
        sji = symbols[ja]
        if sji not in (sj,sk):
            continue
        pj= nsys.poss[ja]
        pij= pj-pi
        pij= pij -np.round(pij)
        vij= np.dot(hmat,pij)
        rij2= np.dot(vij,vij)
        if rij2 >= rcut2:
            continue
        rij= np.sqrt(rij2)
        for ki in range(nsys.nlspr[ia]):
            ka= nsys.lspr[ia,ki]
            if ka == ia or ka <= ja:
                continue
            ski = symbols[ka]
            if set((sji,ski)) != set((sj,sk)):
                continue
            pk= nsys.poss[ka]
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

def adf(nsys,dang,rcut,triplets):

    natm0= nsys.natm

    n1,n2,n3= nsys.get_expansion_num(2.0*rcut)
    if not (n1==1 and n2==1 and n3==1):
        print(' system to be repeated, n1,n2,n3=',n1,n2,n3)
        nsys.repeat(n1,n2,n3)
    nsys.assign_pbc()

    nsys.make_pair_list(rcut=rcut)

    na= int(180.0/dang)+1

    anda= np.zeros((len(triplets),na),dtype=np.float)
    angd= np.array([ dang*ia for ia in range(na) ])
    symbols = nsys.get_symbols()
    for it,t in enumerate(triplets):
        si,sj,sk = t
        for ia in range(natm0):
            if symbols[ia] != si:
                continue
            adfa= adf_atom(ia,dang,rcut,nsys,symbols,sj,sk)
            for iang in range(na):
                anda[it,iang]= anda[it,iang] +adfa[iang]
    return angd,anda,natm0

def adf_average(infiles,dang=1.0,rcut=3.0,triplets=[],no_average=False):
    na= int(180.0/dang) +1
    aadf= np.zeros((len(triplets),na),dtype=float)
    nsum= 0
    for infname in infiles:
        if not os.path.exists(infname):
            print("[Error] File, {0}, does not exist !!!".format(infname))
            sys.exit()
        nsys= NAPSystem(fname=infname,)
        print(' File = ',infname)
        angd,df,n= adf(nsys,dang,rcut,triplets)
        aadf += df
        nsum += n

    if not no_average:
        aadf /= nsum
    return angd,aadf

def plot_figures(angd,agr,triplets):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(context='talk',style='ticks')
    for i,t in enumerate(triplets):
        plt.plot(angd,agr[i],legend='{0:s}-{1:s}-{2:s}'.format(*t))
    plt.xlabel('Angle (degree)')
    plt.ylabel('ADF')
    plt.savefig("graph_adf.png", format='png', dpi=300, bbox_inches='tight')

################################################## main routine

if __name__ == "__main__":

    args= docopt(__doc__)
    
    infiles= args['INFILE']
    triplets = args['--triplets']
    if triplets == 'None':
        raise ValueError('Triplets must be specified.')
    triplets = [ t.split('-') for t in triplets.split(',') ]
    if len(triplets) == 0:
        raise ValueError('There must be at least one triplet.')
    dang= float(args['-w'])
    drad= np.pi *dang/180.0
    rcut= float(args['-r'])
    sigma= int(args['--gsmear'])
    no_average = args['--no-average']
    ofname= args['-o']
    flag_plot= args['--plot']

    na= int(180.0/dang) +1
    angd,agr= adf_average(infiles,dang=dang,
                          rcut=rcut,triplets=triplets,
                          no_average=no_average)

    if not sigma == 0:
        print(' Gaussian smearing...')
        for it,t in enumerate(triplets):
            agr[it,:] = gsmear(angd,agr[it,:],sigma)

    if flag_plot:
        plot_figures(angd,agr,triplets)
        print('')
        print(' RDF graphes are plotted.')
        print(' Check graph_adf.png')
        
    outfile= open(ofname,'w')
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
    print(' Wrote '+ofname)
