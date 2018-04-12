#!/usr/bin/env python
"""
Convert fitpot parameters to in.params.XXX for pmd.

Usage:
  fpvars2prms.py (Morse|BVS) [options] FITPOT_VAR_FILE

Options:
  -h, --help  Show this message and exit.
  -n NSPCS    Number of species. [default: 1]
  --pairs PAIRS
              Specify pairs used in FITPOT_VAR_FILE in the format hyphen-connected
              and comma-separated, e.g.) 1-1,2-1. [default: all] 
  --bvs       Specify BVS parameter set. This sets pairs as 1-1,1-2,...,1-nspeics.
  -i          Inverse conversion, that is in.params.Morse to in.vars.fitpot.
"""
from __future__ import print_function

import os
from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "rev180412"

def ndat2nsp(ndat):
    #...Detect number of species assumed in FITPOT_VAR_FILE
    nd3 = ndat/3
    if nd3 == 1:
        nsp = 1
    elif nd3 == 3:
        nsp = 2
    elif nd3 == 6:
        nsp = 3
    elif nd3 == 10:
        nsp = 4
    elif nd3 == 15:
        nsp = 5
    elif nd3 == 21:
        nsp = 6
    elif nd3 == 28:
        nsp = 7
    elif nd3 == 36:
        nsp = 8
    elif nd3 == 45:
        nsp = 9
    else:
        raise ValueError(' NSP cannot be determined from NDAT.')
    return nsp

def fp2morse(infname,pairs,bvs):
    ds = []
    alps = []
    rs = []
    with open(infname,'r') as f:
        lines = f.readlines()
        ndat = int(lines[0].split()[0])
        if bvs:
            pairs = []
            nsp = ndat/3 +1
            for j in range(nsp):
                ja = j + 1
                if ja == 1:
                    continue
                pairs.append((1,ja))
        elif pairs == 'all':
            nsp = ndat2nsp(ndat)
            pairs = []
            for i in range(nsp):
                ia = i + 1
                for j in range(i,nsp):
                    ja = j + 1
                    pairs.append((ia,ja))
        n = 0
        for i in range(1,len(lines)):
            n += 1
            if n > ndat: break
            ds.append(float(lines[3*(i-1)+1].split()[0]))
            n += 1
            if n > ndat: break
            alps.append(float(lines[3*(i-1)+2].split()[0]))
            n += 1
            if n > ndat: break
            rs.append(float(lines[3*(i-1)+3].split()[0]))

    
    # print(' ds  = ',ds)
    # print(' alps= ',alps)
    # print(' rs  = ',rs)
    with open('in.params.Morse','w') as f:
        f.write('# is, js,  D,      alpha,  rmin\n')
        for l in range(len(ds)):
            f.write(' {0:3d} {1:3d}'.format(pairs[l][0],pairs[l][1]))
            f.write(' {0:7.3f} {1:7.3f} {2:7.3f}\n'.format(ds[l],alps[l],rs[l]))
            
    print(' Wrote in.params.Morse.')
    return

def morse2fp(outfname,bvs):
    with open('in.params.Morse','r') as f:
        lines = f.readlines()
    params = {}
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        isp = int(data[0])
        jsp = int(data[1])
        D = float(data[2])
        alpha = float(data[3])
        rs = float(data[4])
        params[(isp,jsp)] = (D,alpha,rs)
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(len(params)*3))
        for k,v in params.items():
            isp = k[0]
            jsp = k[1]
            D, alpha, rs = v
            f.write(' {0:8.4f}   0.000   8.000\n'.format(D))
            f.write(' {0:8.4f}   1.000   3.000\n'.format(alpha))
            f.write(' {0:8.4f}   1.000   3.000\n'.format(rs))
    print(' Wrote '+outfname)
    return

def get_morse_pairs(nspcs):
    """Compute pairs for BVS-Morse from nspcs, which means only anion-cation pairs.
    """
    pairs = []
    for i in range(2,nspcs+1):
        pairs.append((1,i))
    return pairs

def fp2bvs(infname,nspcs):
    if not os.path.exists('in.params.Coulomb'):
        raise IOError('in.params.Coulomb does not exist.')

    with open(infname,'r') as f:
        lines = f.readlines()
    done_1st_line = False
    il = 0
    while True:
        if lines[il] in ("#", "!"):
            continue
        data = lines[il].split()
        nprms = int(data[0])
        rcut = float(data[1])
        rcut3 = float(data[2])
        done_1st_line = True
        il += 1
        break
    rad_bvs = []
    for i in range(nspcs):
        data = lines[il].split()
        rad_bvs.append(float(data[0]))
        il += 1
    pairs = get_morse_pairs(nspcs)
    prm_morse = {}
    for pair in pairs:
        prm_morse[pair] = []
        for i in range(3):
            data = lines[il].split()
            prm_morse[pair].append(float(data[0]))
            il += 1

    with open('in.params.Morse','w') as f:
        f.write('# is, js,  D,      alpha,  rmin\n')
        for pair in pairs:
            f.write(' {0:3d} {1:3d}'.format(pair[0],pair[1]))
            prm = prm_morse[pair]
            f.write(' {0:7.3f} {1:7.3f} {2:7.3f}\n'.format(prm[0],prm[1],prm[2]))
    print(' Wrote in.params.Morse.')

    os.system('cp in.params.Coulomb in.params.Coulomb.old')
    params = read_params_Coulomb('in.params.Coulomb.old')

    with open('in.params.Coulomb','w') as f:
        for k,v in params.items():
            f.write('{0:s}\n'.format(k))
            if k == 'screened_bvs':
                for iv,entry in enumerate(v):
                    isp,name,vid,rad,npq = entry
                    rad = rad_bvs[iv]
                    f.write(' {0:2d} {1:3s} {2:2d}  {3:6.4f} {4:2d}\n'.format(isp,name,
                                                                              vid,rad,npq))
            elif k == 'interactions':
                for entry in v:
                    isp,jsp = entry
                    f.write(' {0:2d} {1:2d}\n'.format(isp,jsp))
    print(' Wrote in.params.Coulomb')
    return None
            
def read_params_Coulomb(infname='in.params.Coulomb'):
    """
    Read params in in.params.Coulomb and return them in dictionary format.
    """
    if not os.path.exists(infname):
        raise IOError('File does not exists: '+infname)

    keywords = ('screened_bvs',
                'interactions',)
    params = {}
    with open(infname,'r') as f:
        lines = f.readlines()
        mode = None
        for line in lines:
            if line[0] in ('#', '!'):
                continue
            data = line.split()
            kwdline = False
            for kwd in keywords:
                if kwd in line:
                    mode = kwd
                    kwdline = True
                    params[kwd] = []
                    break
            if kwdline:
                continue
            if mode == 'screened_bvs':
                isp = int(data[0])
                name = data[1]
                vid = int(data[2])
                rad = float(data[3])
                npq = int(data[4])
                params[mode].append((isp,name,vid,rad,npq))
            elif mode == 'interactions':
                isp = int(data[0])
                jsp = int(data[0])
                params[mode].append((isp,jsp))
            else:
                raise ValueError('Something wrong with mode: ',mode)
    return params

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['FITPOT_VAR_FILE']
    pairs = args['--pairs']
    bvs = args['--bvs']
    inverse = args['-i']
    nspcs = int(args['-n'])

    if args['Morse']:
        if bvs:
            msg = ' BVS parameters are to be extracted, which means only pairs ' \
                  +'including oxygen are used.'
            print(msg)
        elif pairs == 'all':
            print(' All the pairs are to be extracted.')
        else:
            pairs = [ (pair.split('-')[0],pair.split('-')[1])
                      for pair in pairs.split(',') ]
            print(' Pairs to be extracted:')
            for pair in pairs:
                print('   {0:d}-{1:d}'.format(pair[0],pair[1]))
    
        if not inverse:
            fp2morse(infname,pairs,bvs)
        elif inverse:
            morse2fp(infname,bvs)

    elif args['BVS']:
        """
        The term 'BVS' means that the in.var.fitpot file contains 
        species radius for screened Coulomb and Morse parameters.
        Thus in this case, number of species should be specified and
        the species-1 should be O (anion).
        """
        bvs = True
        if not inverse:
            fp2bvs(infname,nspcs)
        
