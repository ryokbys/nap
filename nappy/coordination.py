#!/usr/bin/env python
"""
Count coordination number of given center-neighbor pairs.
Output histogram of coodination numbers for each center-neighbor pairs.

Usage:
    coordination.py [options] INFILES [INFILES...]

Options:
    -h, --help  Show this help message and exit.
    -c CUTOFF   Cutoff radius for neighbors. [default: 3.0]
    -o OUT      Output file name. [default: out.coord]
    --pairs PAIRS
                Center-neighbor pairs for which coordination is to be counted. 
                Hyphen-connected and comma-separated such as, Si-Si,Si-O. [default: None]
    --maxcoord MAXCOORD
                Max coordination number to be output for histogram. [default: 6]
    --out4fp    Output in the optzer independent-data format to data.coord file.
    --skip=NSKIP
                Skip first NSKIP steps from the statistics. [default: 0]
"""

import os,sys
from docopt import docopt
import numpy as np
from datetime import datetime
from nappy.io import read
from nappy.common import get_key

def pairs_from_pairs_str(pairs_str):
    pairs = []
    print('pairs_str=',pairs_str)
    for pairstr in pairs_str:
        if len(pairstr) < 2:
            continue
        print('  pairstr=',pairstr)
        pair = [ int(i) for i in pairstr.split('-') ]
        pairs.append(pair)
    return pairs

def count(nsys,pairs,rcut=3.0,maxcoord=6):
    """
    Count coordination numbers of neigh around center..
    """
    #print('Making pair list...')
    nsys.make_pair_list(rcut=rcut)
    isps_treated = []
    for pair in pairs:
        if not pair[0] in isps_treated:
            isps_treated.append(pair[0])
    #...Count coordination numbers for given pairs
    coordnums = []
    symbols = nsys.get_symbols()
    centers = [ pair[0] for pair in pairs ]
    neighbors = [ pair[1] for pair in pairs ]
    coords = {}
    counter = {}
    for p in pairs:
        coords[p] = np.zeros(maxcoord+1)  # 0,1,...,maxcoord
        counter[p] = 0
        
        for ia in range(nsys.num_atoms()):
            si = symbols[ia]
            if si not in centers or si != p[0]:
                continue
            # icindex = centers.index(si)
            lspri = nsys.get_atom_attr(ia,'neighbors')
            # for p in pairs:
            #     counter[p] = 0
            counter[p] = 0
            for jj,ja in enumerate(lspri):
                sj = symbols[ja]
                if sj == p[1]:
                    counter[p] += 1
            coords[p][counter[p]] += 1
        # for p in pairs:
        #     nump = counter[p]
        #     if nump <= maxcoord:
        #         coords[p][nump] += 1.0
    return coords

def main(args):
    outfname = args['-o']
    infiles = args['INFILES']
    rcut = float(args['-c'])
    maxcoord = int(args['--maxcoord'])
    out4fp = args['--out4fp']
    nskip = int(args['--skip'])

    pairstr = args['--pairs'].split(',')
    pairs = []
    for pair in pairstr:
        try: 
            ci,cj = pair.split('-')
        except Exception:
            raise
        pairs.append((ci,cj))
    if len(pairs) == 0:
        raise ValueError('len(pairs) == 0 !!')

    if len(infiles) > 1:
        infiles.sort(key=get_key,reverse=True)
    del infiles[:nskip]
    if len(infiles) < 1:
        raise ValueError('No input files to be processed.')
    print(' Number of files to be processed = ',len(infiles))

    try:
        nsys0 = read(fname=infiles[0])
    except Exception:
        raise

    n_per_spec0 = nsys0.natm_per_species()
    specorder = nsys0.specorder
    n_per_spec = {}
    for i,s in enumerate(specorder):
        n_per_spec[s] = n_per_spec0[i]

    coords = {}
    for infile in infiles:
        print(f'   {infile:s}...')
        nsys = read(fname=infile)
        coords_tmp = count(nsys,pairs,rcut=rcut,maxcoord=maxcoord)
        for k,v in coords_tmp.items():
            if k in coords.keys():
                coords[k] += v
            else:
                coords[k] = v

    #...Normalize
    for pair in coords.keys():
        si = pair[0]
        nsi = n_per_spec[si]
        coords[pair] /= nsi*len(infiles)

    cmd = ' '.join(s for s in sys.argv)
    with open(outfname,'w') as f:
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        f.write(f'# Output from {cmd:s}\n')
        f.write(f'#   at {now:s}\n')
        f.write('# Num of files = {0:d}\n'.format(len(infiles)))
        f.write('# Num_per_species = ')
        for s,n in n_per_spec.items():
            f.write(f'  {s:s}:{n:d},')
        f.write('\n')
        f.write('#   center, neighbor,  coordination histogram \n')
        for p in pairs:
            si,sj = p
            coord = coords[p]
            f.write('  {0:s}  {1:s} '.format(si,sj))
            for i in range(0,maxcoord+1):
                f.write(' {0:8.4f}'.format(coord[i]))
            f.write('\n')
    print(f' Wrote {outfname:s}')
    
    if out4fp:
        nperline = 6
        with open('data.coord','w') as f:
            now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            f.write(f'# Output from {cmd:s}\n')
            f.write(f'#   at {now:s}\n')
            f.write('# Num of files = {0:d}\n'.format(len(infiles)))
            f.write('# Num_per_species = ')
            for s,n in n_per_spec.items():
                f.write(f'  {s:s}:{n:d},')
            f.write('\n')
            f.write('# datatype: independent\n')
            f.write('#\n')
            ndat = len(pairs)*(maxcoord+1)
            f.write('  {0:6d}  {1:7.3f}\n'.format(ndat,1.0))
            data = np.zeros(ndat)
            inc = 0
            for p in pairs:
                si,sj = p
                coord = coords[p]
                for i in range(0,maxcoord+1):
                    data[inc] = coord[i]
                    inc += 1
            # j0 = 0
            # while True:
            #     f.write(' '.join('{0:8.1f}'.format(data[j]) for j in range(j0,j0+nperline) if j < ndat))
            #     f.write('\n')
            #     j0 += nperline
            #     if j0 >= ndat:
            #         break
            for i in range(inc):
                d = data[i]
                dw = max(d,0.1)
                f.write(f'  {d:6.4f}  {dw:6.4f}\n')
        print(' Wrote data.coord')
    
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    main(args)
