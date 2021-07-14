#!/usr/bin/env python
"""
Count coordination number of given center-neighbor pairs.
Output histogram of coodination numbers for each center-neighbor pairs.

Usage:
    coordination.py [options] INFILE

Options:
    -h, --help  Show this help message and exit.
    -c CUTOFF   Cutoff radius for neighbors. [default: 3.0]
    -o OUT      Output file name. [default: out.coord]
    --pairs PAIRS
                Center-neighbor pairs for which coordination is to be counted. 
                Hyphen-connected and comma-separated such as, Si-Si,Si-O. [default: None]
    --maxcoord MAXCOORD
                Max coordination number to be output for histogram. [default: 6]
    --out4fp    Output in a format suitable for fp.
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
from datetime import datetime
from nappy.io import read

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
        if si not in centers:
            continue
        icindex = centers.index(si)
        lspri = nsys.get_atom_attr(ia,'lspr')
        for p in pairs:
            counter[p] = 0
        for jj,ja in enumerate(lspri):
            sj = symbols[ja]
            if (si,sj) in counter.keys():
                counter[(si,sj)] += 1
        for p in pairs:
            nump = counter[p]
            if nump <= maxcoord:
                coords[p][nump] += 1.0
        
        # for pair in pairs:
        #     coord_dict = {}
        #     if not nsys.get_atom_attr(ia,'sid') == pair[0]:
        #         continue
        #     coord_dict['ia'] = ia
        #     coord_dict['pair'] = pair
        #     cnum = 0
        #     for jj in range(len(lspri)):
        #         ja = lspri[jj]
        #         if not nsys.get_atom_attr(ja,'sid') == pair[1]:
        #             continue
        #         d = nsys.get_distance(ia,ja)
        #         if d < rcut:
        #             cnum += 1
        #     coord_dict['cnum'] = cnum
        #     if len(coord_dict) > 0:
        #         coordnums.append(coord_dict)
    return coords

def main(args):
    outfname = args['-o']
    infname = args['INFILE']
    rcut = float(args['-c'])
    maxcoord = int(args['--maxcoord'])
    out4fp = args['--out4fp']

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
    
    try:
        nsys = read(fname=infname)
    except Exception:
        raise

    coords = count(nsys,pairs,rcut=rcut,maxcoord=maxcoord)


    cmd = ' '.join(s for s in sys.argv)
    if out4fp:
        nperline = 6
        with open(outfname,'w') as f:
            f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('#  {0:s}\n'.format(cmd))
            ndat = len(pairs)*maxcoord
            f.write('  {0:6d}  {1:7.3f}\n'.format(ndat,1.0))
            data = np.zeros(ndat)
            inc = 0
            for p in pairs:
                si,sj = p
                coord = coords[p]
                for i in range(1,maxcoord+1):
                    data[inc] = coord[i]
                    inc += 1
            j0 = 0
            while True:
                f.write(' '.join('{0:8.1f}'.format(data[j]) for j in range(j0,j0+nperline) if j < ndat))
                f.write('\n')
                j0 += nperline
                if j0 >= ndat:
                    break
        print(' Wrote ',outfname)
    else:
        with open(outfname,'w') as f:
            f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('#  {0:s}\n'.format(cmd))
            f.write('#   center, neighbor,  coordination histogram \n')
            for p in pairs:
                si,sj = p
                coord = coords[p]
                f.write('  {0:s}  {1:s} '.format(si,sj))
                for i in range(1,maxcoord+1):
                    f.write(' {0:8.1f}'.format(coord[i]))
                f.write('\n')
        print(' Wrote ',outfname)
    
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    main(args)
