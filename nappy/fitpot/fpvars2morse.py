#!/usr/bin/env python
"""
Convert fitpot parameters to in.params.Morse for pmd.

Usage:
  fp2morse.py [options] FITPOT_VAR_FILE

Options:
  -h, --help  Show this message and exit.
  --pairs PAIRS
              Specify pairs used in FITPOT_VAR_FILE in the format hyphen-connected
              and comma-separated, e.g.) 1-1,2-1. [default: all] 
  --bvs       Specify BVS parameter set. This sets pairs as 1-1,1-2,...,1-nspeics.
  -i          Inverse conversion, that is in.params.Morse to in.vars.fitpot.
"""
from __future__ import print_function

from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = ""

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


if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['FITPOT_VAR_FILE']
    pairs = args['--pairs']
    bvs = args['--bvs']
    inverse = args['-i']
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
