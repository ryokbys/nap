#!/usr/bin/env python
"""
Conversion between fitpot parameter file `in.vars.fitpot` and `in.params.XX` for pmd.
The fitpot parameter file should include the keyword `fitpot` in its name.
Another parameter file type is determined from the file name, for example,
`in.params.Morse` for Morse parameter file.

Usage:
  conversion.py [options] INFILE OUTFILE

Options:
  -h, --help  Show this message and exit.
  --pairs PAIRS
              Specify pairs by hyphen-connected and comma-separated, 
              e.g.) O-Li,O-P. [default: None]
"""
from __future__ import print_function

from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = "190524"


def fp2morse(infname,outfname,pairs):
    ds = []
    alps = []
    rs = []
    with open(infname,'r') as f:
        lines = f.readlines()
    
    ndat = int(lines[0].split()[0])
    if ndat != len(pairs)*3:
        raise ValueError('Number of parameters in {0:s} is inconsistent with given pairs.'.format(infname))
    # n = 0
    for i in range(1,len(lines)):
        line = lines[i]
        data = line.split()
        if (i-1)%3 == 0:
            ds.append(float(data[0]))
        elif (i-1)%3 == 1:
            alps.append(float(data[0]))
        elif (i-1)%3 == 2:
            rs.append(float(data[0]))
        # n += 1
        # if n > ndat: break
        # ds.append(float(lines[3*(i-1)+1].split()[0]))
        # n += 1
        # if n > ndat: break
        # alps.append(float(lines[3*(i-1)+2].split()[0]))
        # n += 1
        # if n > ndat: break
        # rs.append(float(lines[3*(i-1)+3].split()[0]))

    
    # print(' ds  = ',ds)
    # print(' alps= ',alps)
    # print(' rs  = ',rs)
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    D,      alpha,  rmin\n')
        for l in range(len(ds)):
            f.write('  {0:4s}  {1:4s}'.format(pairs[l][0],pairs[l][1]))
            f.write('   {0:7.3f} {1:7.3f} {2:7.3f}\n'.format(ds[l],alps[l],rs[l]))
            
    print(' Wrote '+outfname)
    return

def fp2bmh(infname,outfname,pairs):
    aijs = []
    alpijs = []
    with open(infname,'r') as f:
        lines = f.readlines()
    
    ndat = int(lines[0].split()[0])
    if ndat != len(pairs)*2:
        raise ValueError('Number of parameters in {0:s} is inconsistent with given pairs.'.format(infname))
    # n = 0
    for i in range(1,len(lines)):
        line = lines[i]
        data = line.split()
        if (i-1)%2 == 0:
            aijs.append(float(data[0]))
        elif (i-1)%2 == 1:
            alpijs.append(float(data[0]))
    
    # print(' aijs  = ',aijs)
    # print(' alpijs= ',alpijs)
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    aij,     alpij\n')
        for l in range(len(aijs)):
            f.write('  {0:4s}  {1:4s}'.format(pairs[l][0],pairs[l][1]))
            f.write('   {0:7.3f}  {1:7.3f}\n'.format(aijs[l],alpijs[l]))
            
    print(' Wrote '+outfname)
    return

def fp2abell(infname,outfname,pairs):
    aijs = []
    alpijs = []
    bijs = []
    betijs = []
    with open(infname,'r') as f:
        lines = f.readlines()
    
    ndat = int(lines[0].split()[0])
    if ndat != len(pairs)*4:
        raise ValueError('Number of parameters in {0:s} is inconsistent with given pairs.'.format(infname))
    # n = 0
    for i in range(1,len(lines)):
        line = lines[i]
        data = line.split()
        if (i-1)%4 == 0:
            aijs.append(float(data[0]))
        elif (i-1)%4 == 1:
            alpijs.append(float(data[0]))
        elif (i-1)%4 == 2:
            bijs.append(float(data[0]))
        elif (i-1)%4 == 3:
            betijs.append(float(data[0]))
    
    # print(' aijs  = ',aijs)
    # print(' alpijs= ',alpijs)
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    aij,     alpij,    bij,    betij\n')
        for l in range(len(aijs)):
            f.write('  {0:4s}  {1:4s}'.format(pairs[l][0],pairs[l][1]))
            f.write('   {0:7.3f}  {1:7.3f}'.format(aijs[l],alpijs[l]))
            f.write('   {0:7.3f}  {1:7.3f}\n'.format(bijs[l],betijs[l]))
            
    print(' Wrote '+outfname)
    return

def morse2fp(infname,outfname):
    with open(infname,'r') as f:
        lines = f.readlines()
    params = {}
    pairs = []
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        if data[0].isdigit():
            raise ValueError('This Morse parameter file is too old.\n'
                             +'Pairs are specified by names, not by integers.')
        cspi = data[0]
        cspj = data[1]
        pairs.append((cspi,cspj))
        D = float(data[2])
        alpha = float(data[3])
        rs = float(data[4])
        params[(cspi,cspj)] = (D,alpha,rs)
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(len(params)*3))
        for k,v in params.items():
            cspi = k[0]
            cspj = k[1]
            D, alpha, rs = v
            f.write(' {0:8.4f}   0.100   8.000  # D for {1:s}-{2:s}\n'.format(D,cspi,cspj))
            f.write(' {0:8.4f}   1.000   3.000  # a for {1:s}-{2:s}\n'.format(alpha,cspi,cspj))
            f.write(' {0:8.4f}   1.000   3.000  # r for {1:s}-{2:s}\n'.format(rs,cspi,cspj))
    print(' Wrote '+outfname)
    print('')
    print(' Following lines should be written in in.fitpot.')
    print('{0:-<72}'.format(' '))
    print(' interactions   {0:d}'.format(len(pairs)))
    for p in pairs:
        print('    {0:s}  {1:s}'.format(p[0],p[1]))
    print('{0:-<72}'.format(' '))
    return

def bmh2fp(infname,outfname):
    with open(infname,'r') as f:
        lines = f.readlines()
    params = {}
    pairs = []
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        if data[0].isdigit():
            raise ValueError('This BMH parameter file seems too old.\n'
                             +'Pairs are specified by names, not by integers.')
        cspi = data[0]
        cspj = data[1]
        pairs.append((cspi,cspj))
        aij = float(data[2])
        alpij = float(data[3])
        params[(cspi,cspj)] = (aij,alpij)
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(len(params)*2))
        for k,v in params.items():
            cspi = k[0]
            cspj = k[1]
            aij, alpij = v
            f.write(' {0:10.3f}   0.100  1000.0  # A     for {1:s}-{2:s}\n'.format(aij,cspi,cspj))
            f.write(' {0:10.3f}   0.200   5.000  # alpha for {1:s}-{2:s}\n'.format(alpij,cspi,cspj))
    print(' Wrote '+outfname)
    print('')
    print(' Following lines should be written in in.fitpot.')
    print('{0:-<72}'.format(' '))
    print(' interactions   {0:d}'.format(len(pairs)))
    for p in pairs:
        print('    {0:s}  {1:s}'.format(p[0],p[1]))
    print('{0:-<72}'.format(' '))
    return

def abell2fp(infname,outfname):
    with open(infname,'r') as f:
        lines = f.readlines()
    params = {}
    pairs = []
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        if data[0].isdigit():
            raise ValueError('This Abell parameter file seems too old.\n'
                             +'Pairs are specified by names, not by integers.')
        cspi = data[0]
        cspj = data[1]
        pairs.append((cspi,cspj))
        aij = float(data[2])
        alpij = float(data[3])
        bij = float(data[4])
        betij = float(data[5])
        params[(cspi,cspj)] = (aij,alpij,bij,betij)
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(len(params)*4))
        for k,v in params.items():
            cspi = k[0]
            cspj = k[1]
            aij,alpij,bij,betij = v
            f.write(' {0:10.3f}   0.100  1000.0  # A     for {1:s}-{2:s}\n'.format(aij,cspi,cspj))
            f.write(' {0:10.3f}   0.200   5.000  # alpha for {1:s}-{2:s}\n'.format(alpij,cspi,cspj))
            f.write(' {0:10.3f}   0.100  1000.0  # B     for {1:s}-{2:s}\n'.format(bij,cspi,cspj))
            f.write(' {0:10.3f}   0.200   5.000  # beta  for {1:s}-{2:s}\n'.format(betij,cspi,cspj))
    print(' Wrote '+outfname)
    print('')
    print(' Following lines should be written in in.fitpot.')
    print('{0:-<72}'.format(' '))
    print(' interactions   {0:d}'.format(len(pairs)))
    for p in pairs:
        print('    {0:s}  {1:s}'.format(p[0],p[1]))
    print('{0:-<72}'.format(' '))
    return


if __name__ == "__main__":

    args = docopt(__doc__,version=__version__)
    pairs = args['--pairs']
    infname = args['INFILE']
    outfname = args['OUTFILE']
    
    if 'fitpot' in infname and \
       ('Morse' in outfname or 'BMH' in outfname or 'Abell' in outfname):
        if pairs == 'None':
            raise ValueError('Morse parameter pairs must be specified.')
        else:
            pairs = [ (pair.split('-')[0],pair.split('-')[1])
                      for pair in pairs.split(',') ]
            print(' Pairs to be extracted:')
            for pair in pairs:
                print('   {0:s}-{1:s}'.format(pair[0],pair[1]))
        if 'Morse' in outfname:
            fp2morse(infname,outfname,pairs)
        elif 'BMH' in outfname:
            fp2bmh(infname,outfname,pairs)
        elif 'Abell' in outfname:
            fp2abell(infname,outfname,pairs)
    elif 'Morse' in infname and 'fitpot' in outfname:
        morse2fp(infname,outfname)
    elif 'BMH' in infname and 'fitpot' in outfname:
        bmh2fp(infname,outfname)
    elif 'Abell' in infname and 'fitpot' in outfname:
        abell2fp(infname,outfname)
    else:
        msg = 'Input and output file names should include ' \
              +'either fitpot or (Morse).'
        raise ValueError(msg)
    
