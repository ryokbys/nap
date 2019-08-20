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
  --specorder SPECORDER
              Species order in comma-separated list. [default: None]
  --bvs BVS   Flag for parametrizing fbvs in fitpot, which take one 
              or more parameters into account when converting between
              Morse and fitpot. [default: 0]
"""
from __future__ import print_function

import os
from docopt import docopt

import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "190524"

def pairs2spcs(pairs):
    spcs = []
    for pair in pairs:
        csp1,csp2 = pair
        if csp1 not in spcs:
            spcs.append(csp1)
        if csp2 not in spcs:
            spcs.append(csp2)
    return spcs


def fp2morse(infname,outfname,pairs,ibvs,specorder=[]):
    ds = []
    alps = []
    rs = []
    with open(infname,'r') as f:
        lines = f.readlines()
    
    ndat = int(lines[0].split()[0])
    nprms = len(pairs)*3
    if ibvs == 3:
        nprms = len(pairs)*2
    if ibvs==1:
        nprms += 1
    elif ibvs in (2,3):
        spcs = pairs2spcs(pairs)
        if len(specorder) > 0:
            if len(spcs) != len(specorder):
                raise ValueError('Length of spcs and specorder are different!')
            else:
                spcs = specorder
        nprms += 1 +len(spcs)
    if ndat != nprms:
        raise ValueError('Number of parameters in {0:s} is inconsistent with given pairs.'.format(infname))

    nl = 0
    if ibvs==1:
        nl += 1
        data = lines[nl].split()
        fbvs = float(data[0])
    elif ibvs in (2,3):
        nl += 1
        data = lines[nl].split()
        fbvs = float(data[0])
        rads = {}
        for i in range(len(spcs)):
            nl += 1
            data = lines[nl].split()
            rads[spcs[i]] = float(data[0])

    nl += 1
    if ibvs == 3:
        for i in range(nl,len(lines)):
            line = lines[i]
            data = line.split()
            if (i-nl)%2 == 0:
                alps.append(float(data[0]))
            elif (i-nl)%2 == 1:
                rs.append(float(data[0]))
        try:
            coul = read_params_Coulomb('in.params.Coulomb')
            bvs = coul['bvs']
        except:
            raise ValueError('Could not read in.params.Coulomb.')
        for ip,pair in enumerate(pairs):
            cspi,cspj = pair
            vidi = bvs[cspi][0]
            npqi = bvs[cspi][2]
            vidj = bvs[cspj][0]
            npqj = bvs[cspj][2]
            c = 1.0
            d0 = c*14.4 *abs(vidi*vidj)**(1.0/c)/rs[ip] /np.sqrt(npqi*npqj) /2 /alps[ip]**2
            ds.append(d0)
    else:
        for i in range(nl,len(lines)):
            line = lines[i]
            data = line.split()
            if (i-nl)%3 == 0:
                ds.append(float(data[0]))
            elif (i-nl)%3 == 1:
                alps.append(float(data[0]))
            elif (i-nl)%3 == 2:
                rs.append(float(data[0]))
    
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    D,      alpha,  rmin\n')
        for l in range(len(ds)):
            f.write('  {0:4s}  {1:4s}'.format(pairs[l][0],pairs[l][1]))
            f.write('   {0:7.3f} {1:7.3f} {2:7.3f}\n'.format(ds[l],alps[l],rs[l]))
    print(' Wrote '+outfname)

    if ibvs==1:
        print('')
        print(' The following line should be added to in.params.Coulomb.')
        print('   fbvs    {0:8.4f}'.format(fbvs))
    elif ibvs in (2,3):
        print('')
        print(' The following line should be added to in.params.Coulomb.')
        print('   fbvs    {0:8.4f}'.format(fbvs))
        print(' And rad_bvs in in.params.Coulomb should be replaced by the followings:')
        for s in spcs:
            print('   {0:<3s}  {1:8.4f}'.format(s,rads[s]))
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

def fp2fpc(infname,outfname,pairs):
    aijs = []
    alpijs = []
    bijs = []
    betijs = []
    with open(infname,'r') as f:
        lines = f.readlines()
    
    ndat = int(lines[0].split()[0])
    if ndat != len(pairs)*4 +1:
        raise ValueError('Number of parameters in {0:s} is inconsistent with given pairs.'.format(infname))
    # n = 0
    line1 = lines[1]
    data = line1.split()
    sclchg = float(data[0])
    for i in range(2,len(lines)):
        line = lines[i]
        data = line.split()
        if (i-2)%4 == 0:
            aijs.append(float(data[0]))
        elif (i-2)%4 == 1:
            alpijs.append(float(data[0]))
        elif (i-2)%4 == 2:
            bijs.append(float(data[0]))
        elif (i-2)%4 == 3:
            betijs.append(float(data[0]))
    
    # print(' pairs = ',pairs)
    # print(' aijs  = ',aijs)
    # print(' alpijs= ',alpijs)
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    aij,     alpij,    bij,    betij\n')
        for l in range(len(aijs)):
            f.write('  {0:4s}  {1:4s}'.format(pairs[l][0],pairs[l][1]))
            f.write('   {0:7.3f}  {1:7.3f}'.format(aijs[l],alpijs[l]))
            f.write('   {0:7.3f}  {1:7.3f}\n'.format(bijs[l],betijs[l]))
            
    print(' Wrote '+outfname)

    if os.path.exists('in.params.Coulomb'):
        coul = read_params_Coulomb()
        chgs = coul['charges']
        print(' Charges in in.params.Coulomb should be replaced to the following:')
        for k,v in chgs.items():
            print('   {0:3s}  {1:9.3f}'.format(k,v*sclchg))
    return

def read_params_Coulomb(fname='in.params.Coulomb'):
    """
    Read in.params.Coulomb. But only for fixed or fixed_bvs charges.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    mode = None
    chgtype = None
    chgs = {}
    coul = {}
    bvs = {}
    for line in lines:
        if line[0] in ('#','!'): continue
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if data[0] == 'charges':
            mode = 'charges'
            chgtype = data[1]
            continue
        elif data[0] == 'terms':
            mode = None
            continue
        elif data[0] == 'interactions':
            mode = None
            continue
        elif data[0] == 'fbvs':
            mode = None
            coul['fbvs'] = float(data[1])
            continue
        if mode == 'charges' and chgtype == 'fixed':
            if len(data) == 2:
                csp = data[0]
                chg = float(data[1])
                chgs[csp] = chg
        elif mode == 'charges' and chgtype == 'fixed_bvs':
            if len(data) == 4:
                csp = data[0]
                vid = float(data[1])
                rad = float(data[2])
                npq = int(data[3])
                bvs[csp] = (vid,rad,npq)
    if chgtype == 'fixed':
        coul['charges'] = chgs
    elif chgtype == 'fixed_bvs':
        coul['bvs'] = bvs
    return coul
    
def morse2fp(infname,outfname,ibvs):
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

    nprms = len(params)*3
    if ibvs == 3:
        nprms = len(params)*2

    if ibvs==1:
        nprms += 1
    elif ibvs in (2,3):
        spcs = pairs2spcs(pairs)
        nprms += 1 +len(spcs)
        try:
            coul = read_params_Coulomb('in.params.Coulomb')
            fbvs = coul['fbvs']
            bvs_data = coul['bvs']
        except:
            fbvs = 0.74
            bvs_data = {}
            for s in spcs:
                bvs_data[s] = (1.0, 1.0, 1)
        
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(nprms))
        if ibvs==1:
            f.write(' {0:8.4f}   0.500   1.000  # fbvs\n'.format(0.74))
        elif ibvs in (2,3):
            f.write(' {0:8.4f}   0.500   1.000  # fbvs\n'.format(fbvs))
            for s in spcs:
                vid,rad,npq = bvs_data[s]
                f.write(' {0:8.4f}   0.500   3.000  # rad for {1:s}\n'.format(rad,s))
        for k,v in params.items():
            cspi = k[0]
            cspj = k[1]
            D, alpha, rs = v
            if ibvs == 3:
                f.write(' {0:8.4f}   1.000   3.000  # a for {1:s}-{2:s}\n'.format(alpha,cspi,cspj))
                f.write(' {0:8.4f}   1.000   3.000  # r for {1:s}-{2:s}\n'.format(rs,cspi,cspj))
            else:
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

# def read_params_Coulomb(fname='in.params.Coulomb'):
#     fbvs = -1.0
#     bvs_data = {}
#     with open(fname,'r') as f:
#         lines = f.readlines()
#     mode = None
#     for line in lines:
#         if line[0] in ('#','!'):
#             continue
#         data = line.split()
#         if data[0] == 'charges':
#             mode = 'charges/'+data[1]
#             continue
#         elif data[0] == 'fbvs':
#             fbvs = float(data[1])
#             mode = None
#         if 'charges' in mode:
#             if 'bvs' in mode and len(data) == 4:
#                 csp = data[0]
#                 vid = float(data[1])
#                 rad = float(data[2])
#                 npq = int(data[3])
#                 bvs_data[csp] = (vid,rad,npq)
#     return fbvs, bvs_data

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

def fpc2fp(infname,outfname):
    with open(infname,'r') as f:
        lines = f.readlines()
    params = {}
    pairs = []
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        if data[0].isdigit():
            raise ValueError('This fpc parameter file seems too old.\n'
                             +'Pairs are specified by names, not by integers.')
        cspi = data[0]
        cspj = data[1]
        pairs.append((cspi,cspj))
        aij = float(data[2])
        alpij = float(data[3])
        bij = float(data[4])
        betij = float(data[5])
        params[(cspi,cspj)] = (aij,alpij,bij,betij)

    #...Scale for charges
    sclchg = 1.0
    #...Write in.vars.fitpot file
    with open(outfname,'w') as f:
        f.write('  {0:d}   6.00   3.00\n'.format(len(params)*4+1))
        f.write(' {0:10.3f}   0.100    5.0   # scale for charges\n'.format(sclchg))
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
    ibvs = int(args['--bvs'])
    specorder = args['--specorder'].split(',')
    if specorder[0] == 'None':
        specorder = None
    
    if 'fitpot' in infname and \
       ('Morse' in outfname or 'BMH' in outfname or 'Abell' in outfname
        or 'fpc' in outfname):
        if pairs == 'None':
            raise ValueError('Pairs must be specified.')
        else:
            pairs = [ (pair.split('-')[0],pair.split('-')[1])
                      for pair in pairs.split(',') ]
            print(' Pairs to be extracted:')
            for pair in pairs:
                print('   {0:s}-{1:s}'.format(pair[0],pair[1]))
        if 'Morse' in outfname:
            fp2morse(infname,outfname,pairs,ibvs,specorder)
        elif 'BMH' in outfname:
            fp2bmh(infname,outfname,pairs)
        elif 'Abell' in outfname:
            fp2abell(infname,outfname,pairs)
        elif 'fpc' in outfname:
            fp2fpc(infname,outfname,pairs)
    elif 'Morse' in infname and 'fitpot' in outfname:
        morse2fp(infname,outfname,ibvs)
    elif 'BMH' in infname and 'fitpot' in outfname:
        bmh2fp(infname,outfname)
    elif 'Abell' in infname and 'fitpot' in outfname:
        abell2fp(infname,outfname)
    elif 'fpc' in infname and 'fitpot' in outfname:
        fpc2fp(infname,outfname)
    else:
        msg = 'Input and output file names should include ' \
              +'either fitpot or (Morse).'
        raise ValueError(msg)
    
