#!/usr/bin/env python
"""
Reduce samples that have too high energies by comparing
between the same group of samples.
The group is defined according to the name of directory before the last 5 digits.
For example, the directory `smpl_XX_YYYYYY_#####` where `#####` is the
last 5 digits and the group name would be `smpl_XX_YYYYY`.

Usage:
  reduce_high_energy_samples.py [options] DIRS...

Options:
  -h,--help  Show this message and exit.
  -o OUT     Output file name. [default: out.high_energy_samples]
  --threshold=THRESHOLD
             Threshold of energy/atom that determines high energy samples.
             [default: 1.0]
"""
import os,sys
from docopt import docopt
from datetime import datetime

from nappy.napsys import NAPSystem

__author__ = "RYO KOBAYASHI"
__version__ = "160727"

def get_obsolete_dirname():
    prefix = "obsolete_"
    today = datetime.today()
    return prefix+today.strftime("%y%m%d")

def get_groups(smpldirs):
    groups = {}
    ns = len(smpldirs)
    if ns < 100:
        ms = 1
    else:
        ms = ns/100
    for i,s in enumerate(smpldirs):
        if i%ms == 0:
            print('.',end=".")
        try:
            with open(s+'/erg.ref','r') as f:
                erg = float(f.readline())
        except:
            print('Failed to read erg.ref, so skip '+s)
            continue
        key = s[:-6]
        if not key in groups:
            groups[key] = []
        groups[key].append([s,erg])
    print('')
    return groups

def get_list_high_energy(gsmpls,threshold):
    emin = 1e+30
    highsmpls = []
    ergs = []
    for i,s in enumerate(gsmpls):
        smpldir = s[0]
        erg = s[1]
        #atoms = read(smpldir+'/POSCAR',format='vasp')
        atoms = NAPSystem(fname=smpldir+"/pos",format='pmd')
        natm = atoms.num_atoms()
        erg /= natm
        ergs.append(erg)
        emin = min(erg,emin)
    for i,s in enumerate(gsmpls):
        smpldir = s[0]
        erg = ergs[i]
        if erg-emin > threshold:
            highsmpls.append(smpldir)
    return highsmpls

if __name__ == "__main__":

    args = docopt(__doc__)
    smpldirs = args['DIRS']
    outfname = args['-o']
    threshold = float(args['--threshold'])
    
    print('grouping samples...')
    groups = get_groups(smpldirs)

    print('looking for high-energy samples...')
    highsmpls = []
    for g,smpls in groups.items():
        print('.',end='')
        highsmpls.extend(get_list_high_energy(smpls,threshold))
    print('')

    with open(outfname,'w') as f:
        for s in highsmpls:
            f.write(s+'\n')
    print('number of samples to be reduced = ',len(highsmpls))
    print('check '+outfname+' and run the following commands:')
    print('')
    # obsdir = get_obsolete_dirname()
    # print('  mkdir '+obsdir)
    # print('  for d in `cat '+outfname+'`; do mv $d '+obsdir 
    #       +'/; done')
    # print('')

