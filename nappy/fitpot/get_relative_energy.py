#!/usr/bin/env python
"""
Get relative energy of the system of DIR w.r.t. the energies ERG.
There must be erg.ref and POSCAR in the DIR.
Relative energy of the system and that per atom will be printed.

Usage:
  get_relative_energy.py [options] DIR

Options:
  -h, --help  Show this message and exit.
  --specorder=SPECORDER
              Order of species. [default: Al,Mg,Si]
  --echems=ECHEMS
              Chemical energy of species.
              [default: -3.5469,-1.5055,-4.6169]
"""
import os,sys
from docopt import docopt
from ase.io import read

author = "RYO KOBAYASHI"
version = ""

def get_relative_energy(atoms,echems,erg):
    symbols = atoms.get_chemical_symbols()
    for s in symbols:
        erg -= echems[s]
    return erg

if __name__ == "__main__":

    args = docopt(__doc__)
    smpldir = args['DIR']
    specorder = args['--specorder'].split(',')
    echems = {}
    for i,e in enumerate(args['--echems'].split(',')):
        echems[specorder[i]] = float(e)

    if len(echems) != len(specorder):
        raise RuntimeError('Num of echems should be {}'.format(len(specorder)))

    # for s,e in echems.items():
    #     print('{0:s}: {1:10.4f}'.format(s,e))
    
    atoms = read(smpldir+'/POSCAR',format="vasp")
    with open(smpldir+'/erg.ref','r') as f:
        erg = float(f.readline())
    erg = get_relative_energy(atoms,echems,erg)
    print('{0:15.7f} {1:15.7f}'.format(erg,erg/len(atoms)))
