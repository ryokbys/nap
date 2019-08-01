#!/bin/env python
"""
Compute 2-body potential energy from 2 atom system.

OUTPUT:
  * out.2body
"""
from __future__ import print_function

import os,sys,math,copy
import optparse
import numpy as np

from atom import Atom
from napsys import NAPSystem

################################################# Functions ############
def get_epot(fname='out.pmd'):
    f= open(fname,'r')
    for line in f.readlines():
        if 'potential energy' in line:
            dat= line.split()
            epot= float(dat[2])
    f.close()
    return epot

def write_banner():
    str="""
          _/_/    _/                        _/
       _/    _/  _/_/_/      _/_/      _/_/_/  _/    _/
          _/    _/    _/  _/    _/  _/    _/  _/    _/
       _/      _/    _/  _/    _/  _/    _/  _/    _/
    _/_/_/_/  _/_/_/      _/_/      _/_/_/    _/_/_/
                                                 _/
                                            _/_/
    """
    print(str)

################################################ Main routine ##########

usage= '%prog [options]'

parser= optparse.OptionParser(usage=usage)
parser.add_option("-n",dest="nsmpl",type="int",
                  default=100,
                  help="number of sampling points.")
parser.add_option("-r",dest="rcut",type="float",
                  default=3.772,
                  help="cutoff radius of the potential.")
parser.add_option("--rmin",dest="rmin",type="float",
                  default=0.5,
                  help="minimum distance in Angstrom.")
parser.add_option("--sid1",dest="sid1",type="int",
                  default=1,
                  help="species ID of atom-1.")
parser.add_option("--sid2",dest="sid2",type="int",
                  default=1,
                  help="species ID of atom-2.")
parser.add_option("--pmdexec",dest="pmdexec",type="string",
                  default='~/src/nap/pmd/pmd',
                  help="path to the pmd executable.")
(options,args)= parser.parse_args()

write_banner()

nsmpl= options.nsmpl
print(' num of points = ',nsmpl)
rcut= options.rcut
print(' rcut          = ',rcut,' Ang.')
rmin= options.rmin
print(' rmin          = ',rmin,' Ang.')
sid1= options.sid1
print(' sid1          = ',sid1)
sid2= options.sid2
print(' sid2          = ',sid2)
pmdexec= options.pmdexec

asys= NAPSystem()
# system size is bigger than 2*rcut
a1= np.array([2.0, 0.0, 0.0])
a2= np.array([0.0, 1.0, 0.0])
a3= np.array([0.0, 0.0, 1.0])
alc= rcut
asys.set_lattice(alc,a1,a2,a3)

atom1= Atom()
atom2= Atom()
atom1.set_pos(0.0,0.0,0.0)
atom1.set_id(1)
atom1.set_sid(sid1)
asys.add_atom(atom1)
atom2.set_pos(0.5,0.0,0.0)
atom2.set_id(2)
atom2.set_sid(sid2)
asys.add_atom(atom2)

hmin= rmin/(2*rcut)
hd  = (0.5-hmin)/nsmpl

os.system('cp pmdini pmdorig')

fout= open('out.2body','w')
for ip in range(nsmpl+1):
    print('.',end='')
    d= hmin +hd*ip
    asys.atoms[1].pos[0] = d
    asys.write_pmd('pmdini')
    os.system(pmdexec+' > out.pmd')
    epot= get_epot('out.pmd')
    fout.write(' {0:12.3f} {1:22.14e}\n'.format(d*2*rcut,epot))

fout.close()
#....restore pmdini
os.system('cp pmdorig pmdini')
print(' program done.')
