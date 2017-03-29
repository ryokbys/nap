#!/bin/env python
"""
Compute 3-body potential energy from 3 atom system.
pmd must be compiled with '-D__PCHECK__'.

OUTPUT:
  * out.3body
"""

import os,sys,math,copy
import optparse
import numpy as np

from atom import Atom
from napsys import NAPSystem

################################################# Functions ############
def get_apot(fname='out.pmd'):
    f= open(fname,'r')
    for line in f.readlines():
        if '3-body term' in line:
            dat= line.split()
            if int(dat[2]) == 1:
                apot= float(dat[3])
    f.close()
    return apot

def write_banner():
    str="""
        _/_/_/    _/                        _/
             _/  _/_/_/      _/_/      _/_/_/  _/    _/
        _/_/    _/    _/  _/    _/  _/    _/  _/    _/
           _/  _/    _/  _/    _/  _/    _/  _/    _/
    _/_/_/    _/_/_/      _/_/      _/_/_/    _/_/_/
                                                 _/
                                            _/_/          
    """
    print str

################################################ Main routine ##########

usage= '%prog [options]'

parser= optparse.OptionParser(usage=usage)
parser.add_option("-n",dest="nsmpl",type="int",
                  default=100,
                  help="number of sampling points.")
parser.add_option("-r",dest="rcut",type="float",
                  default=3.772,
                  help="cutoff radius of the potential.")
parser.add_option("-d",dest="distance",type="float",
                  default=2.0,
                  help="interatomic distance of each pair.")
parser.add_option("--amin",dest="amin",type="float",
                  default= 60,
                  help="minimum angle in degree, [0:180].")
parser.add_option("--amax",dest="amax",type="float",
                  default= 180,
                  help="maximum angle in degree, [0:180].")
parser.add_option("--pmdexec",dest="pmdexec",type="string",
                  default='../pmd/pmd',
                  help="path to the pmd executable.")
(options,args)= parser.parse_args()

write_banner()

nsmpl= options.nsmpl
print ' num of points = ',nsmpl
rcut= options.rcut
print ' rcut          = ',rcut,' Ang.'
amin= options.amin
print ' amin          = ',amin
amax= options.amax
print ' amax          = ',amax
distance= options.distance
print ' distance      = ',distance,' Ang.'
pmdexec= options.pmdexec

asys= NAPSystem()
a1= np.array([2.0, 0.0, 0.0])
a2= np.array([0.0, 2.0, 0.0])
a3= np.array([0.0, 0.0, 1.0])
alc= rcut
asys.set_lattice(alc,a1,a2,a3)

atom1= Atom()
atom1.set_pos(0.0,0.0,0.0)
atom1.set_id(1)
asys.add_atom(atom1)

hd= distance/(alc*2)
atom2= Atom()
atom2.set_pos(hd,0.0,0.0)
atom2.set_id(2)
asys.add_atom(atom2)

atom3= Atom()
atom3.set_pos(0.0,0.0,0.0)
atom3.set_id(3)
asys.add_atom(atom3)

da  = (amax-amin)/nsmpl

fout= open('out.3body','w')
for ip in range(nsmpl+1):
    print '.',
    ang= amin +da*ip
    x= hd*np.cos(ang/180*np.pi)
    y= hd*np.sin(ang/180*np.pi)
    asys.atoms[2].pos[0] = x
    asys.atoms[2].pos[1] = y
    asys.write_pmd('pmdini')
    os.system(pmdexec+' > out.pmd')
    apot= get_apot('out.pmd')
    fout.write(' {0:12.3f} {1:22.14e}\n'.format(ang,apot))

fout.close()
print ' program done.'
