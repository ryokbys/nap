#!/bin/env python
"""
Compute force constants for the calculation of phonopy.
pmdini should be passed to the program as input and
the system will be extended to include all the interaction within
the cutoff range specified.

OUTPUT:
  * POSCAR
  * FORCE_CONSTANTS
"""
from __future__ import print_function

import os,sys,math,copy
import optparse
import numpy as np

from atom import Atom
from napsys import NAPSystem

################################################# Functions ############

def calc_extention_ratio(hmat,rcut):
    a1= np.zeros(3)
    a2= np.zeros(3)
    a3= np.zeros(3)
    # a1[:]= hmat[:,0]
    # a2[:]= hmat[:,1]
    # a3[:]= hmat[:,2]
    a1[:]= hmat[0,:]
    a2[:]= hmat[1,:]
    a3[:]= hmat[2,:]
    vol= abs(np.dot(a1,np.cross(a2,a3)))
    print(' vol of unit cell=',vol)
    a23= np.cross(a2,a3)
    a31= np.cross(a3,a1)
    a12= np.cross(a1,a2)
    area23= abs(np.sqrt(np.dot(a23,a23)))
    area31= abs(np.sqrt(np.dot(a31,a31)))
    area12= abs(np.sqrt(np.dot(a12,a12)))
    n1= int(rcut/vol*area23)+1
    n2= int(rcut/vol*area31)+1
    n3= int(rcut/vol*area12)+1
    return n1,n2,n3

def cart2h(vc,h):
    """
    Convert Cartesian vector (x,y,z) to the vector reduced in hmat.
    """
    hi= np.linalg.inv(h)
    vr= np.dot(hi.T,vc)
    return vr

################################################ Main routine ##########

usage= '%prog [options] datafile'

parser= optparse.OptionParser(usage=usage)
parser.add_option("-d",dest="displace",type="float",
                  default=0.001,
                  help="size of displacement [Ang] given to each atom.")
parser.add_option("-r",dest="rcut",type="float",
                  default=3.772,
                  help="cutoff radius of the potential.")
parser.add_option("--pmdexec",dest="pmdexec",type="string",
                  default='../pmd/pmd',
                  help="path to the pmd executable.")
(options,args)= parser.parse_args()

displace= options.displace
print(' displacement = ',displace,' Ang.')
rcut= options.rcut
print(' rcut         = ',rcut,' Ang.')
pmdexec= options.pmdexec
infname= args[0]

sys0= NAPSystem()
sys0.read_pmd(infname)
sys0.write_POSCAR()
print(' POSCAR was written.')

natm0= sys0.num_atoms()

h0=np.zeros((3,3))
h0[0]= sys0.a1 *sys0.alc
h0[1]= sys0.a2 *sys0.alc
h0[2]= sys0.a3 *sys0.alc
n1,n2,n3= calc_extention_ratio(h0,rcut)

r1= 2*n1+1
r2= 2*n2+1
r3= 2*n3+1
print(' num of cells in each axis=',r1,r2,r3)
print(' num of atoms in extended system=',natm0*r1*r2*r3)

sysext= NAPSystem()
sysext.set_lattice( sys0.alc
                    ,np.multiply(sys0.a1,r1)
                    ,np.multiply(sys0.a2,r2)
                    ,np.multiply(sys0.a3,r3))
hext=np.zeros((3,3))
hext[0]= sysext.a1 *sysext.alc
hext[1]= sysext.a2 *sysext.alc
hext[2]= sysext.a3 *sysext.alc

#...make extended system corresponding to the phonopy definition
for ia in range(natm0):
    ai0= sys0.atoms[ia]
    for i3 in range(r3):
        for i2 in range(r2):
            for i1 in range(r1):
                ai= Atom()
                ai.set_sid(ai0.sid)
                p1= (ai0.pos[0]+i1)/r1
                p2= (ai0.pos[1]+i2)/r2
                p3= (ai0.pos[2]+i3)/r3
                ai.set_pos(p1,p2,p3)
                sysext.add_atom(ai)
sysext.reset_ids()
natme= sysext.num_atoms()
print ' sysext.num_atoms()=',natme

# for ai in sysext.atoms:
#     print ai.pos

os.system('cp pmdini pmdorig')

#...loop for all atoms in the extended system
fcmat= np.zeros((natme,natme,3,3))
for ia in range(natme):
    print '.',
    sys.stdout.flush()
    ppos= np.array(sysext.atoms[ia].pos)
    for ixyz in range(3):
        fam= np.zeros((natme,3))
        fap= np.zeros((natme,3))
        for kd in (-1,1):
            sysext.atoms[ia].pos[:]= ppos[:]
            dx= np.zeros(3)
            dx[ixyz]= kd*displace
            hdx= cart2h(dx,hext)
            sysext.atoms[ia].pos[:] += hdx
            #...perform pmd here
            sysext.write_pmd('pmdini')
            os.system(pmdexec+' > out.pmd')
            #os.system('cp frc0000 frc0000-{0:1d}-{1:1d}-{2:1d}'.format(ia,ixyz,kd))
            #...obtain forces on atoms
            ff= open('frc0000','r')
            ntmp= int(ff.readline().split()[0])
            for ja in range(natme):
                data= ff.readline().split()
                if kd == -1:
                    fam[ja,0]= float(data[0])
                    fam[ja,1]= float(data[1])
                    fam[ja,2]= float(data[2])
                else:
                    fap[ja,0]= float(data[0])
                    fap[ja,1]= float(data[1])
                    fap[ja,2]= float(data[2])
            ff.close()
        for ja in range(natme):
            for jxyz in range(3):
                fcmat[ia,ja,ixyz,jxyz] += -(fap[ja,jxyz]-fam[ja,jxyz])/(2*displace)

fcfile= open('FORCE_CONSTANTS','w')
fcfile.write(' {0:3d}\n'.format(natme))
for ja in range(natme):
    for ia in range(natme):
        fcfile.write(' {0:3d} {1:3d}\n'.format(ja+1,ia+1))
        fcfile.write(' {0:20.12f} {1:20.12f} {2:20.12f}\n'.format(
            fcmat[ia,ja,0,0]
            ,fcmat[ia,ja,0,1]
            ,fcmat[ia,ja,0,2]))
        fcfile.write(' {0:20.12f} {1:20.12f} {2:20.12f}\n'.format(
            fcmat[ia,ja,1,0]
            ,fcmat[ia,ja,1,1]
            ,fcmat[ia,ja,1,2]))
        fcfile.write(' {0:20.12f} {1:20.12f} {2:20.12f}\n'.format(
            fcmat[ia,ja,2,0]
            ,fcmat[ia,ja,2,1]
            ,fcmat[ia,ja,2,2]))
fcfile.close()
os.system('cp pmdorig pmdini')

