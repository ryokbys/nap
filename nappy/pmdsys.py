#!/bin/local/env python
"""
System information used in pmd, including cell information, lattice constant,
number of atoms, atom species, atom positions.

Users can use `pmdsys.py` as a converter program as described below.

Usage:
    pmdsys.py [options] INFILE OUTFILE

Options:
    -h, --help  Show this help message and exit.
    --in-format=INFMT
                Format of the input file. [default: None]
    --out-format=OUTFMT
                Format of the output file. [default: None]
"""

import math
import sys,copy,re
from datetime import datetime

import numpy as np
from docopt import docopt

from atom import Atom

#...constants
_maxnn= 100
_file_formats= ('pmd','akr','POSCAR')

class PMDSystem(object):
    """
    Contains cell information and atoms, and provides some functionalities.
    """

    def __init__(self,fname=None,ffmt=None):
        self.a1= np.zeros(3)
        self.a2= np.zeros(3)
        self.a3= np.zeros(3)
        self.atoms= []

        if not fname == None:
            ftype= self.parse_filename(fname)
            if ftype == 'pmd':
                self.read_pmd(fname)
            elif ftype == 'akr':
                self.read_akr(fname)
            elif ftype == 'POSCAR':
                self.read_POSCAR(fname)

    def set_lattice(self,alc,a1,a2,a3):
        self.alc= alc
        self.a1= a1
        self.a2= a2
        self.a3= a3

    def add_atom(self,atom):
        self.atoms.append(atom)

    def reset_ids(self):
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)

    def num_atoms(self):
        return len(self.atoms)

    def volume(self):
        return self.alc**3 *np.abs(np.dot(self.a1,np.cross(self.a2,self.a3)))

    def write(self,fname="pmd0000"):
        ftype= self.parse_filename(fname)
        if ftype == "pmd":
            self.write_pmd(fname)
        elif ftype == "akr":
            self.write_akr(fname)
        elif ftype == "POSCAR":
            self.write_POSCAR(fname)

    def read(self,fname="pmd0000"):
        ftype= self.parse_filename(fname)
        if ftype == "pmd":
            self.read_pmd(fname)
        elif ftype == "akr":
            self.read_akr(fname)
        elif ftype == "POSCAR":
            self.read_POSCAR(fname)

    def read_pmd(self,fname='pmd0000'):
        f=open(fname,'r')
        # 1st: lattice constant
        self.alc= float(f.readline().split()[0])
        # 2nd-4th: cell vectors
        # for i in range(3):
        #     data= f.readline().split()
        #     self.a1[i]= float(data[0])
        #     self.a2[i]= float(data[1])
        #     self.a3[i]= float(data[2])
        self.a1= np.array([float(x) for x in f.readline().split()])
        self.a2= np.array([float(x) for x in f.readline().split()])
        self.a3= np.array([float(x) for x in f.readline().split()])
        # 5th-7th: velocity of cell vectors
        tmp= f.readline().split()
        tmp= f.readline().split()
        tmp= f.readline().split()
        # 8st: num of atoms
        natm= int(f.readline().split()[0])
        # 9th-: atom positions
        self.atoms= []
        for i in range(natm):
            data= [float(x) for x in f.readline().split()]
            ai= Atom()
            ai.decode_tag(data[0])
            ai.set_pos(data[1],data[2],data[3])
            ai.set_vel(data[4],data[5],data[6])
            self.atoms.append(ai)
        f.close()

    def write_pmd(self,fname='pmd0000'):
        f=open(fname,'w')
        # lattice constant
        f.write(" {0:15.9f}\n".format(self.alc))
        # cell vectors
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(self.a1[0],\
                                                          self.a1[1],\
                                                          self.a1[2]))
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(self.a2[0],\
                                                          self.a2[1],\
                                                          self.a2[2]))
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(self.a3[0],\
                                                          self.a3[1],\
                                                          self.a3[2]))
        # velocities of cell vectors
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
        f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
        # num of atoms
        f.write(" {0:10d}\n".format(len(self.atoms)))
        # atom positions
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)
            f.write(" {0:22.14e} {1:19.15f} {2:19.15f} {3:19.15f}".format(ai.tag(), \
                                                            ai.pos[0],\
                                                            ai.pos[1],\
                                                            ai.pos[2])
                    +"  {0:11.7f}  {1:11.7f}  {2:11.7f}".format(ai.vel[0], 
                                                                ai.vel[1],
                                                                ai.vel[2])
                    +"  {0:.1f}  {1:.1f}".format(0.0, 0.0)
                    +"  {0:.1f}  {1:.1f}  {2:.1f}".format(0.0, 0.0, 0.0)
                    +"  {0:.1f}  {1:.1f}  {2:.1f}".format(0.0, 0.0, 0.0)
                    +"  {0:.1f}  {1:.1f}  {2:.1f}".format(0.0, 0.0, 0.0)
                    +"\n")
        f.close()

    def read_POSCAR(self,fname='POSCAR'):
        with open(fname,'r') as f:
            # 1st line: comment
            f.readline()
            # 2nd: lattice constant
            self.alc= float(f.readline().split()[0])
            # 3rd-5th: cell vectors
            self.a1= np.array([float(x) for x in f.readline().split()])
            self.a2= np.array([float(x) for x in f.readline().split()])
            self.a3= np.array([float(x) for x in f.readline().split()])
            # 6th: species names or number of each species
            buff= f.readline().split()
            if not buff[0].isdigit():
                buff= f.readline().split()
            num_species= np.array([ int(n) for n in buff])
            natm= 0
            for n in num_species:
                natm += n
            #print("Number of atoms = {0:5d}".format(natm))
            # 7th or 8th line: comment
            c7= f.readline()
            if c7[0] in ('s','S'):
                c8= f.readline()
            # Atom positions hereafter
            self.atoms= []
            for i in range(natm):
                buff= f.readline().split()
                ai= Atom()
                sid= 1
                m= 0
                for n in num_species:
                    m += n
                    if i < m:
                        break
                    sid += 1
                ai.set_sid(sid)
                ai.set_pos(float(buff[0]),float(buff[1]),float(buff[2]))
                ai.set_vel(0.0,0.0,0.0)
                self.atoms.append(ai)
                
            

    def write_POSCAR(self,fname='POSCAR'):
        f=open(fname,'w')
        f.write(' Generated by pmdsys.py at {0}.\n'.format(datetime.now().strftime('%Y-%m-%d')))
        # lattice vector
        f.write(' {0:15.7f}\n'.format(self.alc))
        # cell vectors
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a1[0],\
                                                          self.a1[1],\
                                                          self.a1[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a2[0],\
                                                          self.a2[1],\
                                                          self.a2[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a3[0],\
                                                          self.a3[1],\
                                                          self.a3[2]))
        # count num of atoms per specie
        num_species= self.num_species()
        for n in num_species:
            f.write(' {0:4d}'.format(n))
        f.write('\n')
        # comments
        f.write('Selective dynamics\n')
        f.write('Direct\n')
        # atom positions
        for ai in self.atoms:
            f.write(' {0:15.7f} {1:15.7f} {2:15.7f} T T T\n'.format(
                ai.pos[0],ai.pos[1],ai.pos[2]))
        f.close()

    def num_species(self):
        num_species= []
        max_nsp= 0
        for ai in self.atoms:
            max_nsp= max(max_nsp,ai.sid)
        for i in range(max_nsp):
            num_species.append(0)
        for ai in self.atoms:
            num_species[ai.sid-1] += 1
        return num_species

    def read_akr(self,fname='akr0000'):
        f=open(fname,'r')
        # 1st: lattice constant
        self.alc= float(f.readline().split()[0])
        # 2nd-4th: cell vectors
        for i in range(3):
            data= f.readline().split()
            self.a1[i]= float(data[0])
            self.a2[i]= float(data[1])
            self.a3[i]= float(data[2])
        # 5th: num of atoms
        natm= int(f.readline().split()[0])
        # 9th-: atom positions
        self.atoms= []
        for i in range(natm):
            data= [float(x) for x in f.readline().split()]
            ai= Atom()
            ai.set_sid(data[0])
            ai.set_pos(data[1],data[2],data[3])
            ai.set_vel(data[4],data[5],data[6])
            self.atoms.append(ai)
        f.close()

    def write_akr(self,fname='akr0000'):
        f=open(fname,'w')
        # lattice constant
        f.write(" {0:12.4f}\n".format(self.alc))
        # cell vectors
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a1[0],\
                                                          self.a2[0],\
                                                          self.a3[0]))
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a1[1],\
                                                          self.a2[1],\
                                                          self.a3[1]))
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a1[2],\
                                                          self.a2[2],\
                                                          self.a3[2]))
        # num of atoms
        f.write(" {0:10d} {1:4d} {2:4d} {3:4d}\n".format(len(self.atoms),3,0,0))
        # atom positions
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)
            f.write(" {0:4d} {1:10.5f} {2:10.5f} {3:10.5f}".format(ai.sid, \
                                                            ai.pos[0],\
                                                            ai.pos[1],\
                                                            ai.pos[2])
                    +"  {0:10.5f}  {1:10.5f}  {2:10.5f}".format(ai.vel[0],
                                                                ai.vel[1],
                                                                ai.vel[2])
                    +"\n")
        f.close()

    def make_pair_list(self,rcut=3.0):
        rc2= rcut**2
        h= np.zeros((3,3))
        h[0]= self.a1 *self.alc
        h[1]= self.a2 *self.alc
        h[2]= self.a3 *self.alc
        hi= np.linalg.inv(h)
        # print h
        # print hi
        # lcx= int(1.0/math.sqrt(hi[0,0]**2 +hi[0,1]**2 +hi[0,2]**2)/rcut)
        # lcy= int(1.0/math.sqrt(hi[1,0]**2 +hi[1,1]**2 +hi[1,2]**2)/rcut)
        # lcz= int(1.0/math.sqrt(hi[2,0]**2 +hi[2,1]**2 +hi[2,2]**2)/rcut)
        lcx= int(1.0/math.sqrt(hi[0,0]**2 +hi[1,0]**2 +hi[2,0]**2)/rcut)
        lcy= int(1.0/math.sqrt(hi[0,1]**2 +hi[1,1]**2 +hi[2,1]**2)/rcut)
        lcz= int(1.0/math.sqrt(hi[0,2]**2 +hi[1,2]**2 +hi[2,2]**2)/rcut)
        if lcx == 0: lcx= 1
        if lcy == 0: lcy= 1
        if lcz == 0: lcz= 1
        lcyz= lcy*lcz
        lcxyz= lcx*lcy*lcz
        rcx= 1.0/lcx
        rcy= 1.0/lcy
        rcz= 1.0/lcz
        rcxi= 1.0/rcx
        rcyi= 1.0/rcy
        rczi= 1.0/rcz
        lscl= np.zeros((len(self.atoms),),dtype=int)
        lshd= np.zeros((lcxyz,),dtype=int)
        lscl[:]= -1
        lshd[:]= -1
        # print 'lcx,lcy,lcz,lcxyz=',lcx,lcy,lcz,lcxyz
        # print 'rcx,rcy,rcz=',rcx,rcy,rcz

        #...make a linked-cell list
        for i in range(len(self.atoms)):
            pi= self.atoms[i].pos
            # print pi
            #...assign a vector cell index
            mx= int(pi[0]*rcxi)
            my= int(pi[1]*rcyi)
            mz= int(pi[2]*rczi)
            m= mx*lcyz +my*lcz +mz
            # print i,pi,mx,my,mz,m
            lscl[i]= lshd[m]
            lshd[m]= i

        #...make a pair list
        self.nlspr= np.zeros((self.num_atoms(),),dtype=int)
        self.lspr= np.zeros((self.num_atoms(),_maxnn),dtype=int)
        self.lspr[:]= -1
        # self.lspr= []
        # for i in range(len(self.atoms)):
        #     self.lspr.append([])
            
        for ia in range(len(self.atoms)):
            # if ia % 10000 == 0:
            #     print 'ia=',ia
            ai= self.atoms[ia]
            pi= ai.pos
            mx= int(pi[0]*rcxi)
            my= int(pi[1]*rcyi)
            mz= int(pi[2]*rczi)
            m= mx*lcyz +my*lcz +mz
            #print 'ia,pi,mx,my,mz,m=',ia,pi[0:3],mx,my,mz,m
            for kuz in range(-1,2):
                m1z= mz +kuz
                if m1z < 0: m1z += lcz
                if m1z >= lcz: m1z -= lcz
                for kuy in range(-1,2):
                    m1y= my +kuy
                    if m1y < 0: m1y += lcy
                    if m1y >= lcy: m1y -= lcy
                    for kux in range(-1,2):
                        m1x= mx +kux
                        if m1x < 0: m1x += lcx
                        if m1x >= lcx: m1x -= lcx
                        m1= m1x*lcyz +m1y*lcz +m1z
                        ja= lshd[m1]
                        if ja== -1: continue
                        self.scan_j_in_cell(ia,pi,ja,lscl,h,rc2)
        #...after makeing lspr
        # for ia in range(len(self.atoms)):
        #     print ia,self.lspr[ia]

    def scan_j_in_cell(self,ia,pi,ja,lscl,h,rc2):
        if ja == ia: ja = lscl[ja]
        if ja == -1: return 0
        if not ja in self.lspr[ia]:
            aj= self.atoms[ja]
            pj= aj.pos
            xij= pj-pi
            xij= xij -np.round(xij)
            rij= np.dot(h,xij)
            rij2= rij[0]**2 +rij[1]**2 +rij[2]**2
            if rij2 < rc2:
                n= self.nlspr[ia]
                self.lspr[ia,n]= ja
                self.nlspr[ia] += 1
                if self.nlspr[ia] >= _maxnn:
                    print ' [Error] self.nlspr[{0}] >= _maxnn !!!'.format(ia)
                    print self.nlspr[ia]
                    print self.lspr[ia]
                    sys.exit()
        ja= lscl[ja]
        self.scan_j_in_cell(ia,pi,ja,lscl,h,rc2)

    def _pbc(self,x):
        if x < 0.:
            return x +1.0
        elif x >= 1.0:
            return x -1.0
        else:
            return x

    def assign_pbc(self):
        for ai in self.atoms:
            ai.pos[0]= self._pbc(ai.pos[0])
            ai.pos[1]= self._pbc(ai.pos[1])
            ai.pos[2]= self._pbc(ai.pos[2])

    def get_expansion_num(self,length):
        """
        Compute expansion digits along a1, a2, and a3 from a given *length*.
        System size should be larger than the *length*.
        """
        h= np.zeros((3,3),dtype=np.float)
        h[0]= self.a1 *self.alc
        h[1]= self.a2 *self.alc
        h[2]= self.a3 *self.alc
        vol= np.abs(np.dot(h[0],np.cross(h[1],h[2])))
        l1= np.abs(vol/(np.linalg.norm(np.cross(h[1],h[2]))))
        l2= np.abs(vol/(np.linalg.norm(np.cross(h[2],h[0]))))
        l3= np.abs(vol/(np.linalg.norm(np.cross(h[0],h[1]))))
        n1= int(np.ceil(length/l1))
        n2= int(np.ceil(length/l2))
        n3= int(np.ceil(length/l3))
        return n1,n2,n3

    def expand(self,n1,n2,n3):
        #...expand unit vectors
        self.a1= self.a1*n1
        self.a2= self.a2*n2
        self.a3= self.a3*n3
        n123= n1*n2*n3
        nsid= 0
        for ai in self.atoms:
            nsid= max(nsid,ai.sid)
        natm0= self.num_atoms()
        atoms0= copy.copy(self.atoms)
        self.atoms= []
        aid= 0
        for i1 in range(n1):
            for i2 in range(n2):
                for i3 in range(n3):
                    for ai0 in atoms0:
                        aid += 1
                        ai= Atom()
                        ai.sid= ai0.sid
                        ai.ifmv= ai0.ifmv
                        x= ai0.pos[0]/n1 +1.0/n1*i1
                        y= ai0.pos[1]/n2 +1.0/n2*i2
                        z= ai0.pos[2]/n3 +1.0/n3*i3
                        ai.set_pos(x,y,z)
                        ai.set_vel(ai0.vel[0],ai0.vel[1],ai0.vel[2])
                        ai.set_id(aid)
                        self.atoms.append(ai)

    @classmethod
    def parse_filename(self,filename):
        for format in _file_formats:
            if re.search(format,filename):
                return format
        
if __name__ == "__main__":

    args= docopt(__doc__)

    infmt= args['--in-format']
    outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']

    psys= PMDSystem(fname=infname,ffmt=infmt)

    if not outfmt == None:
        outfmt= psys.parse_filename(outfname)
    
    if outfmt == 'pmd':
        psys.write_pmd(outfname)
    elif outfmt == 'akr':
        psys.write_akr(outfname)
    elif outfmt == 'POSCAR':
        psys.write_POSCAR(outfname)
