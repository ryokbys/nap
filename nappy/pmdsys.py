#!/usr/bin/env python
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
  --specorder=SPECORDER
              Order of species. [default: Al,Mg,Si]
  --scale=SCALE
              Scale the cell. [default: None]
"""

import math
import sys,copy,re
from datetime import datetime

import numpy as np
from docopt import docopt

from atom import Atom,get_symbol_from_number,get_number_from_symbol

#...constants
_maxnn= 100
_file_formats= ('pmd',
                'smd',
                'akr',
                'POSCAR',
                'dump',
                'xsf')

class PMDSystem(object):
    """
    Contains cell information and atoms, and provides some functionalities.
    """

    def __init__(self,fname=None,ffmt=None,specorder=[]):
        self.alc= 1.0
        self.a1= np.zeros(3)
        self.a2= np.zeros(3)
        self.a3= np.zeros(3)
        self.atoms= []
        self.specorder= specorder

        if not fname == None:
            if ffmt == None or \
               ffmt not in _file_formats:
                ftype= parse_filename(fname)
            else:
                ftype = ffmt
            if ftype == 'pmd':
                self.read_pmd(fname)
            elif ftype == 'smd':
                self.read_pmd(fname)
            elif ftype == 'akr':
                self.read_akr(fname)
            elif ftype == 'POSCAR':
                self.read_POSCAR(fname)
            elif ftype == 'dump':
                self.read_dump(fname)
            elif ftype == 'xsf':
                self.read_xsf(fname)
            else:
                raise IOError('Cannot detect input file format.')

    def set_lattice(self,alc,a1,a2,a3):
        self.alc= alc
        self.a1[:]= a1[:]
        self.a2[:]= a2[:]
        self.a3[:]= a3[:]


    def get_hmat(self):
        """
        Be careful about the definition of H-matrix here.
        It is transpose of the definition of that in POSCAR.
        Here hmat = [a1,a2,a3] so that
        pos = hmat * spos
        where pos and spos are real Cartessian position and scaled position.
        """
        hmat = np.zeros((3,3),dtype=float)
        hmat[:,0] = self.a1 *self.alc
        hmat[:,1] = self.a2 *self.alc
        hmat[:,2] = self.a3 *self.alc
        return hmat

    def set_hmat(self,hmat):
        self.alc = 1.0
        self.a1[:] = hmat[:,0]
        self.a2[:] = hmat[:,1]
        self.a3[:] = hmat[:,2]

    def get_lattice_vectors(self):
        return self.a1*self.alc, self.a2*self.alc, self.a3*self.alc

    def get_lattice_lengths(self):
        a = np.linalg.norm(self.a1*self.alc)
        b = np.linalg.norm(self.a2*self.alc)
        c = np.linalg.norm(self.a3*self.alc)
        return a,b,c

    def get_lattice_angles(self):
        a = np.linalg.norm(self.a1)
        b = np.linalg.norm(self.a2)
        c = np.linalg.norm(self.a3)
        bc = np.cross(self.a2,self.a3)
        ac = np.cross(self.a1,self.a3)
        ab = np.cross(self.a1,self.a2)
        bc = np.linalg.norm(bc)
        ac = np.linalg.norm(ac)
        ab = np.linalg.norm(ab)
        # make it inside the range of arcsin
        ta = min(1.0-1.0e-10,bc/b/c)
        tb = min(1.0-1.0e-10,ac/a/c)
        tc = min(1.0-1.0e-10,ab/a/b)
        alpha = np.arcsin(ta)
        beta  = np.arcsin(tb)
        gamma = np.arcsin(tc)
        return alpha,beta,gamma

    def add_atom(self,atom):
        if self.specorder and atom.symbol:
            atom.set_sid(self.specorder.index(atom.symbol)+1)
        self.atoms.append(atom)

    def remove_atom(self,ia):
        """
        Remove an atom of ID ia from the atoms list.
        """
        self.atoms.pop(ia)

    def reset_ids(self):
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)

    def num_atoms(self,sid=0):
        if sid == 0:
            return len(self.atoms)
        else:
            n= 0
            for ai in self.atoms:
                if ai.sid == sid:
                    n += 1
            return n


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

    def volume(self):
        return self.alc**3 *np.abs(np.dot(self.a1,np.cross(self.a2,self.a3)))

    def get_real_positions(self):
        hmat = self.get_hmat()
        spos = self.get_scaled_positions()
        pos = np.zeros((self.num_atoms(),3),dtype=float)
        for ia,a in enumerate(self.atoms):
            pos[ia,:] = np.dot(hmat,spos[ia])
        return pos

    def get_scaled_positions(self):
        spos = np.zeros((self.num_atoms(),3),dtype=float)
        for ia,a in enumerate(self.atoms):
            spos[ia,:] = a.pos
        return spos


    def get_symbols(self):
        symbols = []
        for a in self.atoms:
            symbols.append(a.symbol)
        return symbols

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
        elif ftype == "dump":
            self.read_dump(fname)
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
        symbol = None
        for i in range(natm):
            data= [float(x) for x in f.readline().split()]
            ai= Atom()
            ai.decode_tag(data[0])
            if self.specorder:
                symbol = self.specorder[ai.sid-1]
            if symbol and ai.symbol != symbol:
                ai.set_symbol(symbol)
            ai.set_pos(data[1],data[2],data[3]) # position
            ai.set_vel(data[4],data[5],data[6]) # velocity
            ai.set_ekin(data[7])
            ai.set_epot(data[8])
            ai.set_strs(data[9],data[10],data[11],
                        data[12],data[13],data[14])
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
                    +"  {0:8.4f}  {1:8.4f}  {2:8.4f}".format(ai.vel[0], 
                                                             ai.vel[1],
                                                             ai.vel[2])
                    +"  {0:4.1f}  {1:4.1f}".format(0.0, 0.0)
                    +"  {0:4.1f}  {1:4.1f}  {2:4.1f}".format(0.0, 0.0, 0.0)
                    +"  {0:4.1f}  {1:4.1f}  {2:4.1f}".format(0.0, 0.0, 0.0)
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
                spcs = copy.deepcopy(buff)
                buff= f.readline().split()
                if not self.specorder:
                    self.specorder = spcs
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
                sindex=0
                symbol = None
                for n in num_species:
                    m += n
                    if i < m:
                        if spcs and self.specorder:
                            sid = self.specorder.index(spcs[sindex]) + 1
                            symbol = spcs[sindex]
                        break
                    sid += 1
                    sindex += 1
                ai.set_id(i+1)
                ai.set_sid(sid)
                if symbol:
                    ai.symbol = symbol
                ai.set_pos(float(buff[0]),float(buff[1]),float(buff[2]))
                ai.set_vel(0.0,0.0,0.0)
                self.atoms.append(ai)
                
            

    def write_POSCAR(self,fname='POSCAR'):
        f=open(fname,'w')
        f.write('Generated by pmdsys.py at {0}.'.format(datetime.now().strftime('%Y-%m-%d')))
        if self.specorder:
            for s in self.specorder:
                f.write(' '+s)
        f.write('\n')
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
        # if specorder is defined, write 
        if self.specorder:
            for i,ns in enumerate(num_species):
                s = self.specorder[i]
                f.write(' {0:>4s}'.format(s))
            f.write('\n')
        for n in num_species:
            f.write(' {0:4d}'.format(n))
        f.write('\n')
        # comments
        f.write('Selective dynamics\n')
        f.write('Direct\n')
        # atom positions
        # before writing out, the order should be sorted
        outorder = []
        for spc in self.specorder:
            for ia,ai in enumerate(self.atoms):
                if ai.symbol == spc:
                    outorder.append(ia)
        if len(outorder) != len(self.atoms):
            raise RuntimeError(' len(outorder) != len(self.atoms)')
        for ia in outorder:
            ai = self.atoms[ia]
            f.write(' {0:15.7f} {1:15.7f} {2:15.7f} T T T\n'.format(
                ai.pos[0],ai.pos[1],ai.pos[2]))
        f.close()


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
        symbol = None
        for i in range(natm):
            data= [float(x) for x in f.readline().split()]
            ai= Atom()
            ai.set_sid(data[0])
            ai.set_pos(data[1],data[2],data[3])
            ai.set_vel(data[4],data[5],data[6])
            if self.specorder:
                symbol = self.specorder[ai.sid-1]
            if symbol and ai.symbol != symbol:
                ai.set_symbol(symbol)
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

    def read_dump(self,fname="dump"):
        f=open(fname,'r')
        mode= 'None'
        ixyz= 0
        iatm= 0
        symbol = None
        self.atoms= []
        for line in f.readlines():
            if 'ITEM: NUMBER OF ATOMS' in line:
                mode= 'NUMBER OF ATOMS'
                continue
            elif 'ITEM: BOX BOUNDS' in line:
                mode= 'BOX BOUNDS'
                continue
            elif 'ITEM: ATOMS' in line:
                mode= 'ATOMS'
                continue
            
            if mode == 'NUMBER OF ATOMS':
                natm= int(line.split()[0])
            elif mode == 'BOX BOUNDS':
                if ixyz == 0:
                    xlo= float(line.split()[0])
                    xhi= float(line.split()[1])
                    xlen= xhi -xlo
                elif ixyz == 1:
                    ylo= float(line.split()[0])
                    yhi= float(line.split()[1])
                    ylen= yhi -ylo
                elif ixyz == 2:
                    zlo= float(line.split()[0])
                    zhi= float(line.split()[1])
                    zlen= zhi -zlo
                ixyz += 1
            elif mode == 'ATOMS':
                if iatm < natm:
                    data= line.split()
                    ai= Atom()
                    ai.set_sid(int(data[1]))
                    if self.specorder:
                        symbol = self.specorder[ai.sid-1]
                    if symbol and ai.symbol != symbol:
                        ai.set_symbol(symbol)
                    xi= float(data[2])
                    yi= float(data[3])
                    zi= float(data[4])
                    xi= xi-xlo -int((xi-xlo)/xlen)*xlen
                    yi= yi-ylo -int((yi-ylo)/ylen)*ylen
                    zi= zi-zlo -int((zi-zlo)/zlen)*zlen
                    if xi < 0.0:
                        xi= xi +xlen
                    if yi < 0.0:
                        yi= yi +ylen
                    if zi < 0.0:
                        zi= zi +zlen
                    ai.set_pos(xi/xlen,yi/ylen,zi/zlen)
                    ai.set_vel(0.0,0.0,0.0)
                    self.atoms.append(ai)
                iatm += 1
        self.alc= 1.0
        self.a1[0]= xlen
        self.a2[1]= ylen
        self.a3[2]= zlen
        # print self.alc
        # print self.a1[:]
        # print self.a2[:]
        # print self.a3[:]
        f.close()

    def write_dump(self,fname='dump'):
        """
        Write LAMMPS dump format file.
        Only applicable to orthogonal system.
        """
        f= open(fname,'w')
        f.write("ITEM: TIMESTEP\n")
        f.write("0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write("{0:d}\n".format(len(self.atoms)))
        f.write("ITEM: BOX BOUNDS xy xz yz\n")
        a,b,c = hmat_to_lammps(self.get_hmat())
        xlo = ylo = zlo = 0.0
        xhi = a[0]
        xy  = b[0]
        yhi = b[1]
        xz  = c[0]
        yz  = c[1]
        zhi = c[2]
        xlo_bound = xlo +min(0.0, xy, xz, xy+xz)
        xhi_bound = xhi +max(0.0, xy, xz, xy+xz)
        ylo_bound = ylo +min(0.0, yz)
        yhi_bound = yhi +max(0.0, yz)
        zlo_bound = zlo
        zhi_bound = zhi
        # f.write("{0:15.4f}  {1:15.4f}\n".format(0.0, self.a1[0]))
        # f.write("{0:15.4f}  {1:15.4f}\n".format(0.0, self.a2[1]))
        # f.write("{0:15.4f}  {1:15.4f}\n".format(0.0, self.a3[2]))
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(xlo_bound,
                                                         xhi_bound,
                                                         xy))
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(ylo_bound,
                                                         yhi_bound,
                                                         xz))
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(zlo_bound,
                                                         zhi_bound,
                                                         yz))
        f.write("ITEM: ATOMS id type x y z vx vy vz"
                +" ekin epot sxx syy szz syz sxz sxy\n")
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            x= ai.pos[0] *self.a1[0]
            y= ai.pos[1] *self.a2[1]
            z= ai.pos[2] *self.a3[2]
            vx= ai.vel[0]
            vy= ai.vel[1]
            vz= ai.vel[2]
            ekin= ai.ekin
            epot= ai.epot
            sti= ai.strs
            f.write("{0:8d} {1:3d} ".format(i+1,ai.sid))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(x,y,z))
            f.write("{0:8.3f} {1:8.3f} {2:8.3f} ".format(vx,vy,vz))
            f.write("{0:8.3f} {1:8.3f} ".format(ekin,epot))
            f.write("{0:8.3f} {1:8.3f} {2:8.3f} ".format(sti[0],
                                                         sti[1],
                                                         sti[2]))
            f.write("{0:8.3f} {1:8.3f} {2:8.3f} ".format(sti[3],
                                                         sti[4],
                                                         sti[5]))
            f.write("\n")
        f.close()

    def read_xsf(self,fname="xsf"):
        f=open(fname,'r')
        mode= 'None'
        ixyz= 0
        iatm= 0
        self.atoms= []
        natm= 0
        for line in f.readlines():
            if 'CRYSTAL' in line:
                mode= 'CRYSTAL'
                continue
            elif 'PRIMVEC' in line:
                mode= 'PRIMVEC'
                continue
            elif 'PRIMCOORD' in line:
                mode= 'PRIMCOORD'
                # Before going further, create inversed h-matrix
                hi = unitvec_to_hi(self.a1,self.a2,self.a3)
                print 'Inversed h-matrix:'
                print hi
                continue
            
            if mode == 'CRYSTAL':
                pass
            elif mode == 'PRIMVEC':
                if ixyz == 0:
                    arr = [ float(x) for x in line.split() ]
                    self.a1[0] = arr[0]
                    self.a1[1] = arr[1]
                    self.a1[2] = arr[2]
                elif ixyz == 1:
                    arr = [ float(x) for x in line.split() ]
                    self.a2[0] = arr[0]
                    self.a2[1] = arr[1]
                    self.a2[2] = arr[2]
                elif ixyz == 2:
                    arr = [ float(x) for x in line.split() ]
                    self.a3[0] = arr[0]
                    self.a3[1] = arr[1]
                    self.a3[2] = arr[2]
                ixyz += 1
            elif mode == 'PRIMCOORD':
                data = line.split()
                if len(data) < 4:
                    natm= int(data[0])
                    continue
                else:
                    if iatm >= natm:
                        continue
                    symbol = get_symbol_from_number(int(data[0]))
                    if symbol not in self.specorder:
                        self.specorder.append(symbol)
                    sid = self.specorder.index(symbol) +1
                    ai= Atom()
                    ai.set_sid(sid)
                    xc= float(data[1])
                    yc= float(data[2])
                    zc= float(data[3])
                    xi,yi,zi = cartessian_to_scaled(hi,xc,yc,zc)
                    ai.set_pos(xi,yi,zi)
                    ai.set_vel(0.0,0.0,0.0)
                    self.atoms.append(ai)
                iatm += 1
        self.alc= 1.0
        print self.alc
        print self.a1[:]
        print self.a2[:]
        print self.a3[:]
        f.close()

    def write_xsf(self,fname='xsf'):
        """
        Write XCrysden xsf format.
        Only applicable to orthogonal system.
        """
        if not self.specorder:
            raise ValueError('Specorder has to be defined to write'
                             +' xsf format file.')
        h = np.zeros((3,3),dtype=float)
        h[0,:] = self.a1[:]
        h[1,:] = self.a2[:]
        h[2,:] = self.a3[:]
        f= open(fname,'w')
        f.write("CRYSTAL\n")
        f.write("PRIMVEC\n")
        f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(self.a1[0],
                                                      self.a1[1],
                                                      self.a1[2]))
        f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(self.a2[0],
                                                      self.a2[1],
                                                      self.a2[2]))
        f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(self.a3[0],
                                                      self.a3[1],
                                                      self.a3[2]))
        f.write("PRIMCOORD\n")
        f.write("{0:>8d}  1\n".format(len(self.atoms)))
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            x,y,z = scaled_to_cartessian(h,ai.pos[0],ai.pos[1],ai.pos[2])
            vx= ai.vel[0]
            vy= ai.vel[1]
            vz= ai.vel[2]
            symbol = self.specorder[ai.sid-1]
            number = get_number_from_symbol(symbol)
            f.write(" {0:3d} ".format(number))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(x,y,z))
            #f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(vx,vy,vz))
            f.write("\n")
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

    def repeat(self,n1,n2,n3):
        #...unit vectors to be repeated
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


    def to_ase_atoms(self):
        """
        Convert PMDSystem object to ASE atoms.
        Note that some information will be abandonned.
        """
        try:
            from ase import Atoms
        except ImportError:
            raise ImportError('ASE Atoms cannot be loaded.')
            
        cell = [self.a1, self.a2, self.a3]
        spos = [ a.pos for a in self.atoms ]
        symbols = [ a.symbol for a in self.atoms ]
        atoms = Atoms(symbols=symbols,
                      cell=cell,
                      scaled_positions=spos,
                      pbc=True)
        return atoms
        
    def from_ase_atoms(self,atoms):
        """
        Convert ASE Atoms object to PMDSystem object.
        """
        self.a1 = np.array(atoms.cell[0])
        self.a2 = np.array(atoms.cell[1])
        self.a3 = np.array(atoms.cell[2])
        spos = atoms.get_scaled_positions()
        symbols = atoms.get_chemical_symbols()
        #...initialize and remake self.specorder
        self.specorder = []
        for s in symbols:
            if s not in self.specorder:
                self.specorder.append(s)
        #...first, initialize atoms array
        self.atoms = []
        #...append each atom from ASE-Atoms
        for ia,spi in enumerate(spos):
            si = symbols[ia]
            ai = Atom()
            sid = self.specorder.index(si)+1
            ai.set_id(ia+1)
            ai.set_sid(sid)
            ai.set_symbol(si)
            ai.set_pos(spi[0],spi[1],spi[2])
            ai.set_vel(0.,0.,0.)
            self.atoms.append(ai)
        return


def parse_filename(filename):
    for fmt in _file_formats:
        if fmt in filename:
            return fmt
    return None

def cartessian_to_scaled(hi,xc,yc,zc):
    """
    Convert an atomic position in Cartessian coordinate
    to scaled position using inversed h-matrix.
    Inversed h-matrix has to be given.
    """
    x = 0.0
    y = 0.0
    z = 0.0
    x += hi[0,0]*xc +hi[0,1]*yc +hi[0,2]*zc
    y += hi[1,0]*xc +hi[1,1]*yc +hi[1,2]*zc
    z += hi[2,0]*xc +hi[2,1]*yc +hi[2,2]*zc
    return x,y,z

def scaled_to_cartessian(h,xs,ys,zs):
    """
    Convert a scaled positions to Cartessian coordinate.
    H-matrix has to be given.
    """
    xc = 0.0
    yc = 0.0
    zc = 0.0
    xc += h[0,0]*xs +h[0,1]*ys +h[0,2]*zs
    yc += h[1,0]*xs +h[1,1]*ys +h[1,2]*zs
    zc += h[2,0]*xs +h[2,1]*ys +h[2,2]*zs
    return xc,yc,zc


def get_axis_and_angle(v,u):
    """
    Get rotation axis and angle between given two vectors v and u.
    """
    lv = np.linalg.norm(v)
    lu = np.linalg.norm(u)
    cs = np.dot(v,u)/lv/lu
    vxu = np.cross(v,u)
    axis = vxu/lv/lu
    sn = np.linalg.norm(vxu)/lv/lu
    ang = np.arccos(cs)
    return axis, ang

def rotate(vector,axis,ang):
    """
    Rotate the given *vector* around the *axis* by *ang*.
    *axis* should be normalized vector.
    """
    rmat = np.zeros((3,3),dtype=float)
    nx,ny,nz = axis[:]
    rmat[0,:] = [ 0., -nz, ny]
    rmat[1,:] = [ nz,  0.,-nx]
    rmat[2,:] = [-ny,  nx, 0.]
    mmat = np.zeros((3,3),dtype=float)
    imat = np.identity(3)
    rmat2 = np.dot(rmat,rmat)
    mmat[:,:] = imat[:,:] +np.sin(ang)*rmat[:,:] \
                +(1.0 -np.cos(ang))*rmat2[:,:]
    return np.dot(mmat,vector)


def hmat_to_lammps(hmat):
    """
    Convert h-matrix to LAMMPS cell vectors.
    LAMMPS cell is defined as,
      a = ( xhi-xlo,       0,       0 )
      b = (      xy, yhi-hlo,       0 )
      c = (      xz,      yz, zhi-zlo )
    """
    a0 = hmat[:,0]
    b0 = hmat[:,1]
    c0 = hmat[:,2]
    #...rotate a0 to x-axis
    xaxis = np.array([1.0, 0.0, 0.0])
    ax0,ang0 = get_axis_and_angle(a0,xaxis)
    a1 = rotate(a0,ax0,ang0)
    b1 = rotate(b0,ax0,ang0)
    c1 = rotate(c0,ax0,ang0)
    #...rotate b1 to xy plane
    b1yz = copy.deepcopy(b1)
    b1yz[0] = 0.0
    yaxis = np.array([0.0, 1.0, 0.0])
    ax1,ang1 = get_axis_and_angle(b1,yaxis)
    a2 = rotate(a1,ax1,ang1)
    b2 = rotate(b1,ax1,ang1)
    c2 = rotate(c1,ax1,ang1)
    return a2,b2,c2

def unitvec_to_hi(a1,a2,a3):
    """
    Convert 3 unitvectors to inversed h-matrix via h-matrix.
    """
    h = np.zeros((3,3),dtype=float)
    for i in range(3):
        h[0,i] = a1[i]
        h[1,i] = a2[i]
        h[2,i] = a3[i]
    return np.linalg.inv(h)


if __name__ == "__main__":

    args= docopt(__doc__)

    infmt= args['--in-format']
    outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']
    specorder= args['--specorder'].split(',')
    scalefactor= args['--scale']

    psys= PMDSystem(fname=infname,ffmt=infmt,specorder=specorder)

    if outfmt == 'None':
        outfmt= parse_filename(outfname)

    if scalefactor != "None":
        psys.alc *= float(scalefactor)

    if outfmt == 'pmd':
        psys.write_pmd(outfname)
    elif outfmt == 'smd':
        psys.write_pmd(outfname)
    elif outfmt == 'akr':
        psys.write_akr(outfname)
    elif outfmt == 'POSCAR':
        psys.write_POSCAR(outfname)
    elif outfmt == 'dump':
        psys.write_dump(outfname)
    elif outfmt == 'xsf':
        psys.write_xsf(outfname)
    else:
        print 'Cannot detect output file format.'
