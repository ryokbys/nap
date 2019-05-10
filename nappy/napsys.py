#!/usr/bin/env python
"""
System information used in pmd, including cell information, lattice constant,
number of atoms, atom species, atom positions.

Users can use `napsys.py` as a converter program as described below.
Available formats are,
  pmd, akr, POSCAR, dump, xsf

Usage:
  napsys.py convert [options] INFILE OUTFILE
  napsys.py analyze [options] INFILE

Options:
  -h, --help  Show this help message and exit.
  --in-format=INFMT
              Format of the input file. [default: None]
  --out-format=OUTFMT
              Format of the output file. [default: None]
  --specorder=SPECORDER
              Order of species. [default: None]
  --scale=SCALE
              Scale the cell. [default: None]
  --shift=SHIFT
              Shift original atom positions within a periodic system. [default: 0.0,0.0,0.0]
  --periodic-copy=COPIES
              Number of copies of a,b,and c directions. [default: 1,1,1]
  --cycle-coord=NCYCLE
              Number of cyclic shift of coordinates, e.g. x-y-z to y-z-x. [default: 0]
  --charges=CHARGES
              Charges of species. Charges should be numbers separated by comma.
              If it is None, no charge is set.
              If it is set and the output format accepts, charges are written.
              [default: None]
"""
from __future__ import print_function
from __future__ import division

import math
import sys, copy
from datetime import datetime

import numpy as np
from docopt import docopt

from atom import Atom, get_symbol_from_number, get_number_from_symbol
from units import kB

#...constants
_maxnn = 100
_file_formats = ('pmd',
                 'akr',
                 'POSCAR',
                 'dump',
                 'xsf',
                 'lammps')


class NAPSystem(object):
    """
    Contains cell information and atoms, and provides some functionalities.
    """

    def __init__(self, fname=None, ffmt=None, specorder=[], ase_atoms=None,
                 charges=[]):
        self.alc = 1.0
        self.a1 = np.zeros(3)
        self.a2 = np.zeros(3)
        self.a3 = np.zeros(3)
        self.atoms = []
        self.specorder = specorder

        specorder_good = False
        for s in self.specorder:
            if len(s) > 0:
                specorder_good = True
                break
        if not specorder_good:
            self.specorder = None
        
        if fname is not None:
            self.read(fname=fname,fmt=ffmt)

        if ase_atoms is not None:
            self.from_ase_atoms(ase_atoms)

        self.set_charges(charges)

    def set_lattice(self, alc, a1, a2, a3):
        """
        Set lattice by specifying lattice constant and 3 vectors.

        Args:
          alc (float): lattice constant
          a1 (float vector of 3 components): a1 vector
          a2 (float vector of 3 components): a2 vector
          a3 (float vector of 3 components): a3 vector
        """
        self.alc = alc
        self.a1[:] = a1[:]
        self.a2[:] = a2[:]
        self.a3[:] = a3[:]
        return None

    def set_lattice_parameters(self,a,b,c,alpha,beta,gamma):
        """
        Set lattice by specifying lattice parameters, a,b,c,alpha,beta,gamma.

        Args:
          a (float): lattice parameter *a*
          b (float): lattice parameter *b*
          c (float): lattice parameter *c*
          alpha (float): angle in degree
          beta (float): angle in degree
          gamma (float): angle in degree
        """
        self.alc = 1.0
        alpha_r = np.radians(alpha)
        beta_r = np.radians(beta)
        gamma_r = np.radians(gamma)
        val = (np.cos(alpha_r) * np.cos(beta_r) - np.cos(gamma_r))\
            / (np.sin(alpha_r) * np.sin(beta_r))
        val = max(abs(val),1.0)
        gamma_star = np.arccos(val)
        self.a1[:] = [float(a), 0.0, 0.0]
        self.a2[:] = [b*np.cos(gamma_r),
                      b*np.sin(gamma_r), 0.0]
        self.a3[:] = [c*np.cos(beta_r),
                      -c*np.sin(beta_r)*np.cos(gamma_star),
                      c*np.sin(beta_r)*np.sin(gamma_star)]
        return None
        
    def get_hmat(self):
        """
        Be careful about the definition of H-matrix here.
        It is transpose of the definition of that in POSCAR.
        Here hmat = [a1,a2,a3] so that
        pos = hmat * spos
        where pos and spos are real Cartesian position and scaled position.
        """
        hmat = np.zeros((3, 3), dtype=float)
        hmat[:, 0] = self.a1 * self.alc
        hmat[:, 1] = self.a2 * self.alc
        hmat[:, 2] = self.a3 * self.alc
        return hmat
    
    def set_hmat(self, hmat):
        self.alc = 1.0
        self.a1[:] = hmat[:, 0]
        self.a2[:] = hmat[:, 1]
        self.a3[:] = hmat[:, 2]

    def get_hmat_inv(self):
        hmat = self.get_hmat()
        return np.linalg.inv(hmat)

    def set_specorder(self,*specorder):
        #...Check if the number of species is relevant
        maxsid = 0
        for ai in self.atoms:
            maxsid = max(maxsid,ai.sid)
        if len(specorder) < maxsid:
            txt = 'Number of species is not sufficient.' \
                  +'Must be greater than {0:d}'.format(maxsid)
            raise ValueError(txt)
        self.specorder = specorder
        for ai in self.atoms:
            ai.symbol = specorder[ai.sid-1]
            # try:
            #     ai.sid = self.specorder.index(ai.symbol) +1
            # except ValueError:
            #     ai.sid = 1
                
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

    def get_reciprocal_vectors(self):
        a1,a2,a3 = self.get_lattice_vectors()
        v = np.dot(a1,np.cross(a2,a3))
        b1 = 2.0 *np.pi *np.cross(a2,a3)/v
        b2 = 2.0 *np.pi *np.cross(a3,a1)/v
        b3 = 2.0 *np.pi *np.cross(a1,a2)/v
        return b1,b2,b3

    def add_atom(self,atom):
        if self.specorder and atom.symbol:
            if atom.symbol not in self.specorder:
                self.specorder.append(atom.symbol)
            atom.set_sid(self.specorder.index(atom.symbol)+1)
        self.atoms.append(atom)

    def remove_atom(self,ia):
        """
        Remove an atom of ID ia from the atoms list.
        """
        self.atoms.pop(ia)
        return None

    def remove_atoms(self,dellist):
        """
        Remove atoms of given list of IDs from the system.
        """
        dellist.sort(reverse=True)
        for ia in dellist:
            if ia >= len(self.atoms):
                continue
            self.atoms.pop(ia)
        return None

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
        try:
            if len(self.specorder) > 0:
                for i in range(len(self.specorder)):
                    num_species.append(0)
                for ai in self.atoms:
                    num_species[self.specorder.index(ai.symbol)] += 1
            else:
                for ai in self.atoms:
                    max_nsp= max(max_nsp,ai.sid)
                for i in range(max_nsp):
                    num_species.append(0)
                for ai in self.atoms:
                    num_species[ai.sid-1] += 1
        except:
            num_species.append(len(self.atoms))
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

    def set_symbols(self,symbols):
        """
        Set symbols of all atoms and append new symbols
        if they are not in the current specorder.
        """
        if len(symbols) != len(self.atoms):
            raise RuntimeError('len(symbols) != len(self.atoms)')
        for i,a in enumerate(self.atoms):
            a.set_symbol(symbols[i])
        for s in symbols:
            if s not in self.specorder:
                self.specorder.append(s)
        #...Reset sids
        for i,a in enumerate(self.atoms):
            s = a.symbol
            a.set_sid(self.specorder.index(s)+1)
        return None

    def get_charges(self):
        return self.charges

    def set_charges(self,charges):
        self.charges = charges
        if len(self.charges) > 0:
            if self.specorder is not None \
               and len(self.charges) < len(self.specorder):
                lenc = len(self.charges)
                for i in range(len(self.specorder)-lenc):
                    self.charges.append(0.0)
        return None
        
        
    def write(self,fname="pmdini",fmt=None):
        if fmt in (None,'None'):
            fmt= parse_filename(fname)

        if fmt == 'pmd':
            self.write_pmd(fname)
        elif fmt == 'akr':
            self.write_akr(fname)
        elif fmt == 'POSCAR':
            self.write_POSCAR(fname)
        elif fmt == 'dump':
            self.write_dump(fname)
        elif fmt == 'xsf':
            self.write_xsf(fname)
        elif fmt == 'lammps':
            if len(self.charges) > 0:
                self.write_lammps_data(fname,atom_style='charge')
            else:
                self.write_lammps_data(fname)
        else:
            raise ValueError('Cannot detect output file format: '+fmt)

    def read(self,fname="pmdini",fmt=None):
        if fmt in (None, 'None'):
            fmt= parse_filename(fname)
        
        if fmt == 'pmd':
            self.read_pmd(fname)
        elif fmt == 'akr':
            self.read_akr(fname)
        elif fmt == 'POSCAR':
            self.read_POSCAR(fname)
        elif fmt == 'dump':
            self.read_dump(fname)
        elif fmt == 'xsf':
            self.read_xsf(fname)
        elif fmt == 'lammps':
            self.read_lammps_data(fname)
        else:
            raise IOError('Cannot detect input file format: '+fmt)

    def read_pmd(self,fname='pmdini'):
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

    def write_pmd(self,fname='pmdini'):
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
                else:
                    for s in spcs:
                        if s not in self.specorder:
                            self.specorder.append(s)
            num_species= np.array([ int(n) for n in buff])
            try:
                spcs
            except NameError:
                spcs = self.specorder
            #...Check number of species in POSCAR file and in specorder
            if len(num_species) > len(self.specorder):
                msg = '''
Numbers of species in POSCAR is greater than the one in specorder, which should be the same or less.
Number of species in POSCAR = {0:d}
You need to specify the species order correctly with --specorder option.
                '''.format(len(num_species))
                raise ValueError(msg)
            natm= 0
            for n in num_species:
                natm += n
            #print("Number of atoms = {0:5d}".format(natm))
            # 7th or 8th line: comment
            c7= f.readline()
            if c7[0] in ('s','S'):
                c7= f.readline()
            if c7[0] in ('c','C'): # positions are in Cartesian coordinate
                hi = unitvec_to_hi(self.a1,self.a2,self.a3)
                coord = 'cartesian'
            else:
                coord = 'scaled'
            
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
                pos = [ float(buff[0]), float(buff[1]), float(buff[2])]
                if coord == 'cartesian':
                    x1,x2,x3 = cartesian_to_scaled(hi,pos[0],pos[1],pos[2])
                elif coord == 'scaled':
                    x1,x2,x3 = pos[0],pos[1],pos[2]
                ai.set_pos(x1,x2,x3)
                ai.set_vel(0.0,0.0,0.0)
                self.atoms.append(ai)
                
    def write_POSCAR(self,fname='POSCAR'):
        f=open(fname,'w')
        f.write('Generated by napsys.py at {0}.'.format(datetime.now().strftime('%Y-%m-%d')))
        if self.specorder:
            for s in self.specorder:
                f.write(' '+s)
        f.write('\n')
        # lattice vector
        f.write(' {0:15.7f}\n'.format(self.alc))
        # cell vectors
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a1[0],
                                                          self.a1[1],
                                                          self.a1[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a2[0],
                                                          self.a2[1],
                                                          self.a2[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a3[0],
                                                          self.a3[1],
                                                          self.a3[2]))
        # count num of atoms per specie
        num_species= self.num_species()
        # if specorder is defined, write species names
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
            # print 'spc =',spc
            for ia,ai in enumerate(self.atoms):
                # print ia,ai.symbol,ai.sid
                if ai.symbol == spc:
                    outorder.append(ia)
        if len(outorder) != len(self.atoms):
            print('len(outorder),len(self.atoms)=', len(outorder),len(self.atoms))
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
        self.alc= 1.0
        xy = 0.0
        xz = 0.0
        yz = 0.0
        aux_exists = {
            'x': -1, 'y': -1, 'z': -1,
            'xu': -1, 'yu': -1, 'zu': -1,
            'fx': -1, 'fy': -1, 'fz': -1,
            'ekin': -1,
            'epot': -1,
            'sxx': -1,
            'syy': -1,
            'szz': -1,
            'syz': -1,
            'sxz': -1,
            'sxy': -1,
            'chg': -1,
            'chi': -1
        }
        for line in f.readlines():
            if 'ITEM' in line:
                if 'NUMBER OF ATOMS' in line:
                    mode= 'NUMBER OF ATOMS'
                    continue
                elif 'BOX BOUNDS' in line:
                    mode= 'BOX BOUNDS'
                    continue
                elif 'ATOMS' in line:
                    mode= 'ATOMS'
                    aux_names = [ name for i,name in enumerate(line.split()) if i > 1 ]
                    for ia,a in enumerate(aux_names):
                        if a in aux_exists.keys():
                            aux_exists[a] = ia
                    continue
                elif 'TIMESTEP' in line:
                    mode= 'TIMESTEP'
                    continue
                
            if mode == 'TIMESTEP':
                timestep = int(line.split()[0])
            elif mode == 'NUMBER OF ATOMS':
                natm= int(line.split()[0])
            elif mode == 'BOX BOUNDS':
                data = line.split()
                if ixyz == 0:
                    xlo_bound= float(data[0])
                    xhi_bound= float(data[1])
                    if len(data) > 2:
                        xy = float(data[2])
                elif ixyz == 1:
                    ylo_bound= float(line.split()[0])
                    yhi_bound= float(line.split()[1])
                    if len(data) > 2:
                        xz = float(data[2])
                elif ixyz == 2:
                    zlo_bound= float(line.split()[0])
                    zhi_bound= float(line.split()[1])
                    if len(data) > 2:
                        yz = float(data[2])
                ixyz += 1
                if ixyz > 2:
                    xlo = xlo_bound -min(0.0,xy,xz,xy+yz)
                    xhi = xhi_bound -max(0.0,xy,xz,xy+yz)
                    ylo = ylo_bound -min(0.0,yz)
                    yhi = yhi_bound -max(0.0,yz)
                    zlo = zlo_bound
                    zhi = zhi_bound
                    self.a1 = np.array([xhi-xlo, 0., 0.],dtype=float)
                    self.a2 = np.array([xy, yhi-ylo, 0.],dtype=float)
                    self.a3 = np.array([xz, yz, zhi-zlo],dtype=float)
                    hmat = self.get_hmat()
                    hmati= np.linalg.inv(hmat)
            elif mode == 'ATOMS':
                if iatm < natm:
                    data= line.split()
                    ai= Atom()
                    if data[1].isdigit():
                        ai.set_sid(int(data[1]))
                    else:
                        symbol = data[1]
                        if symbol not in self.specorder:
                            ValueError('Symbol {0:s} is not in specorder.'.format(symbol))
                        ai.set_sid(self.specorder.index(symbol)+1)
                    if self.specorder:
                        symbol = self.specorder[ai.sid-1]
                    if symbol and ai.symbol != symbol:
                        ai.set_symbol(symbol)
                    ix = max( aux_exists['x'], aux_exists['xu'])
                    iy = max( aux_exists['y'], aux_exists['yu'])
                    iz = max( aux_exists['z'], aux_exists['zu'])
                    if ix < 0 or iy < 0 or iz < 0:
                        raise ValueError('Not enough coordinate info.\nPlease check the dump file format.')
                    x0= float(data[ix])
                    y0= float(data[iy])
                    z0= float(data[iz])
                    x = hmati[0,0]*x0 +hmati[0,1]*y0 +hmati[0,2]*z0
                    y = hmati[1,0]*x0 +hmati[1,1]*y0 +hmati[1,2]*z0
                    z = hmati[2,0]*x0 +hmati[2,1]*y0 +hmati[2,2]*z0
                    x = self._pbc(x)
                    y = self._pbc(y)
                    z = self._pbc(z)
                    ai.set_pos(x,y,z)
                    ai.set_vel(0.0,0.0,0.0)
                    for k,v in aux_exists.items():
                        if v >= 0:
                            if k == 'ekin':
                                eki = float(data[v])
                                ai.set_ekin(float(data[v]))
                                ai.set_aux(key='temperature',
                                           value=eki*2.0/3/kB)
                            elif k == 'epot':
                                ai.set_epot(float(data[v]))
                            elif k == 'chg':
                                ai.set_aux(key=k,value=float(data[v]))
                            elif k == 'chi':
                                ai.set_aux(key=k,value=float(data[v]))
                            elif k == 'sxx':
                                ai.strs[0] = float(data[v])
                            elif k == 'syy':
                                ai.strs[1] = float(data[v])
                            elif k == 'szz':
                                ai.strs[2] = float(data[v])
                            elif k == 'syz' or k == 'szy':
                                ai.strs[3] = float(data[v])
                            elif k == 'sxz' or k == 'szx':
                                ai.strs[4] = float(data[v])
                            elif k == 'sxy' or k == 'syx':
                                ai.strs[5] = float(data[v])
                    self.atoms.append(ai)
                iatm += 1
        # print self.alc
        # print self.a1[:]
        # print self.a2[:]
        # print self.a3[:]
        f.close()

    def write_dump(self,fname='dump'):
        """
        Write LAMMPS dump format file.
        """
        f= open(fname,'w')
        f.write("ITEM: TIMESTEP\n")
        f.write("0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write("{0:d}\n".format(len(self.atoms)))

        hmat = self.get_hmat()
        poss = np.zeros((len(self.atoms),3),dtype=float)
        for i in range(len(self.atoms)):
            poss[i,:] = self.atoms[i].pos[:]
        xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,poss = to_lammps(hmat,poss)
        xlo_bound = xlo +min(0.0, xy, xz, xy+xz)
        xhi_bound = xhi +max(0.0, xy, xz, xy+xz)
        ylo_bound = ylo +min(0.0, yz)
        yhi_bound = yhi +max(0.0, yz)
        f.write("ITEM: BOX BOUNDS xy xz yz\n")
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(xlo_bound,
                                                         xhi_bound,
                                                         xy))
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(ylo_bound,
                                                         yhi_bound,
                                                         xz))
        f.write("{0:15.4f} {1:15.4f} {2:15.4f}\n".format(zlo,
                                                         zhi,
                                                         yz))
        
        
        f.write("ITEM: ATOMS id type x y z vx vy vz"
                +" ekin epot sxx syy szz syz sxz sxy\n")
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            # x= ai.pos[0] *self.a1[0] *self.alc
            # y= ai.pos[1] *self.a2[1] *self.alc
            # z= ai.pos[2] *self.a3[2] *self.alc
            pos = poss[i]
            vx= ai.vel[0]
            vy= ai.vel[1]
            vz= ai.vel[2]
            ekin= ai.ekin
            epot= ai.epot
            sti= ai.strs
            # f.write("{0:8d} {1:3d} ".format(i+1,ai.sid))
            f.write("{0:8d} {1:3s} ".format(i+1,ai.symbol))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(pos[0],pos[1],pos[2]))
            f.write("{0:8.3f} {1:8.3f} {2:8.3f} ".format(vx,vy,vz))
            f.write("{0:11.3e} {1:11.3e} ".format(ekin,epot))
            f.write("{0:11.3e} {1:11.3e} {2:11.3e} ".format(sti[0],
                                                            sti[1],
                                                            sti[2]))
            f.write("{0:11.3e} {1:11.3e} {2:11.3e} ".format(sti[3],
                                                            sti[4],
                                                            sti[5]))
            f.write("\n")
        f.close()

    def read_lammps_data(self,fname="data.lammps",atom_style='atomic'):
        f=open(fname,'r')
        mode= 'None'
        iatm= 0
        symbol = None
        self.atoms= []
        self.alc= 1.0
        xy = 0.0
        xz = 0.0
        yz = 0.0
        for line in f.readlines():
            data = line.split()
            if mode == 'None':
                if 'atoms' in line:
                    natm = int(data[0])
                elif 'atom types' in line:
                    nspcs = int(data[0])
                elif 'xlo' in line:
                    xlo = float(data[0])
                    xhi = float(data[1])
                elif 'ylo' in line:
                    ylo = float(data[0])
                    yhi = float(data[1])
                elif 'zlo' in line:
                    zlo = float(data[0])
                    zhi = float(data[1])
                elif 'xy' in line:
                    xy = float(data[0])
                    xz = float(data[1])
                    yz = float(data[2])
                elif 'Atoms' in line:
                    mode = 'Atoms'
                    #...Cell info (xhi,xlo,...) should already be read
                    # self.a1 = np.array([xhi-xlo,xy,xz],dtype=float)
                    # self.a2 = np.array([0.0,yhi-ylo,yz],dtype=float)
                    # self.a3 = np.array([0.0,0.0,zhi-zlo],dtype=float)
                    self.a1 = np.array([xhi-xlo,0.0,0.0],dtype=float)
                    self.a2 = np.array([xy,yhi-ylo,0.0],dtype=float)
                    self.a3 = np.array([xz,yz,zhi-zlo],dtype=float)
                    hmat = self.get_hmat()
                    hmati= np.linalg.inv(hmat)
                    continue
            elif mode == 'Atoms':
                if len(data) >= 5 and iatm < natm:
                    idat = 0
                    ai = Atom()
                    idat += 1
                    ai.set_sid(int(data[idat]))
                    if self.specorder:
                        symbol = self.specorder[ai.sid-1]
                    if symbol and ai.symbol != symbol:
                        ai.set_symbol(symbol)
                    if atom_style == 'charge': 
                        idat += 1
                        chg = float(data[idat])
                        ai.set_aux('charge',chg)
                    idat += 1
                    x0= float(data[idat])
                    idat += 1
                    y0= float(data[idat])
                    idat += 1
                    z0= float(data[idat])
                    x = hmati[0,0]*x0 +hmati[0,1]*y0 +hmati[0,2]*z0
                    y = hmati[1,0]*x0 +hmati[1,1]*y0 +hmati[1,2]*z0
                    z = hmati[2,0]*x0 +hmati[2,1]*y0 +hmati[2,2]*z0
                    x = self._pbc(x)
                    y = self._pbc(y)
                    z = self._pbc(z)
                    ai.set_pos(x,y,z)
                    ai.set_vel(0.0,0.0,0.0)
                    self.atoms.append(ai)
                    iatm += 1
        f.close()

    def write_lammps_data(self,fname='data.lammps',atom_style='atomic'):
        """
        Write LAMMPS data format file.
        The definition of cell vectors is a bit tricky, see the following page
        http://lammps.sandia.gov/doc/Section_howto.html#howto-12
        And also the format of Atoms entry could change depending on 
        `atom_style` which is not given in the same file.
        """
        f= open(fname,'w')
        f.write("LAMMPS data format file written by napsys.py\n")
        f.write("\n")
        f.write("{0:d}  atoms\n".format(len(self.atoms)))
        f.write("{0:d}  atom types\n".format(len(self.num_species())))
        f.write('\n')
        hmat = self.get_hmat()
        poss = np.zeros((len(self.atoms),3),dtype=float)
        for i in range(len(self.atoms)):
            poss[i,:] = self.atoms[i].pos[:]
        xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,poss = to_lammps(hmat,poss)
        f.write("{0:20.10f} {1:20.10f} xlo xhi\n".format(xlo,xhi))
        f.write("{0:20.10f} {1:20.10f} ylo yhi\n".format(ylo,yhi))
        f.write("{0:20.10f} {1:20.10f} zlo zhi\n".format(zlo,zhi))
        if abs(xy) > 1e-8 or abs(xz) > 1e-8 or abs(yz) > 1e-8:
            f.write("{0:20.10f} {1:20.10f} {2:20.10f} xy xz yz\n".format(xy,xz,yz))
        f.write("\n")
        f.write("Atoms\n")
        f.write("\n")
        # pos = np.zeros(3,dtype=float)
        if self.atoms[0].aux.has_key('charge'):
            atom_style = 'charge'
        #print poss
        # poss = spos_to_lammps_pos(hmat,poss)
        #print poss
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            # print hmat
            # print ai.pos
            # pos = np.dot(hmat,ai.pos)
            # print pos
            pos = poss[i]
            f.write("{0:8d} {1:3d} ".format(i+1,ai.sid))
            if atom_style == 'charge':
                f.write('{0:10.4f} '.format(self.charges[ai.sid-1]))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(pos[0],
                                                            pos[1],
                                                            pos[2]))
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
                # print 'Inversed h-matrix:'
                # print hi
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
                if len(data) == 1:
                    natm= int(data[0])
                    continue
                elif len(data) == 2:
                    natm= int(data[0])
                    nspcs= int(data[1])
                    continue
                elif len(data) == 4 or len(data) == 7:
                    if iatm >= natm:
                        continue
                    symbol = get_symbol_from_number(int(data[0]))
                    if symbol not in self.specorder:
                        self.specorder.append(symbol)
                    sid = self.specorder.index(symbol) +1
                    ai= Atom()
                    ai.symbol = symbol
                    ai.set_sid(sid)
                    xc= float(data[1])
                    yc= float(data[2])
                    zc= float(data[3])
                    xi,yi,zi = cartesian_to_scaled(hi,xc,yc,zc)
                    ai.set_pos(xi,yi,zi)
                    ai.set_vel(0.0,0.0,0.0)
                    self.atoms.append(ai)
                    # print 'iatm,symbol,sid,xc,yc,zc = ',iatm,symbol,sid,xc,yc,zc
                else:
                    continue
                iatm += 1
        self.alc= 1.0
        # print self.alc
        # print self.a1[:]
        # print self.a2[:]
        # print self.a3[:]
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
            x,y,z = scaled_to_cartesian(h,ai.pos[0],ai.pos[1],ai.pos[2])
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

    def get_distance(self,ia,ja):
        """
        Compute distance between atoms ia and ja taking the periodic boundary
        condition into account.
        """
        if ia > len(self.atoms):
            raise ValueError('ia > natms, ia,natms = ',ia,len(self.atoms))
        if ja > len(self.atoms):
            raise ValueError('ja > natms, ja,natms = ',ja,len(self.atoms))
        xi = self.atoms[ia].pos
        xj = self.atoms[ja].pos
        xij = xj-xi -np.round(xj-xi)
        hmat = self.get_hmat()
        rij = np.dot(hmat,xij)
        rij2 = rij[0]**2 +rij[1]**2 +rij[2]**2
        return np.sqrt(rij2)

    def get_angle(self,i,j,k):
        """
        Compute angle in degree between bonds i-j and i-k.
        """
        if i > len(self.atoms):
            raise ValueError('i > natms, i,natms = ',i,len(self.atoms))
        if j > len(self.atoms):
            raise ValueError('j > natms, j,natms = ',j,len(self.atoms))
        if k > len(self.atoms):
            raise ValueError('k > natms, k,natms = ',k,len(self.atoms))
        xi = self.atoms[i].pos
        xj = self.atoms[j].pos
        xk = self.atoms[k].pos
        xij = xj-xi -np.round(xj-xi)
        xik = xk-xi -np.round(xk-xi)
        hmat = self.get_hmat()
        rij = np.dot(hmat,xij)
        rik = np.dot(hmat,xik)
        dij = np.linalg.norm(rij)
        dik = np.linalg.norm(rik)
        angle = np.arccos(np.dot(rij,rik)/dij/dik) /np.pi *180.0
        return angle
        
    def make_pair_list(self,rcut=3.0):
        rc2= rcut**2
        h= np.zeros((3,3))
        h[:,0]= self.a1 *self.alc
        h[:,1]= self.a2 *self.alc
        h[:,2]= self.a3 *self.alc
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
        self.assign_pbc()
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
                    print(' [Error] self.nlspr[{0}] >= _maxnn !!!'.format(ia))
                    print(self.nlspr[ia])
                    print(self.lspr[ia])
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

    def shift_atoms(self,sx,sy,sz):
        """
        Shift all the atoms by (sx,sy,sz) in scaled unit.
        """
        for ai in self.atoms:
            newpos = [ ai.pos[0]+sx, ai.pos[1]+sy, ai.pos[2]+sz ]
            ai.set_pos(*newpos)
        return

    def cycle_coord(self,ncycle=0):
        """
        Cyclic shift of coordinates, e.g. x-y-z to y-z-x.
        """
        ncycle = ncycle % 3
        if ncycle == 0:
            return None
        for icycle in range(ncycle):
            a1,a2,a3 = self.get_lattice_vectors()
            spos = self.get_scaled_positions()
            a1 = [ a1[1], a1[2], a1[0] ]
            a2 = [ a2[1], a2[2], a2[0] ]
            a3 = [ a3[1], a3[2], a3[0] ]
            self.set_lattice(1.0, a2, a3, a1)
            for i,ai in enumerate(self.atoms):
                p = spos[i]
                ai.set_pos(p[1],p[2],p[0])
        return None
        
    def repeat(self,n1,n2,n3,n1m=0,n2m=0,n3m=0):
        if n1 == n2 == n3 == 1:
            return None
        #...unit vectors to be repeated
        m1 = n1-n1m
        m2 = n2-n2m
        m3 = n3-n3m
        self.a1= self.a1*m1
        self.a2= self.a2*m2
        self.a3= self.a3*m3
        #n123= m1*m2*m3
        nsid= 0
        for ai in self.atoms:
            nsid= max(nsid,ai.sid)
        #natm0= self.num_atoms()
        atoms0= copy.copy(self.atoms)
        self.atoms= []
        aid= 0
        for i1 in range(n1m,n1):
            for i2 in range(n2m,n2):
                for i3 in range(n3m,n3):
                    for ai0 in atoms0:
                        aid += 1
                        ai= Atom()
                        ai.sid= ai0.sid
                        ai.symbol= ai0.symbol
                        ai.ifmv= ai0.ifmv
                        x= ai0.pos[0]/m1 +1.0/m1*i1
                        y= ai0.pos[1]/m2 +1.0/m2*i2
                        z= ai0.pos[2]/m3 +1.0/m3*i3
                        ai.set_pos(x,y,z)
                        ai.set_vel(ai0.vel[0],ai0.vel[1],ai0.vel[2])
                        ai.set_id(aid)
                        self.atoms.append(ai)

    def add_vacuum(self,va,vb,vc):
        """
        Add vacuum of the given height.
        And atoms are shifted along each direction so that they are placed at the center.
        """
        a,b,c = self.get_lattice_lengths()
        aratio = (a+va)/a
        bratio = (b+vb)/b
        cratio = (c+vc)/c
        self.assign_pbc()
        rpos = self.get_real_positions()
        self.a1 *= aratio
        self.a2 *= bratio
        self.a3 *= cratio
        hmati = self.get_hmat_inv()
        spos = np.zeros((len(self.atoms),3),dtype=float)
        mins = np.array((1.0,1.0,1.0),dtype=float)
        maxs = np.zeros(3,dtype=float)
        for ia in range(len(self.atoms)):
            xi = rpos[ia]
            spos[ia] = np.dot(hmati,xi)
            for l in range(3):
                mins[l] = min(mins[l],spos[ia,l])
                maxs[l] = max(maxs[l],spos[ia,l])
        cntrs = (maxs +mins)/2
        shfts = -(cntrs-0.5)
        for ia in range(len(self.atoms)):
            spos[ia] += shfts
            self.atoms[ia].set_pos(spos[ia,0],spos[ia,1],spos[ia,2])
        return None
        
    def to_ase_atoms(self):
        """
        Convert NAPSystem object to ASE atoms.
        Note that some information will be abandonned.
        """
        try:
            from ase import Atoms
        except ImportError:
            raise ImportError('Cannot load ase Atoms.')

        cell = np.array([self.a1, self.a2, self.a3])
        cell *= self.alc
        a = np.linalg.norm(cell[0,:])
        b = np.linalg.norm(cell[1,:])
        c = np.linalg.norm(cell[2,:])
        spos = [ at.pos for at in self.atoms ]
        vels = [ [at.vel[0]*a, at.vel[1]*b, at.vel[2]*c]  for at in self.atoms ]
        symbols = [ at.symbol for at in self.atoms ]
        atoms = Atoms(symbols=symbols,
                      cell=cell,
                      scaled_positions=spos,
                      pbc=True)
        atoms.set_velocities(vels)
        return atoms

    @classmethod
    def from_ase_atoms(cls,ase_atoms,specorder=None):
        """
        Convert ASE Atoms object to NAPSystem object.
        """
        spcorder = []
        if specorder is not None:
            spcorder = specorder
        symbols = ase_atoms.get_chemical_symbols()
        spos = ase_atoms.get_scaled_positions()
        #...initialize and remake self.specorder
        for s in symbols:
            if s not in spcorder:
                spcorder.append(s)
        nap = cls(specorder=spcorder)
        nap.alc= 1.0
        nap.a1[:] = ase_atoms.cell[0]
        nap.a2[:] = ase_atoms.cell[1]
        nap.a3[:] = ase_atoms.cell[2]
        #...first, initialize atoms array
        nap.atoms = []
        #...append each atom from ASE-Atoms
        for ia,spi in enumerate(spos):
            si = symbols[ia]
            ai = Atom()
            sid = nap.specorder.index(si)+1
            ai.set_id(ia+1)
            ai.set_sid(sid)
            ai.set_symbol(si)
            ai.set_pos(spi[0],spi[1],spi[2])
            ai.set_vel(0.,0.,0.)
            nap.atoms.append(ai)
        return nap

    def change_unitcell(self,a,b,c):
        """
        Change the current unitcell to the new one with 
        given unit vectors, a, b, and c.
        And the atoms are reduced to those within the new unitcell.
        The new unitcell should be included or the same size as the 
        current unitcell, otherwise there appear vacuum region 
        in the new unitcell.
        """
        #...Repeat the system in order to fill the outer space
        #   arround the current unitcell but in the new unitcell...
        self.repeat(2,2,2,-1,-1,-1)
        rpos = self.get_real_positions()
        #...Set new unitcell 
        self.alc = 1.0
        self.a1 = np.array(a)
        self.a2 = np.array(b)
        self.a3 = np.array(c)
        newhmat = np.zeros((3,3),dtype=float)
        newhmat[:,0] = a[:]
        newhmat[:,1] = b[:]
        newhmat[:,2] = c[:]
        newhmati = np.linalg.inv(newhmat)
        for i in range(len(rpos)):
            spos = np.dot(newhmati,rpos[i])
            self.atoms[i].pos[:]  = spos[:]

        #...Remove atoms outside the unitcell
        remove_ids = []
        tiny = 1.0e-5
        for i in range(len(self.atoms)):
            pi = self.atoms[i].pos
            if pi[0] < 0.0 or pi[0] >= 1.0-tiny or \
               pi[1] < 0.0 or pi[1] >= 1.0-tiny or \
               pi[2] < 0.0 or pi[2] >= 1.0-tiny:
                remove_ids.append(i)
        self.atoms = [ atom for i,atom in enumerate(self.atoms)
                       if i not in remove_ids ]
        return

    def remove_overlapping_atoms(self,criterion=0.01):
        """
        Remove overlapping atoms in the system.
        Atoms with bigger indices in overlapping atoms are removed.
        Judgement of overlap is done by comparing the distance between
        atoms and *criterion*.
        """
        try:
            self.lspr
        except:
            self.make_pair_list(rcut=1.0)
        remove_ids = []
        h = np.zeros((3,3),dtype=float)
        h[:,0] = self.a1 *self.alc
        h[:,1] = self.a2 *self.alc
        h[:,2] = self.a3 *self.alc
        cr2 = criterion*criterion
        for i in range(len(self.atoms)):
            ai = self.atoms[i]
            pi = ai.pos
            for jj in range(self.nlspr[i]):
                j = self.lspr[i,jj]
                aj = self.atoms[j]
                pj = aj.pos
                xij = pj-pi
                xij = xij -np.round(xij)
                rij = np.dot(h,xij)
                rij2= np.linalg.norm(rij)
                if rij2 < cr2:
                    remove_ids.append(max(i,j))
        # print remove_ids
        self.atoms = [ atom for i,atom in enumerate(self.atoms)
                       if i not in remove_ids ]
        return
        
def parse_filename(filename):
    for fmt in _file_formats:
        if fmt in filename:
            return fmt
    return None


def cartesian_to_scaled(hi,xc,yc,zc):
    """
    Convert an atomic position in Cartesian coordinate
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


def scaled_to_cartesian(h,xs,ys,zs):
    """
    Convert a scaled positions to Cartesian coordinate.
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


##def hmat_to_lammps(hmat):
##    """
##    Convert h-matrix to LAMMPS cell vectors.
##    Parameters to be output:
##      xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
##    LAMMPS cell should be defined as,
##      a = ( xhi-xlo,       0,       0 )
##      b = (      xy, yhi-hlo,       0 )
##      c = (      xz,      yz, zhi-zlo )
##    See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
##    """
##    import numpy as np
##    a0 = hmat[:,0]
##    b0 = hmat[:,1]
##    c0 = hmat[:,2]
##    xlo = 0.0
##    ylo = 0.0
##    zlo = 0.0
##    a = np.linalg.norm(a0)
##    b = np.linalg.norm(b0)
##    c = np.linalg.norm(c0)
##    alpha = np.arccos(np.dot(b0,c0)/b/c)
##    beta  = np.arccos(np.dot(a0,c0)/a/c)
##    gamma = np.arccos(np.dot(a0,b0)/a/b)
##    # print 'hmat=',hmat
##    # print 'a,b,c = ',a,b,c
##    # print 'alpha,beta,gamma = ',alpha,beta,gamma
##    xhi = a
##    xy = b*np.cos(gamma)
##    xz = c*np.cos(beta)
##    yhi = np.sqrt(b*b -xy*xy)
##    yz = (b*c*np.cos(alpha) -xy*xz)/yhi
##    zhi = np.sqrt(c*c -xz*xz -yz*yz)
##    # print 'xhi-xlo,yhi-ylo,zhi-zlo= ',xhi-xlo,yhi-ylo,zhi-zlo
##    # print 'xy,     xz,     yz     = ',xy/xhi,xz/xhi,yz/yhi
##
##    if xy > xhi/2:
##        xy -= xhi
##    elif xy < -xhi/2:
##        xy += xhi
##
##    if xz > xhi/2:
##        xz -= xhi
##    elif xz < -xhi/2:
##        xz += xhi
##
##    if yz > yhi/2:
##        yz -= yhi
##    elif yz < -yhi/2:
##        yz += yhi
##    # print 'xy,     xz,     yz     = ',xy/xhi,xz/xhi,yz/yhi
##
##    return xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
##
##def spos_to_lammps_pos(hmat,spos):
##    """
##    Scaled positions in hmat-representation to
##    positions in the lammps basis.
##    """
##    if isinstance(hmat,list):
##        hmat = np.array(hmat)
##        
##    if not isinstance(spos,np.ndarray):
##        if isinstance(spos,list):
##            spos = np.array(spos)
##        else:
##            raise TypeError('spos should be list or numpy.ndarray.')
##    a1 = np.array(hmat[:,0])
##    a2 = np.array(hmat[:,1])
##    a3 = np.array(hmat[:,2])
##    vol = abs(np.dot(a1,np.cross(a2,a3)))
##    a23 = np.cross(a2,a3)
##    a31 = np.cross(a3,a1)
##    a12 = np.cross(a1,a2)
##    amat = np.zeros((3,3),dtype=float)
##    amat[0,:] = a23[:]
##    amat[1,:] = a31[:]
##    amat[2,:] = a12[:]
##    # print 'hmat=',hmat
##    # print 'vol=',vol
##    # print 'a1=',a1
##    # print 'a2=',a2
##    # print 'a3=',a3
##    # print 'amat=',amat
##
##    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz = hmat_to_lammps(hmat)
##    b1 = np.array((xhi-xlo,0.0,0.0))
##    b2 = np.array((xy,yhi-ylo,0.0))
##    b3 = np.array((xz,yz,zhi-zlo))
##    bmat = np.zeros((3,3),dtype=float)
##    bmat[:,0] = b1[:]
##    bmat[:,1] = b2[:]
##    bmat[:,2] = b3[:]
##    # print 'bmat=',bmat
##    if len(spos.shape) == 1:  # only one atom
##        pos = np.zeros(spos.shape,dtype=float)
##        pos = np.dot(hmat,spos)
##        pos = np.dot(bmat,np.dot(amat,pos))/vol
##    elif len(spos.shape) == 2:  # array of atoms
##        pos = np.zeros(spos.shape,dtype=float)
##        for i,sp in enumerate(spos):
##            pos[i] = np.dot(hmat,sp)
##            pos[i] = np.dot(bmat,np.dot(amat,pos[i]))/vol
##    return pos

def to_lammps(hmat,spos):
    """
    Convert h-matrix and scaled positions in napsys to 
    LAMMPS representation.
    Parameters to be output:
      xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, pos
    LAMMPS cell should be defined as,
      a = ( xhi-xlo,       0,       0 )
      b = (      xy, yhi-hlo,       0 )
      c = (      xz,      yz, zhi-zlo )
    See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
    Note that xy, xz, and yz have limitations in their ranges.
    """
    import numpy as np
    if isinstance(hmat,list):
        hmat = np.array(hmat)
        
    if not isinstance(spos,np.ndarray):
        if isinstance(spos,list):
            spos = np.array(spos)
        else:
            raise TypeError('spos should be list or numpy.ndarray.')
    a0 = hmat[:,0]
    b0 = hmat[:,1]
    c0 = hmat[:,2]
    xlo = 0.0
    ylo = 0.0
    zlo = 0.0
    a = np.linalg.norm(a0)
    b = np.linalg.norm(b0)
    c = np.linalg.norm(c0)
    alpha = np.arccos(np.dot(b0,c0)/b/c)
    beta  = np.arccos(np.dot(a0,c0)/a/c)
    gamma = np.arccos(np.dot(a0,b0)/a/b)
    xhi = a
    xy = b*np.cos(gamma)
    xz = c*np.cos(beta)
    yhi = np.sqrt(b*b -xy*xy)
    yz = (b*c*np.cos(alpha) -xy*xz)/yhi
    zhi = np.sqrt(c*c -xz*xz -yz*yz)
    x = xhi-xlo
    y = yhi-ylo
    z = zhi-zlo
    
    lxy = 0
    if xy > xhi/2:
        xy -= xhi
        lxy = -1
    elif xy < -xhi/2:
        xy += xhi
        lxy = 1
    lxz = 0
    if xz > xhi/2:
        xz -= xhi
        lxz = -1
    elif xz < -xhi/2:
        xz += xhi
        lxz = 1
    lyz = 0
    if yz > yhi/2:
        yz -= yhi
        lyz = -1
    elif yz < -yhi/2:
        yz += yhi
        lyz = 1
    a1 = np.array(hmat[:,0])
    a2 = np.array(hmat[:,1])
    a3 = np.array(hmat[:,2])
    vol = abs(np.dot(a1,np.cross(a2,a3)))
    a23 = np.cross(a2,a3)
    a31 = np.cross(a3,a1)
    a12 = np.cross(a1,a2)
    amat = np.zeros((3,3),dtype=float)
    amat[0,:] = a23[:]
    amat[1,:] = a31[:]
    amat[2,:] = a12[:]
    b1 = np.array((x,0.0,0.0))
    b2 = np.array((xy,y,0.0))
    b3 = np.array((xz,yz,z))
    bmat = np.zeros((3,3),dtype=float)
    bmat[:,0] = b1[:]
    bmat[:,1] = b2[:]
    bmat[:,2] = b3[:]
    if len(spos.shape) == 1:  # only one atom
        pos = np.zeros(spos.shape,dtype=float)
        newspos = shift_spos_for_lammps(spos,lxy,lxz,lyz,x,y,z,yz,xz,xy)
        pos = np.dot(hmat,newspos)
        pos = np.dot(bmat,np.dot(amat,pos))/vol
    elif len(spos.shape) == 2:  # array of atoms
        pos = np.zeros(spos.shape,dtype=float)
        for i,sp in enumerate(spos):
            newspos = shift_spos_for_lammps(sp,lxy,lxz,lyz,x,y,z,yz,xz,xy)
            pos[i] = np.dot(hmat,newspos)
            pos[i] = np.dot(bmat,np.dot(amat,pos[i]))/vol
        
    return xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,pos

def pbc(x):
    if x < 0.:
        return x -int(x) +1.0
    elif x >= 1.0:
        return x -int(x)
    else:
        return x

def shift_spos_for_lammps(spos,lxy,lxz,lyz,x,y,z,yz,xz,xy):
    import copy
    xyp = xy -lxy*x
    new_spos = copy.deepcopy(spos)
    new_spos[1] -= lyz*spos[2]
    new_spos[0] = new_spos[0] -lxz*spos[2] \
                  +(spos[1]*xyp -new_spos[1]*xy)/x
    for i in range(3):
        new_spos[i] = pbc(new_spos[i])
    return new_spos

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


def analyze(nsys):
    import elements
    from units import amu_to_g, Ang_to_cm
    a1 = nsys.a1 *nsys.alc
    a2 = nsys.a2 *nsys.alc
    a3 = nsys.a3 *nsys.alc
    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    c = np.linalg.norm(a3)
    vol = nsys.volume()
    alpha = np.arccos(np.dot(a2,a3)/b/c)/np.pi*180.0
    beta  = np.arccos(np.dot(a1,a3)/a/c)/np.pi*180.0
    gamma = np.arccos(np.dot(a1,a2)/a/b)/np.pi*180.0
    print(' a1 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]'.format(a1[0],
                                                                  a1[1],
                                                                  a1[2]))
    print(' a2 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]'.format(a2[0],
                                                                  a2[1],
                                                                  a2[2]))
    print(' a3 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]'.format(a3[0],
                                                                  a3[1],
                                                                  a3[2]))
    print(' a = {0:10.3f} A'.format(a))
    print(' b = {0:10.3f} A'.format(b))
    print(' c = {0:10.3f} A'.format(c))
    print(' alpha = {0:7.2f} deg.'.format(alpha))
    print(' beta  = {0:7.2f} deg.'.format(beta))
    print(' gamma = {0:7.2f} deg.'.format(gamma))
    print(' volume= {0:10.3f} A^3'.format(vol))
    print(' number of atoms   = ',nsys.num_atoms())
    if nsys.specorder:
        print(' number of atoms per species:')
        nspcs = nsys.num_species()
        mass = 0.0
        for i,s in enumerate(nsys.specorder):
            print('   {0:s}: {1:d}'.format(s,nspcs[i]))
            mass += nspcs[i]*elements.elements[s]['mass']
        print(' density = {0:7.2f} g/cm^3'.format(mass*amu_to_g
                                                 /(vol*Ang_to_cm**3) ))

if __name__ == "__main__":

    args= docopt(__doc__)

    infmt= args['--in-format']
    outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']
    scalefactor= args['--scale']
    shift= [ float(s) for s in args['--shift'].split(',') ]
    ncycle= int(args['--cycle-coord'])
    specorder= args['--specorder'].split(',')
    if specorder == 'None' or 'None' in specorder:
        specorder = []
    copies= [ int(i) for i in args['--periodic-copy'].split(',') ]
    charges= args['--charges']
    if charges == 'None':
        charges = []
    else:
        charges = [ float(c) for c in charges.split(',') ]

    nsys= NAPSystem(fname=infname,ffmt=infmt,specorder=specorder,
                    charges=charges)

    nsys.shift_atoms(*shift)
    if ncycle > 0:
        nsys.cycle_coord(ncycle)
    
    #...Periodic copy if needed
    copy_needed = False
    for c in copies:
        if c != 1:
            copy_needed = True
            break
        elif c < 1:
            raise ValueError(' Periodic copy was wrong. It should be >= 1.')
    if copy_needed:
        nsys.repeat(copies[0],copies[1],copies[2])

    if args['analyze']:
        analyze(nsys)

    elif args['convert']:
        if scalefactor != "None":
            nsys.alc *= float(scalefactor)

        nsys.write(fname=outfname,fmt=outfmt)

    else:
        raise NotImplementedError()
