#!/usr/bin/env python
"""
IO functions for NAPSystem object.

Usage:
  io.py [options]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
import copy

from nappy.napsys import NAPSystem
from nappy.util import get_tag, decode_tag, pbc, \
    scaled_to_cartesian, cartesian_to_scaled

__author__ = "RYO KOBAYASHI"
__version__ = ""

FILE_FORMATS = ('pmd','POSCAR','CONTCAR','dump','xsf','lammps',
                'cube','CHGCAR','pdb')

def write(nsys,fname="pmdini",format=None,**kwargs):
    if format in (None,'None'):
        format= parse_filename(fname)

    if format == 'pmd':
        write_pmd(nsys,fname)
    elif format == 'POSCAR':
        write_POSCAR(nsys,fname)
    elif format == 'dump':
        write_dump(nsys,fname)
    elif format == 'xsf':
        write_xsf(nsys,fname)
    elif format == 'lammps':
        if hasattr(nsys, 'charges') and len(nsys.charges) > 0:
            write_lammps_data(nsys,fname,atom_style='charge')
        else:
            write_lammps_data(nsys,fname)
    elif format == 'cube':
        write_cube(nsys,fname,**kwargs)
    elif format in ('pdb','PDB'):
        import ase.io
        ase.io.write(filename=fname,images=nsys.to_ase_atoms,
                     format='proteindatabank')
    else:
        raise IOError('Cannot write out in the given format: '+format)

    return None

def read(fname="pmdini",format=None,specorder=None):
    if format in (None, 'None'):
        format= parse_filename(fname)

    if format == 'pmd':
        nsys = read_pmd(fname,specorder=specorder)
    elif format in ('POSCAR','CONTCAR','vasp','VASP'):
        nsys = read_POSCAR(fname,specorder=specorder)
    elif format == 'CHGCAR':
        nsys = read_CHGCAR(fname,specorder=specorder)
    elif format == 'dump':
        nsys = read_dump(fname,specorder=specorder)
    elif format == 'xsf':
        nsys = read_xsf(fname,specorder=specorder)
    elif format == 'lammps':
        nsys = read_lammps_data(fname,specorder=specorder)
    elif format == 'cube':
        nsys = read_cube(fname,specorder=specorder)
    else:
        print('Since the file format is unknown, try to read the file using ASE.')
        try:
            import ase.io
            atoms = ase.io.read(fname)
            nsys = from_ase(atoms)
        except Exception as e:
            print(' Failed to load input file even with ase.')
            raise
    return nsys

def read_pmd(fname='pmdini',specorder=None):
    nsys = NAPSystem()
    if specorder is not None:
        nsys.specorder = specorder
    incatm = 0
    with open(fname,'r') as f:
        iline = 0
        symbol = None
        for line in f.readlines():
            if line[0] in ('#','!'):  # comment line
                if 'specorder:' in line:  # overwrite specorder if specified in file
                    data = line.split()
                    specorder = [ d for d in data[2:len(data)]]
                    if nsys.specorder and set(nsys.specorder)!=set(specorder):
                        print(' WARNING: specorders are inconsistent, '
                              +'use one in the file.')
                    nsys.specorder = specorder
            else:
                if nsys.specorder is None or len(nsys.specorder) == 0:
                    raise ValueError('Specorder must be specified via the file or an argument.')
                iline = iline +1
                data = line.split()
                # 1st: lattice constant
                if iline == 1:
                    nsys.alc= float(data[0])
                # 2nd-4th: cell vectors
                elif iline == 2:
                    nsys.a1= np.array([float(x) for x in data])
                elif iline == 3:
                    nsys.a2= np.array([float(x) for x in data])
                elif iline == 4:
                    nsys.a3= np.array([float(x) for x in data])
                # 5th-7th: velocity of cell vectors
                elif 5 <= iline <= 7:
                    pass
                # 8th: num of atoms
                elif iline == 8:
                    natm = int(data[0])
                    sids = [ 0 for i in range(natm) ]
                    # poss = [ np.zeros(3) for i in range(natm) ]
                    # vels = [ np.zeros(3) for i in range(natm) ]
                    # frcs = [ np.zeros(3) for i in range(natm) ]
                    poss = np.zeros((natm,3))
                    vels = np.zeros((natm,3))
                    frcs = np.zeros((natm,3))
                # 9th-: atom positions
                else:
                    if incatm > natm:
                        break
                    fdata = [float(x) for x in data]
                    tag = fdata[0]
                    sid,ifmv,num = decode_tag(tag)
                    # poss[incatm][:] = fdata[1:4]
                    # vels[incatm][:] = fdata[4:7]
                    poss[incatm,:] = fdata[1:4]
                    vels[incatm,:] = fdata[4:7]
                    sids[incatm] = sid
                    incatm += 1
    nsys.atoms[['x','y','z']] = poss
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs
    nsys.atoms['sid'] = sids
    return nsys

def write_pmd(nsys,fname='pmdini'):
    f=open(fname,'w')
    if nsys.specorder and len(nsys.specorder)> 0:
        f.write("!\n")
        f.write("!  specorder: ")
        for s in nsys.specorder:
            f.write(" {0:<3s}".format(s))
        f.write("\n")
        f.write("!\n")
    # lattice constant
    f.write(" {0:15.9f}\n".format(nsys.alc))
    # cell vectors
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(*nsys.a1))
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(*nsys.a2))
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(*nsys.a3))
    # velocities of cell vectors
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
    f.write(" {0:19.15f} {1:19.15f} {2:19.15f}\n".format(0.0, 0.0, 0.0))
    # num of atoms
    f.write(" {0:10d}\n".format(len(nsys.atoms)))
    # atom positions
    poss = nsys.get_scaled_positions()
    vels = nsys.get_velocities()
    sids = nsys.atoms.sid
    if 'ifmv' not in nsys.atoms.columns:
        ifmvs = [ 1 for i in range(len(poss)) ]
    else:
        ifmvs = nsys.atoms.ifmv
    for i in range(len(nsys.atoms)):
        pi = poss[i]
        vi = vels[i]
        sid = sids[i]
        ifmv = ifmvs[i]
        tag = get_tag(sid,ifmv,i+1)  # assuming ifmv=1
        f.write(" {0:22.14e}".format(tag)
                +"  {0:19.15f} {1:19.15f} {2:19.15f}".format(*pi)
                +"  {0:8.4f}  {1:8.4f}  {2:8.4f}".format(*vi)
                +"\n")
    f.close()
    return None

def read_POSCAR(fname='POSCAR',specorder=None):
    nsys = NAPSystem()
    with open(fname,'r') as f:
        # 1st line: comment
        f.readline()
        # 2nd: lattice constant
        nsys.alc= float(f.readline().split()[0])
        # 3rd-5th: cell vectors
        nsys.a1= np.array([float(x) for x in f.readline().split()])
        nsys.a2= np.array([float(x) for x in f.readline().split()])
        nsys.a3= np.array([float(x) for x in f.readline().split()])
        # 6th: species names or number of each species
        buff= f.readline().split()
        if not buff[0].isdigit():
            spcs = copy.deepcopy(buff)
            buff= f.readline().split()
            if specorder is None:
                nsys.specorder = spcs
            else:
                nsys.specorder = specorder
                for s in spcs:
                    if s not in nsys.specorder:
                        nsys.specorder.append(s)
        num_species= np.array([ int(n) for n in buff])
        try:
            spcs
        except NameError:
            spcs = nsys.specorder
        #...Check number of species in POSCAR file and in specorder
        if len(num_species) > len(nsys.specorder):
            msg = '''
ers of species in POSCAR is greater than the one in specorder, which should be the same or less.
er of species in POSCAR = {0:d}
need to specify the species order correctly with --specorder option.
            '''.format(len(num_species))
            raise ValueError(msg)
        natm = np.sum(num_species)
        sids = [ 0 for i in range(natm) ]
        # poss = [ np.zeros(3) for i in range(natm) ]
        # vels = [ np.zeros(3) for i in range(natm) ]
        # frcs = [ np.zeros(3) for i in range(natm) ]
        poss = np.zeros((natm,3))
        vels = np.zeros((natm,3))
        frcs = np.zeros((natm,3))
        #print("Number of atoms = {0:5d}".format(natm))
        # 7th or 8th line: comment
        c7= f.readline()
        if c7[0] in ('s','S'):
            c7= f.readline()
        if c7[0] in ('c','C'):  # positions are in Cartesian coordinate
            hi = nsys.get_hmat_inv()
            coord = 'cartesian'
        else:  # such as "Direct"
            coord = 'scaled'
        
        #...Atom positions
        for i in range(natm):
            buff= f.readline().split()
            sid= 1
            m= 0
            sindex=0
            for n in num_species:
                m += n
                if i < m:
                    if spcs and nsys.specorder:
                        sid = nsys.specorder.index(spcs[sindex]) + 1
                    break
                sid += 1
                sindex += 1
            sids[i] = sid
            pos = [ float(buff[0]), float(buff[1]), float(buff[2])]
            if coord == 'cartesian':
                x1,x2,x3 = cartesian_to_scaled(hi,pos[0],pos[1],pos[2])
            elif coord == 'scaled':
                x1,x2,x3 = pos[0],pos[1],pos[2]
            poss[i,:] = [x1,x2,x3]

    nsys.atoms[['x','y','z']] = poss
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs
    nsys.atoms['sid'] = sids
    return nsys
            
def write_POSCAR(nsys,fname='POSCAR'):
    from datetime import datetime
    f=open(fname,'w')
    f.write('Generated by napsys.py at {0}.'.format(datetime.now().strftime('%Y-%m-%d')))
    if nsys.specorder:
        for s in nsys.specorder:
            f.write(' '+s)
    f.write('\n')
    # lattice vector
    f.write(' {0:15.7f}\n'.format(nsys.alc))
    # cell vectors
    f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(*nsys.a1))
    f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(*nsys.a2))
    f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(*nsys.a3))
    # count num of atoms per species
    num_species= nsys.natm_per_species()
    # if specorder is defined, write species names
    if nsys.specorder:
        for i,ns in enumerate(num_species):
            s = nsys.specorder[i]
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
    outorder = np.zeros(len(nsys.atoms),dtype=int)
    inc = 0
    for isp,spc in enumerate(nsys.specorder):
        for ia in range(len(nsys.atoms)):
            sid = nsys.atoms.sid[ia]
            if sid == isp +1:
                outorder[inc] = ia
                inc += 1
    if len(outorder) != len(nsys.atoms):
        print('len(outorder),natm=', len(outorder),len(nsys.atoms))
        raise ValueError(' len(outorder) != natm')
    spos = nsys.get_scaled_positions()
    for ia in outorder:
        pi = spos[ia]
        f.write(' {0:15.7f} {1:15.7f} {2:15.7f} T T T\n'.format(*pi))
    f.close()
    return None

def read_dump(fname="dump",specorder=None):
    nsys = NAPSystem()
    f=open(fname,'r')
    mode= 'None'
    ixyz= 0
    iatm= 0
    natm= -1
    symbol = None
    if specorder is None:
        nsys.specorder = []
    else:
        nsys.specorder = specorder
    nsys.alc = 1.0
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
    ivx = -1
    ivy = -1
    ivz = -1
    ifx = -1
    ify = -1
    ifz = -1
    for line in f.readlines():
        data = line.split()
        if 'ITEM' in line:
            if 'NUMBER OF ATOMS' in line:
                mode= 'NUMBER OF ATOMS'
                continue
            elif 'BOX BOUNDS' in line:
                mode= 'BOX BOUNDS'
                continue
            elif 'ATOMS' in line:
                mode= 'ATOMS'
                aux_names = [ name for i,name in enumerate(data) if i > 1 ]
                aux_names.remove('id')
                aux_names.remove('type')
                if ('x' not in aux_names and 'xu' not in aux_names) or \
                   ('y' not in aux_names and 'zu' not in aux_names) or \
                   ('z' not in aux_names and 'zu' not in aux_names):
                    raise ValueError('Not enough coordinate info.\nCheck the dump file format.')
                try:
                    ix = aux_names.index('x') +2
                except Exception:
                    ix = aux_names.index('xu') +2
                try:
                    iy = aux_names.index('y') +2
                except Exception:
                    iy = aux_names.index('yu') +2
                try:
                    iz = aux_names.index('z') +2
                    # iauxstart = 5
                except Exception:
                    iz = aux_names.index('zu') +2
                    # iauxstart = 5
                try:
                    ivx = aux_names.index('vx') +2
                    ivy = aux_names.index('vy') +2
                    ivz = aux_names.index('vz') +2
                    # iauxstart = 8
                except Exception:
                    pass
                try:
                    ifx = aux_names.index('fx') +2
                    ify = aux_names.index('fy') +2
                    ifz = aux_names.index('fz') +2
                except Exception:
                    pass
                # for s in ('x','xu','y','yu','z','zu','vx','vy','vz'):
                #     if s in aux_names:
                #         aux_names.remove(s)
                if len(aux_names)>0:
                    auxs = np.zeros((natm,len(aux_names)))
                continue
            elif 'TIMESTEP' in line:
                mode= 'TIMESTEP'
                continue
            
        if mode == 'TIMESTEP':
            timestep = int(data[0])
        elif mode == 'NUMBER OF ATOMS':
            natm= int(data[0])
            sids = [ 0 for i in range(natm) ]
            # poss = [ np.zeros(3) for i in range(natm) ]
            # vels = [ np.zeros(3) for i in range(natm) ]
            # frcs = [ np.zeros(3) for i in range(natm) ]
            poss = np.zeros((natm,3))
            vels = np.zeros((natm,3))
            frcs = np.zeros((natm,3))
        elif mode == 'BOX BOUNDS':
            if ixyz == 0:
                xlo_bound= float(data[0])
                xhi_bound= float(data[1])
                if len(data) > 2:
                    xy = float(data[2])
            elif ixyz == 1:
                ylo_bound= float(data[0])
                yhi_bound= float(data[1])
                if len(data) > 2:
                    xz = float(data[2])
            elif ixyz == 2:
                zlo_bound= float(data[0])
                zhi_bound= float(data[1])
                if len(data) > 2:
                    yz = float(data[2])
            ixyz += 1
            if ixyz > 2:
                xlo = xlo_bound -min(0.0,xy,xz,xy+xz)
                xhi = xhi_bound -max(0.0,xy,xz,xy+xz)
                ylo = ylo_bound -min(0.0,yz)
                yhi = yhi_bound -max(0.0,yz)
                zlo = zlo_bound
                zhi = zhi_bound
                #...Original definition of lattice vectors could be different
                #   from this, because the definition in dump format
                #   requires y,z-components of vector a1 to be zero.
                nsys.a1 = np.array([xhi-xlo, 0., 0.],dtype=float)
                nsys.a2 = np.array([xy, yhi-ylo, 0.],dtype=float)
                nsys.a3 = np.array([xz, yz, zhi-zlo],dtype=float)
                hmat = nsys.get_hmat()
                hmati= nsys.get_hmat_inv()
        elif mode == 'ATOMS':
            if iatm < natm:
                symbol = None
                if data[1].isdigit():
                    sid = int(data[1])
                    sids[iatm] = sid
                    symbol = nsys.specorder[sid-1]
                else:
                    symbol = data[1]
                    if symbol not in nsys.specorder:
                        nsys.specorder.append(symbol)
                    sid = nsys.specorder.index(symbol)+1
                    sids[iatm] = sid
                r0 = [ float(data[ix]), float(data[iy]), float(data[iz]) ]
                if ivx > 0 and ivy > 0 and ivz > 0:
                    v0 = [ float(data[ivx]),float(data[ivy]),float(data[ivz])]
                else:
                    v0 = [ 0., 0., 0. ]
                if ifx > 0 and ify > 0 and ifz > 0:
                    f0 = [ float(data[ifx]),float(data[ify]),float(data[ifz])]
                else:
                    f0 = [ 0., 0., 0. ]
                sr = np.dot(hmati,r0)
                sv = np.dot(hmati,v0)
                sr[0] = pbc(sr[0])
                sr[1] = pbc(sr[1])
                sr[2] = pbc(sr[2])
                # poss[iatm][:] = sr[:]
                # vels[iatm][:] = sv[:]
                poss[iatm,:] = sr[:]
                vels[iatm,:] = sv[:]
                frcs[iatm,:] = f0[:]
                
                if len(aux_names)>0:
                    # auxs[iatm,:] = [ float(x) for x in data[iauxstart:] ]
                    auxs[iatm,:] = [ float(x) for x in data[2:] ]

            iatm += 1

    nsys.atoms[['x','y','z']] = poss
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs
    nsys.atoms['sid'] = sids
    for ia in range(len(aux_names)):
        name = aux_names[ia]
        if name in ('x','xu','y','yu','z','zu','vx','vy','vz','fx','fy','fz'):
            continue
        aux = auxs[:,ia]
        nsys.atoms[name] = aux.tolist()
    f.close()
    return nsys

def write_dump(nsys,fname='dump',auxs=['vx','vy','vz']):
    """
    Write LAMMPS dump format file.
    """
    f= open(fname,'w')
    f.write("ITEM: TIMESTEP\n")
    f.write("0\n")
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write("{0:d}\n".format(len(nsys.atoms)))

    hmat = nsys.get_hmat()
    poss = nsys.get_scaled_positions()
    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,rposs = to_lammps(hmat,poss)
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
    
    auxs_exist = nsys.get_aux_names()
    aux_names = []
    for aux in auxs:
        if aux in ('vx','vy','vz',*auxs_exist):
            aux_names.append(aux)
    f.write("ITEM: ATOMS id type x y z")
    if len(aux_names) > 0:
        for aux in aux_names:
            f.write(' {0:s}'.format(aux))
    f.write("\n")

    vels = nsys.get_scaled_velocities()
    frcs = nsys.get_scaled_forces()

    if len(aux_names)>0:
        for i in range(len(nsys.atoms)):
            rpos = rposs[i]
            #...NOTE: velocity is scaled value here,
            #   if one wants to get real velocity, one needs to convert
            #   is in to_lammps function, not here.
            vel = vels[i]
            sid = nsys.atoms.sid[i]
            symbol = nsys.specorder[sid-1]
            f.write("{0:8d} {1:3s} ".format(i+1,symbol))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(*rpos))
            for name in aux_names:
                if name == 'vx':
                    f.write("{0:8.3f} ".format(vel[0]))
                elif name == 'vy':
                    f.write("{0:8.3f} ".format(vel[1]))
                elif name == 'vz':
                    f.write("{0:8.3f} ".format(vel[2]))
                else:
                    f.write("{0:11.3e}".format(nsys.atoms[name][i]))
            f.write("\n")
    else:
        for i in range(len(nsys.atoms)):
            rpos = rposs[i]
            vel = vels[i]
            sid = nsys.atoms.sid[i]
            symbol = nsys.specorder[sid-1]
            f.write("{0:8d} {1:3s} ".format(i+1,symbol))
            f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(*rpos))
            f.write("{0:8.3f} {1:8.3f} {2:8.3f} ".format(*vel))
            f.write("\n")
    f.close()
    return None

def read_lammps_data(fname="data.lammps",atom_style='atomic',specorder=None):
    nsys = NAPSystem()
    if specorder is None:
        nsys.specorder = []
    else:
        nsys.specorder = specorder
    f=open(fname,'r')
    mode= 'None'
    iatm= 0
    symbol = None
    nsys.alc= 1.0
    xy = 0.0
    xz = 0.0
    yz = 0.0
    for line in f.readlines():
        data = line.split()
        if mode == 'None':
            if 'atoms' in line:
                natm = int(data[0])
                sids = [ 0 for i in range(natm) ]
                poss = np.zeros((natm,3))
                vels = np.zeros((natm,3))
                frcs = np.zeros((natm,3))
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
                # nsys.a1 = np.array([xhi-xlo,xy,xz],dtype=float)
                # nsys.a2 = np.array([0.0,yhi-ylo,yz],dtype=float)
                # nsys.a3 = np.array([0.0,0.0,zhi-zlo],dtype=float)
                nsys.a1 = np.array([xhi-xlo,0.0,0.0],dtype=float)
                nsys.a2 = np.array([xy,yhi-ylo,0.0],dtype=float)
                nsys.a3 = np.array([xz,yz,zhi-zlo],dtype=float)
                hmat = nsys.get_hmat()
                hmati= np.linalg.inv(hmat)
                continue
        elif mode == 'Atoms':
            if len(data) >= 5 and iatm < natm:
                idat = 0
                # ai = Atom()
                idat += 1
                # ai.set_sid(int(data[idat]))
                sid = int(data[idat])
                sids[iatm] = sid
                if nsys.specorder:
                    symbol = nsys.specorder[sid-1]
                # if symbol and ai.symbol != symbol:
                #     ai.set_symbol(symbol)
                # if atom_style == 'charge': 
                #     idat += 1
                #     chg = float(data[idat])
                #     # ai.set_aux('charge',chg)
                #     if aux_names and 'charge' not in aux_names:
                        
                idat += 1
                x0= float(data[idat])
                idat += 1
                y0= float(data[idat])
                idat += 1
                z0= float(data[idat])
                x = hmati[0,0]*x0 +hmati[0,1]*y0 +hmati[0,2]*z0
                y = hmati[1,0]*x0 +hmati[1,1]*y0 +hmati[1,2]*z0
                z = hmati[2,0]*x0 +hmati[2,1]*y0 +hmati[2,2]*z0
                x = pbc(x)
                y = pbc(y)
                z = pbc(z)
                poss[iatm,:] = [x,y,z]
                iatm += 1
    f.close()
    nsys.atoms[['x','y','z']] = poss
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs
    nsys.atoms['sid'] = sids
    return nsys

def write_lammps_data(nsys,fname='data.lammps',atom_style='atomic'):
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
    f.write("{0:d}  atoms\n".format(len(nsys.atoms)))
    f.write("{0:d}  atom types\n".format(nsys.num_species()))
    f.write('\n')
    hmat = nsys.get_hmat()
    poss = nsys.get_scaled_positions()
    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,rposs = to_lammps(hmat,poss)
    f.write("{0:20.10f} {1:20.10f} xlo xhi\n".format(xlo,xhi))
    f.write("{0:20.10f} {1:20.10f} ylo yhi\n".format(ylo,yhi))
    f.write("{0:20.10f} {1:20.10f} zlo zhi\n".format(zlo,zhi))
    if abs(xy) > 1e-8 or abs(xz) > 1e-8 or abs(yz) > 1e-8:
        f.write("{0:20.10f} {1:20.10f} {2:20.10f} xy xz yz\n".format(xy,xz,yz))
    f.write("\n")
    f.write("Atoms\n")
    f.write("\n")
    rpos = nsys.get_real_positions()
    for i in range(len(nsys.atoms)):
        pi = rpos[i]
        sid = nsys.atoms.sid[i]
        f.write("{0:8d} {1:3d} ".format(i+1,sid))
        # if atom_style == 'charge':
        #     f.write('{0:10.4f} '.format(nsys.charges[ai.sid-1]))
        f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(*pi))
        f.write("\n")
    f.close()
    return None

def read_xsf(fname="xsf",specorder=None):
    from nappy.elements import get_symbol_from_number
    nsys = NAPSystem()
    if specorder is None:
        nsys.specorder = []
    else:
        nsys.specorder = specorder
    f=open(fname,'r')
    mode= 'None'
    ixyz= 0
    iatm= 0
    natm= -1
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
            hi = nsys.get_hmat_inv()
            # print 'Inversed h-matrix:'
            # print hi
            continue
        
        if mode == 'CRYSTAL':
            pass
        elif mode == 'PRIMVEC':
            if ixyz == 0:
                arr = [ float(x) for x in line.split() ]
                nsys.a1[0] = arr[0]
                nsys.a1[1] = arr[1]
                nsys.a1[2] = arr[2]
            elif ixyz == 1:
                arr = [ float(x) for x in line.split() ]
                nsys.a2[0] = arr[0]
                nsys.a2[1] = arr[1]
                nsys.a2[2] = arr[2]
            elif ixyz == 2:
                arr = [ float(x) for x in line.split() ]
                nsys.a3[0] = arr[0]
                nsys.a3[1] = arr[1]
                nsys.a3[2] = arr[2]
            ixyz += 1
        elif mode == 'PRIMCOORD':
            data = line.split()
            if len(data) == 1:
                natm= int(data[0])
                sids = [ 0 for i in range(natm) ]
                poss = np.zeros((natm,3))
                vels = np.zeros((natm,3))
                frcs = np.zeros((natm,3))
                continue
            elif len(data) == 2:
                natm= int(data[0])
                nspcs= int(data[1])
                sids = [ 0 for i in range(natm) ]
                poss = np.zeros((natm,3))
                vels = np.zeros((natm,3))
                frcs = np.zeros((natm,3))
                continue
            elif len(data) == 4 or len(data) == 7:
                if iatm >= natm:
                    continue
                symbol = get_symbol_from_number(int(data[0]))
                if symbol not in nsys.specorder:
                    nsys.specorder.append(symbol)
                sid = nsys.specorder.index(symbol) +1
                sids[iatm] = sid
                # ai.set_sid(sid)
                xc= float(data[1])
                yc= float(data[2])
                zc= float(data[3])
                xi,yi,zi = cartesian_to_scaled(hi,xc,yc,zc)
                poss[iatm,:] = [xi,yi,zi]
                # print 'iatm,symbol,sid,xc,yc,zc = ',iatm,symbol,sid,xc,yc,zc
            else:
                continue
            iatm += 1
    nsys.alc= 1.0
    nsys.atoms[['x','y','z']] = poss
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs
    nsys.atoms['sid'] = sids
    f.close()
    return nsys

def write_xsf(nsys,fname='xsf'):
    """
    Write XCrysden xsf format.
    Only applicable to orthogonal system.
    """
    from nappy.elements import get_number_from_symbol
    if not nsys.specorder:
        raise ValueError('Specorder has to be defined to write'
                         +' xsf format file.')
    h = np.zeros((3,3),dtype=float)
    h[0,:] = nsys.a1[:]
    h[1,:] = nsys.a2[:]
    h[2,:] = nsys.a3[:]
    f= open(fname,'w')
    f.write("CRYSTAL\n")
    f.write("PRIMVEC\n")
    f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(*nsys.a1))
    f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(*nsys.a2))
    f.write("{0:9.3f} {1:9.3f} {2:9.3f}\n".format(*nsys.a3))
    f.write("PRIMCOORD\n")
    f.write("{0:>8d}  {1:2d}\n".format(len(nsys.atoms),nsys.num_species()))
    spos = nsys.get_scaled_positions()
    for i in range(len(nsys.atoms)):
        pos = spos[i]
        x,y,z = scaled_to_cartesian(h,*pos)
        sid = nsys.atoms.sid[i]
        symbol = nsys.specorder[sid-1]
        number = get_number_from_symbol(symbol)
        f.write(" {0:3d} ".format(number))
        f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(x,y,z))
        #f.write("{0:12.5f} {1:12.5f} {2:12.5f} ".format(vx,vy,vz))
        f.write("\n")
    f.close()
    return None

def read_CHGCAR(fname='CHGCAR',specorder=None):
    """
    Read CHGCAR file and get information of cell, atoms, and volumetric data.

    Parameters
    ----------
    fname : str
         File name to be read.
    """
    nsys = NAPSystem()
    with open(fname,'r') as f:
        # 1st line: comment
        f.readline()
        # 2nd: lattice constant
        nsys.alc= float(f.readline().split()[0])
        # 3rd-5th: cell vectors
        nsys.a1= np.array([float(x) for x in f.readline().split()])
        nsys.a2= np.array([float(x) for x in f.readline().split()])
        nsys.a3= np.array([float(x) for x in f.readline().split()])
        # 6th: species names or number of each species
        buff= f.readline().split()
        if not buff[0].isdigit():
            spcs = copy.deepcopy(buff)
            buff= f.readline().split()
            if specorder is None:
                nsys.specorder = spcs
            else:
                nsys.specorder = specorder
                for s in spcs:
                    if s not in nsys.specorder:
                        nsys.specorder.append(s)
        num_species= np.array([ int(n) for n in buff])
        try:
            spcs
        except NameError:
            spcs = nsys.specorder
        #...Check number of species in POSCAR file and in specorder
        if len(num_species) > len(nsys.specorder):
            msg = '''
ers of species in POSCAR is greater than the one in specorder, which should be the same or less.
er of species in POSCAR = {0:d}
need to specify the species order correctly with --specorder option.
            '''.format(len(num_species))
            raise ValueError(msg)
        natm = np.sum(num_species)
        sids = [ 0 for i in range(natm) ]
        poss = np.zeros((natm,3))
        vels = np.zeros((natm,3))
        frcs = np.zeros((natm,3))
        #print("Number of atoms = {0:5d}".format(natm))
        # 7th or 8th line: comment
        c7= f.readline()
        if c7[0] in ('s','S'):
            c7= f.readline()
        if c7[0] in ('c','C'):  # positions are in Cartesian coordinate
            hi = nsys.get_hmat_inv()
            coord = 'cartesian'
        else:  # such as "Direct"
            coord = 'scaled'
        
        #...Atom positions
        for i in range(natm):
            buff= f.readline().split()
            sid= 1
            m= 0
            sindex=0
            for n in num_species:
                m += n
                if i < m:
                    if spcs and nsys.specorder:
                        sid = nsys.specorder.index(spcs[sindex]) + 1
                    break
                sid += 1
                sindex += 1
            sids[i] = sid
            pos = [ float(buff[0]), float(buff[1]), float(buff[2])]
            if coord == 'cartesian':
                x1,x2,x3 = cartesian_to_scaled(hi,pos[0],pos[1],pos[2])
            elif coord == 'scaled':
                x1,x2,x3 = pos[0],pos[1],pos[2]
            poss[i][:] = [x1,x2,x3]

        #...Up to here, code should be the same as read_POSCAR, in case of CHGCAR lines follow
        nsys.atoms[['x','y','z']] = poss
        nsys.atoms[['vx','vy','vz']] = vels
        nsys.atoms[['fx','fy','fz']] = frcs
        nsys.atoms['sid'] = sids

        #...Read volumetric data
        f.readline()  # should be one blank line
        ndiv = [ int(x) for x in f.readline().split() ]
        ndata = ndiv[0]*ndiv[1]*ndiv[2]
        voldata = np.zeros(ndata)
        inc = 0
        while True:
            line = [ float(d) for d in f.readline().split() ]
            for d in line:
                voldata[inc] = d
                inc += 1
            if inc >= ndata:
                break
        #...Stop reading CHGCAR file here, ignoring the rest of file
        nsys.ndiv = copy.copy(ndiv)
        nsys.voldata = np.reshape(voldata,ndiv,order='F')
    return nsys

def read_cube(fname, specorder=None):
    """Read Gaussian cube format file."""
    from nappy.elements import get_symbol_from_number
    from nappy.units import Bohr_to_Ang
    nsys = NAPSystem()
    if specorder is None:
        nsys.specorder = []
    else:
        nsys.specorder = specorder

    with open(fname, 'r') as f:
        lines = f.readlines()
    natm = -1
    nvdata = -1

    # 1,2-th lines: comments

    # 3rd line: natm, vorig(x,y,x)
    d = lines[2].split()
    natm = int(d[0])
    vorig = np.array([float(o) for o in d[1:4] ])
    # 4,5,6-th line: ndiv, ax, ay, az
    ndiv = np.zeros(3,dtype=int)
    dhmat = np.zeros((3,3),dtype=float)
    hmat = np.zeros((3,3),dtype=float)
    for i in range(3):
        d = lines[3+i].split()
        ndiv[i] = int(d[0])
        dhmat[i,:] = [ float(v) for v in d[1:4] ]
        hmat[i,:] = dhmat[i,:] *ndiv[i]
    nvoldat = ndiv[0]*ndiv[1]*ndiv[2]
    hmati = np.linalg.inv(hmat)
    # 7-(7+natm)-th lines: id, elem-ID, rh[0:3]
    rpos = np.zeros(3,dtype=float)  # real pos in Bohr
    sposs = np.zeros((natm,3),dtype=float)  # scaled poss
    sids = np.zeros(natm, dtype=int)
    for i in range(natm):
        d = lines[6+i].split()
        symbol = get_symbol_from_number(int(eval(d[1])))
        if symbol not in nsys.specorder:
            nsys.specorder.append(symbol)
        sid = nsys.specorder.index(symbol) +1
        sids[i] = sid
        rpos = [ float(p)*Bohr_to_Ang for p in d[2:5] ]
        sposs[i] = np.dot(hmati,rpos)
    nsys.alc = 1.0
    nsys.set_hmat(hmat)
    nsys.atoms[['x','y','z']] = sposs
    nsys.atoms[['vx','vy','vz']] = np.zeros((natm,3))
    nsys.atoms[['fx','fy','fz']] = np.zeros((natm,3))
    nsys.atoms[['sid']] = sids
    # hereafter, volumetric data if exists
    if nvoldat > 0:
        vdata = np.zeros(nvoldat, dtype=float)
        nlines = int(nvoldat/6) +1
        inc = 0
        for il in range(nlines):
            d = lines[6+natm+il].split()
            vdata[inc:inc+len(d)] = [ float(v) for v in d ]
            inc += len(d)
        nsys.voldata = vdata
        nsys.nvoldiv = ndiv
        nsys.volorig = vorig
    return nsys

def write_cube(nsys, fname='cube', origin=[0.,0.,0.]):
    """
    Output the system info including volumetric data to a file in Gaussian cube format.
    """
    txt = get_cube_txt(nsys,origin=origin)
    with open(fname,'w') as f:
        f.write(txt)
    return None

def get_cube_txt(nsys,origin=[0.,0.,0.]):
    """
    Conver the system info to Gaussian cube format for the purpose of visualization using py3Dmol.

    Inputs
    ------
    nsys : NAPSystem
        NAPSystem object to be written out.
    origin : list or tuple of length of 3
        Shift of origin in Angstrom.

    Returns
    -------
    txt : str
        Test string of system information including volumetric data in cube format.
    """
    from datetime import datetime
    from nappy.units import Ang_to_Bohr
    txt = ''
    txt += 'Gaussian cube format for {0},\n'.format(nsys.get_chemical_formula())
    txt += 'generated by napsys.py at {0}.\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    #...Num of atoms, origin
    org_bohr = [ o*Ang_to_Bohr for o in origin ]
    txt += ' {0:8d} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(nsys.num_atoms(),*org_bohr)
    #...Num of divisions, lattice vectors
    a,b,c = nsys.get_lattice_vectors()
    if hasattr(nsys,'nvoldiv'):
        ndiv = nsys.nvoldiv
        ndata = ndiv[0]*ndiv[1]*ndiv[2]
        voldata = nsys.voldata.reshape((ndata,),order='C')
    else:
        ndiv = [ 1,1,1 ]
        voldata = [ 0. ]
    # Convert the unit from Angstrom to Bohr (Bohr is the default length unit in cube?)
    a = a /ndiv[0] *Ang_to_Bohr
    b = b /ndiv[1] *Ang_to_Bohr
    c = c /ndiv[2] *Ang_to_Bohr
        
    txt += ' {0:8d} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(ndiv[0],*a)
    txt += ' {0:8d} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(ndiv[1],*b)
    txt += ' {0:8d} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(ndiv[2],*c)
    

    #...Atoms
    from nappy.elements import elements
    rposs = nsys.get_real_positions()
    for i in range(len(nsys.atoms)):
        sid = nsys.atoms.loc[i,'sid']
        spc = nsys.specorder[sid-1]
        try:
            element = elements[spc]
            anum = element['number']
            val = float(element['valence'])
        except Exception as e:
            anum = 0
            val = 0.0
        txt += ' {0:5d} {1:12.6f}'.format(anum,val)
        pi = rposs[i]
        pi_bohr = [ pi[l]*Ang_to_Bohr +org_bohr[l] for l in range(3) ]
        txt += ' {0:12.6f} {1:12.6f} {2:12.6f}\n'.format(*pi_bohr)

    #...Volumetric data: 6 entries per line
    if hasattr(nsys,'nvoldiv'):
        inc = 0
        for vd in voldata:
            txt += ' {0:12.4e}'.format(vd)
            inc += 1
            if inc % 6 == 0:
                txt += '\n'
                inc = 0
        txt += '\n'
    else:
        txt += ' {0:12.4e}\n'.format(0.)
    return txt

def get_PDB_txt(nsys,**kwargs):
    """
    Convert the system info to Protein Data Bank (PDB) format
    for the purpose of visualization using py3Dmol (3Dmol.js).

    Inputs
    ------
    nsys : NAPSystem
        NAPSystem object to be output.

    Returns
    -------
    txt : str
        Text string of system information in PDB format.
    """
    max_num_atoms = 100000
    if 'max_num_atoms' in kwargs:
        max_num_atoms = kwargs['max_num_atoms']
    
    if len(nsys.atoms) >= max_num_atoms:
        raise ValueError('Number of atoms is too large for PDB format...')
    
    from datetime import datetime
    txt = ''
    txt += '{0:6s}    Generated by napsys.py at {1}\n'.format('COMPND',
                                                              datetime.now().strftime('%Y-%m-%d'))
    #...CRYST1
    a1,a2,a3 = nsys.get_lattice_vectors()
    a,b,c = nsys.get_lattice_lengths()
    alpha,beta,gamma = nsys.get_lattice_angles(unit='degree')
    txt += 'CRYST1' +'{0:9.3f}{1:9.3f}{2:9.3f}'.format(a,b,c) \
           +'{0:7.2f}{1:7.2f}{2:7.2f}'.format(alpha,beta,gamma) \
           +' P 1\n'
    #...HETATM
    rposs = nsys.get_real_positions()
    for i in range(len(nsys.atoms)):
        pi = rposs[i]
        sid= nsys.atoms.loc[i,'sid']
        spc= nsys.specorder[sid-1]
        txt += 'HETATM' +'{0:5d}'.format(i+1) +' ' \
               +'{0:4s}'.format(spc) +' ' \
               +'UNL     1    ' \
               +'{0:8.3f}{1:8.3f}{2:8.3f}'.format(*pi) \
               +'{0:6.2f}{1:6.2f}'.format(1.0,0.0) \
               +'          {0:2s}'.format(spc)
        if 'chg' in nsys.atoms.columns:
            chgi = nsys.atoms.loc[i,'chg']
            txt += '{0:2d}'.format(chgi)
        txt += '\n'

    #...MASTER
    txt += 'MASTER        0    0    0    0    0    0    0    0' \
           +'{0:5d}'.format(len(nsys.atoms))+'    0    0    0\n'
    txt += 'END\n'

    return txt

def from_ase(atoms,specorder=None):
    """
    Convert ASE Atoms object to NAPSystem object.
    """
    spcorder = []
    if specorder is not None:
        spcorder = specorder
    symbols = atoms.get_chemical_symbols()
    spos = atoms.get_scaled_positions()
    vels = atoms.get_velocities()
    cell = atoms.get_cell()
    celli = np.linalg.inv(cell)
    if spos is None:
        raise ValueError('ASE atoms object has no atom in it.')
    #...Initialize and remake self.specorder
    for s in symbols:
        if s not in spcorder:
            spcorder.append(s)
    nsys = NAPSystem(specorder=spcorder)
    # nsys = cls(specorder=spcorder)
    nsys.alc= 1.0
    nsys.a1[:] = atoms.cell[0]
    nsys.a2[:] = atoms.cell[1]
    nsys.a3[:] = atoms.cell[2]
    #...First, initialize arrays
    natm = len(atoms)
    sids = [ 0 for i in range(natm) ]
    poss = np.array(spos)
    if vels is None:
        vels = np.zeros((natm,3))
    else:
        vels = np.array(vels)
    frcs = np.zeros((natm,3))

    #...Create arrays to be installed into nsys.atoms
    sids = [ nsys.specorder.index(si)+1 for si in symbols ]
    nsys.atoms.sid = sids
    nsys.atoms[['x','y','z']] = poss
    for i in range(len(vels)):
        vels[i] = np.dot(celli,vels[i])
    nsys.atoms[['vx','vy','vz']] = vels
    nsys.atoms[['fx','fy','fz']] = frcs

    return nsys

def parse_filename(filename):
    for format in FILE_FORMATS:
        if format in filename:
            return format
    return None

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
    b1 = np.array((x , 0.0, 0.0))
    b2 = np.array((xy, y  , 0.0))
    b3 = np.array((xz, yz , z  ))
    bmat = np.zeros((3,3),dtype=float)
    bmat[:,0] = b1[:]
    bmat[:,1] = b2[:]
    bmat[:,2] = b3[:]
    if (spos == None).any() or len(spos) == 0:
        pos = None
    elif len(spos.shape) == 1:  # only one atom
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

def get_nglview(nsys):
    """
    Retern a nglview object via ase_atoms.
    """
    import nglview as nv
    return nv.show_ase(nsys.to_ase_atoms())

if __name__ == "__main__":

    args = docopt(__doc__)
    print(__file__)
