"""
Reading and writing of pmd files for ASE Atoms object.
"""

import os
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms, FixScaled

def read_pmd(fname='pmd00000',specorder=[],fmvs=[]):
    """
    Import pmd format file.

    fname
        Name of the file to be read.
    specorder
        Order of species appeared in pmd file.
    fmvs
        Constraints of motion in bool of x,y,z direction.
        e.g., fmv = ((True,True,True),(True,True,False),)
        True means fix and False means free.
    """
    
    f = open(fname,'r')
    # 1st line: lattice_constant
    lattice_constant = float(f.readline().split()[0])
    # 2-4 lines: lattice vectors
    a = []
    for ii in range(3):
        s = f.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        a.append(floatvect)
    basis_vectors = np.array(a) * lattice_constant
    # 5-7 lines: velocities of lattice vectors
    tmp = f.readline()
    tmp = f.readline()
    tmp = f.readline()
    natm = int(f.readline().split()[0])
    positions = np.zeros((natm,3),dtype=float)
    symbols = []
    ifmvs = np.zeros((natm),dtype=int)
    maxifmv = 0
    for ia in range(natm):
        data = f.readline().split()
        tag = float(data[0])
        for ii in range(3):
            positions[ia,ii] = float(data[ii+1])
        sid,symbol,ifmv,aid = decode_tag(tag,specorder)
        symbols.append(symbol)
        ifmvs[ia] = ifmv
        maxifmv = max(maxifmv,ifmv)
    f.close()
    if maxifmv > len(fmvs)+1:
        raise ValueError('Length of fmvs are too short.')

    atoms = Atoms(symbols=symbols,cell=basis_vectors,pbc=True)
    atoms.set_scaled_positions(positions)
    #set constraints
    constraints = []
    indices = []
    for ia in range(natm):
        if ifmvs[ia] == 0:
            indices.append(ia)
        else:
            ifmv = ifmvs[ia]
            fmv = fmvs[ifmv-1]
            if all(fmv):
                indices.append(ia)
            elif any(fmv):
                constraints.append(FixScaled(atoms.get_cell(),ia,fmv))
    if indices:
        constraints.append(FixAtoms(indices))
    if constraints:
        atoms.set_constraint(constraints)
    return atoms

def write_pmd(atoms,fname='pmd00000',specorder=[]):
    """
    Write pmd format file from ASE Atoms object.
    """
    dname = os.path.dirname(fname)
    if not os.path.exists(dname):
        os.makedirs(dname)
    with open(fname,'w') as f:
        f.write(get_atom_conf_txt(atoms,specorder))


def get_atom_conf_txt(atoms,specorder=[]):
    txt= ''
    # no lattice constant in ASE
    txt+='  1.00000  \n'
    # cell vectors
    cell= atoms.cell
    txt += ' {0:12.7f}'.format(cell[0,0]) \
           +' {0:12.7f}'.format(cell[0,1]) \
           +' {0:12.7f}\n'.format(cell[0,2])
    txt += ' {0:12.7f}'.format(cell[1,0]) \
           +' {0:12.7f}'.format(cell[1,1]) \
           +' {0:12.7f}\n'.format(cell[1,2])
    txt += ' {0:12.7f}'.format(cell[2,0]) \
           +' {0:12.7f}'.format(cell[2,1]) \
           +' {0:12.7f}\n'.format(cell[2,2])
    txt += ' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(0.0,0.0,0.0)
    txt += ' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(0.0,0.0,0.0)
    txt += ' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(0.0,0.0,0.0)
    # num of atoms
    txt += ' {0:10d}\n'.format(len(atoms))
    # extract unique constraints from atoms.constraints
    fmvs,ifmvs = get_fmvs(atoms)
    # atom positions
    spos= atoms.get_scaled_positions()
    if not specorder:
        specorder = uniq(atoms.get_chemical_symbols())
        specorder.sort()
    for i in range(len(atoms)):
        atom= atoms[i]
        ifmv = ifmvs[i]
        txt += ' {0:s}'.format(get_tag(specorder,atom.symbol,i+1,ifmv))
        txt += ' {0:12.7f} {1:12.7f} {2:12.7f}'.format(spos[i,0],
                                                       spos[i,1],
                                                       spos[i,2])
        txt += ' 0.0 0.0 0.0'
        txt += ' 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n'
    return txt


def get_tag(specorder,symbol,atom_id,ifmv):
    sid= specorder.index(symbol)+1
    tag= float(sid) +ifmv*0.1 +atom_id*1e-14
    return '{0:16.14f}'.format(tag)

def decode_tag(tag,specorder):
    sid = int(tag)
    symbol = specorder[sid-1]
    ifmv = int((tag - sid)*10)
    atom_id = int((tag - sid - ifmv*0.1)*1e+14)
    return sid,symbol,ifmv,atom_id

def constraint2fmv(constraint):
    fmv = [1.0,1.0,1.0]
    for ii in range(3):
        if constraint[ii]:
            fmv[ii] = 0.0
    return fmv

def fmv2constraint(fmv):
    constraint = [False,False,False]
    for ii in range(3):
        if fmv[ii] < 0.01:
            constraint[ii] = True
    return constraint

def get_fmvs(atoms):
    """
    Extract unique constraints from atoms.constraints and
    return fmvs and ifmvs.
    """
    # extract unique constraints from atoms.constraints
    constraints = []
    constraints.append([False,False,False])
    ifmvs = np.zeros((len(atoms)),dtype=int)
    ifmvs[:] = 1
    if atoms.constraints:
        for cnst in atoms.constraints:
            if isinstance(cnst, FixScaled):
                mask= cnst.mask
            for i,c in enumerate(constraints):
                matched = mask == c
                if all(matched):
                    ifmvs[cnst.a] = i+1
                    break
            if not all(matched):
                constraints.append(mask)
                ifmvs[cnst.a] = len(constraints)
    # convert constraints to fmv
    fmvs = []
    for c in constraints:
        fmv = np.array((1.0, 1.0, 1.0))
        for ii in range(3):
            if c[ii]:
                fmv[ii] = 0.0
        fmvs.append(fmv)    
    return fmvs,ifmvs

def uniq(lst):
    newlst= []
    for l in lst:
        if not l in newlst:
            newlst.append(l)
    return newlst

