"""
Reading and writing of pmd files for ASE Atoms object.
"""

import numpy as np
from ase.constraints import FixScaled

def get_atom_conf_txt(atoms,specorder=[]):
    if not specorder:
        specorder = uniq(atoms.get_chemical_symbols())
        specorder.sort()
    # print 'atoms = ',atoms
    txt= ''
    #...specorder info as comment lines
    txt= '!\n'
    txt+='!  specorder '
    for s in specorder:
        txt += ' {0:s}'.format(s)
    txt += '\n'
    txt += '!\n'
    # no lattice constant in ASE
    txt+='  1.00000  \n'
    # cell vectors
    cell= atoms.get_cell()
    a = np.linalg.norm(cell[0,:])
    b = np.linalg.norm(cell[1,:])
    c = np.linalg.norm(cell[2,:])
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
    spos = atoms.get_scaled_positions()
    vels = atoms.get_velocities()
    if np.size(vels) != 3*len(atoms):
        vels = np.zeros((len(atoms),3))
    if not specorder:
        specorder = uniq(atoms.get_chemical_symbols())
        specorder.sort()
    for i in range(len(atoms)):
        atom= atoms[i]
        ifmv = ifmvs[i]
        txt += ' {0:s}'.format(get_tag(specorder,atom.symbol,i+1,ifmv))
        #...Scaled positions
        txt += ' {0:23.14e} {1:23.14e} {2:23.14e}'.format(spos[i,0],
                                                          spos[i,1],
                                                          spos[i,2])
        #...Scaled velocities
        txt += ' {0:15.7e} {1:15.7e} {2:15.7e}'.format(vels[i,0]/a,
                                                       vels[i,1]/b,
                                                       vels[i,2]/c)
        txt += '  0.0  0.0'
        txt += '  0.0  0.0  0.0  0.0  0.0  0.0\n'
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
            #...If the mask does not exist in the constraints list, add it
            constraints.append(mask)
            ifmvs[cnst.a] = len(constraints)
    
    #...Convert constraints to fmv
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

