#!/usr/bin/env python
"""
Convert Quantum Espresso output to fitpot data.

Usage:
  qeout2fp.py [options] FILES [FILES...]

Options:
  -h, --help  Show this message and help.
  --specorder=SPECORDER
             Specify the order of species needed to convert to pos. [default: Al,Mg,Si]
  --index=INDEX
             Convert a snapshot of INDEX. [default: -1]

"""
from __future__ import print_function

import os,sys
from glob import glob
from docopt import docopt
import numpy as np

sys.path.append(os.environ['HOME']+'/src/nap/nappy')
from atom import Atom
from pmdsys import PMDSystem,unitvec_to_hi,cartessian_to_scaled

__author__ = "Ryo KOBAYASHI"
__version__ = "160620"

def get_tag(symbol,atom_id):
    sid= _specorder.index(symbol)+1
    tag= float(sid) +0.1 +atom_id*1e-14
    return '{0:16.14f}'.format(tag)

def write_pos(atoms,fname="pos"):
    cell= atoms.cell
    pos= atoms.get_scaled_positions()
    with open(fname,'w') as f:
        f.write('   1.000  \n')
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(cell[0,0],cell[0,1],cell[0,2]))
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(cell[1,0],cell[1,1],cell[1,2]))
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(cell[2,0],cell[2,1],cell[2,2]))
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' {0:10d}\n'.format(len(atoms)))
        for i in range(len(atoms)):
            atom= atoms[i]
            f.write(' {0:s}'.format(get_tag(atom.symbol,i+1)))
            f.write(' {0:12.8f} {1:12.8f} {2:12.8f}'.format(pos[i,0],pos[i,1],pos[i,2]))
            f.write(' 0.0 0.0 0.0 ')
            f.write(' 0.0 0.0 ' 
                    +' 0.0 0.0 0.0 0.0 0.0 0.0\n')

def write_frcref(fname='frc.ref',frcs=[]):
    with open(fname,'w') as f:
        f.write('{0:8d}\n'.format(len(frcs)))
        for i in range(len(frcs)):
            f.write('{0:24.14e}  {1:24.14e}  {2:24.14e}\n'.format(frcs[i,0],frcs[i,1],frcs[i,2]))


def write_ergref(fname='erg.ref',erg=0.0):
    with open(fname,'w') as f:
        f.write('{0:24.14e}\n'.format(erg))


def read_espresso_out(fname,):
    """
    Read cell info, atom coordinates, energy, forces from Quantum Espresso output.
    """
    
    f= open(fname,'r')
    lines = f.readlines()
    f.close()
    
    for il in range(len(lines)):
        if '   atomic species ' in lines[il]:
            il_spcs = il
        elif '!    total energy' in lines[il]:
            erg = float(lines[il].split()[4])
        elif 'Forces acting on atoms' in lines[il]:
            il_forces= il
        elif 'CELL_PARAMETERS' in lines[il]:
            il_cell = il
        elif 'ATOMIC_POSITIONS' in lines[il]:
            il_pos = il
        elif 'number of atoms/cell' in lines[il]:
            natm = int(lines[il].split()[4])
        elif 'number of atomic types' in lines[il]:
            nspcs = int(lines[il].split()[5])
            
    
    # read atom species
    il= il_spcs +1
    spcs = []
    for isp in range(nspcs):
        l = lines[il+isp].split()
        spcs.append(l[0])
    
    # read cell parameters
    il= il_cell
    cell = np.zeros((3,3),dtype=float)
    ixyz = 0
    while True:
        il += 1
        l = lines[il].split()
        if len(l) == 0:
            break
        cell[:,ixyz] = [ float(x) for x in l ]
        ixyz += 1
    
    # read atom coordinates
    elems = []
    pos = np.zeros((natm,3),dtype=float)
    il = il_pos + 1
    for ia in range(natm):
        l = lines[il+ia].split()
        elems.append(l[0])
        pos[ia,:] = [ float(x) for x in l[1:4] ]
    
    # read forces
    il = il_forces +2
    frcs = np.zeros((natm,3),dtype=float)
    for ia in range(natm):
        l = lines[il+ia].split()
        frcs[ia,:] = [ float(x) for x in l[6:9] ]

    return natm,nspcs,spcs,cell,pos,elems,erg,frcs


def convert(fname,specorder,index):
    """
    Convert data in fname to fitpot format.
    """

    ry2ev = 13.605698066
    au2ang = 0.529177249
    ryau2evang = ry2ev/au2ang

    if not os.path.exists(fname):
        raise RuntimeError('No '+fname+' exists.')
        
    #atoms= read('POSCAR',index=0,format='vasp')
    try:
        natm,nspcs,spcs,cell,pos,elems,erg,frcs = read_espresso_out(fname)
    except:
        raise IOError('Failed to read '+fname)

    psys = PMDSystem(specorder=spcs)
    psys.set_hmat(cell)
    hi = unitvec_to_hi(cell[0,:],cell[1,:],cell[2,:])
    frcs[:,:] *= ryau2evang
    for ia in range(natm):
        ai = Atom()
        pi = pos[ia,:]
        sx,sy,sz = cartessian_to_scaled(hi,pi[0],pi[1],pi[2])
        ai.set_pos(sx,sy,sz)
        ai.set_frc(frcs[ia,0],frcs[ia,1],frcs[ia,2]) # Ry/au at the moment
        ai.set_symbol(elems[ia])
        psys.add_atom(ai)
    psys.assign_pbc()
    psys.write_POSCAR()
    psys.write_pmd(fname='pos')
    write_ergref(fname='erg.ref',erg=erg)
    write_frcref(fname='frc.ref',frcs=frcs)


if __name__ == "__main__":
    
    args = docopt(__doc__)
    files = args['FILES']
    specorder= args['--specorder'].split(',')
    index= int(args['--index'])

    print('specorder = ',specorder)
    print('index = ',index)

    nfiles = len(files)
    print('number of files = ',nfiles)
    cwd = os.getcwd()
    for i in range(len(files)):
        f = files[i]
        path = os.path.dirname(f)
        fname = os.path.basename(f)
        print('{0:5d}/{1:d}: '.format(i+1,nfiles)+fname)
        if path == '':
            convert(fname,specorder,index)
        else:
            os.chdir(path)
            convert(fname,specorder,index)
        os.chdir(cwd)


