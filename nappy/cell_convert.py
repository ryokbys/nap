#!/usr/bin/env python
"""
Change unit vectors of the system from the original a1, a2, and a3
to ax1, ax2, and ax3, each of them is written as (x1,x2,x3) that are
multiplying factors to a1, a2, and a3, respectively.
For example, if ax1 = (x1,x2,x3), the new a1 vector is given as:
  a1_new = x1*a1 +x2*a2 +x3*a3

Usage:
  cell_convert.py [options] INFILE

Options:
  -h,--help  Show this message and exit.
  --ortho    Convert to a nearly orthogonal cell where a1->x, a2->y, a3~z.
             Requires gamma=90 deg (a1 perp a2). Overrides --ax1/2/3.
  --ax1 AX1  Factors to be multiplied to a1,a2,a3, comma separated. [default: 1.0,0.0,0.0]
  --ax2 AX2  Factors to be multiplied to a1,a2,a3, comma separated. [default: 0.0,1.0,0.0]
  --ax3 AX3  Factors to be multiplied to a1,a2,a3, comma separated. [default: 0.0,0.0,1.0]
"""
import os,sys,copy
from docopt import docopt
import numpy as np

#sys.path.append(__file__)
import nappy
from nappy.napsys import NAPSystem

__author__ = 'Ryo KOBAYASHI'
__version__ = '250521'

def convert_by_multiply(nsys,ax1,ax2,ax3):
    """Convert the input nsys by using ax1,ax2,ax3, each of them is written as (x1,x2,x3) that are
multiplying factors to a1, a2, and a3, respectively.
For example, if ax1 = (x1,x2,x3), the new a1 vector is given as:
  a1_new = x1*a1 +x2*a2 +x3*a3
    """
    nsys.assign_pbc()
    #...Reset alc to 1.0
    nsys.a1 = nsys.a1 *nsys.alc
    nsys.a2 = nsys.a2 *nsys.alc
    nsys.a3 = nsys.a3 *nsys.alc
    nsys.alc = 1.0
    print(' a1  = ',nsys.a1)
    print(' a2  = ',nsys.a2)
    print(' a3  = ',nsys.a3)
    
    pos = nsys.get_real_positions()
    spos = nsys.get_scaled_positions()
    symbols = nsys.get_symbols()
    for i in range(min(len(nsys.atoms),10)):
        print(' {0:5d} {1:s}'.format(i,symbols[i])
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(spos[i,0],
                                                       spos[i,1],
                                                       spos[i,2])
              +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(pos[i,0],
                                                       pos[i,1],
                                                       pos[i,2]))

    
    hmat = nsys.get_hmat()
    hmati= nsys.get_hmat_inv()
    a1new = ax1[0]*nsys.a1 +ax1[1]*nsys.a2 +ax1[2]*nsys.a3
    a2new = ax2[0]*nsys.a1 +ax2[1]*nsys.a2 +ax2[2]*nsys.a3
    a3new = ax3[0]*nsys.a1 +ax3[1]*nsys.a2 +ax3[2]*nsys.a3
    #...Check if a1.(a2xa3) > 0
    if np.dot(a1new,np.cross(a2new,a3new)) < 0.0:
        raise ValueError('a1.(a2xa3) < 0 of new axes, which is now allowed in VASP.')
    sa1new = np.dot(hmati,a1new)
    sa2new = np.dot(hmati,a2new)
    sa3new = np.dot(hmati,a3new)
    print(' new a1 =',a1new)
    print(' new a2 =',a2new)
    print(' new a3 =',a3new)
    psnew = NAPSystem(specorder=nsys.specorder)
    psnew.set_lattice(nsys.alc,a1new,a2new,a3new)

    # Expand the original system for the search of atoms to be included 
    # in the new system.
    # First, compute how much we have to expand the original system
    # Compute the expansion range by summing positive/negative components
    # of new cell vectors in each old direction. Using max/min of individual
    # components is insufficient when multiple new vectors share a direction.
    sa_vecs = [sa1new, sa2new, sa3new]
    irange_lo = [0, 0, 0]
    irange_hi = [0, 0, 0]
    for i in range(3):
        comps = [v[i] for v in sa_vecs]
        pos_sum = sum(c for c in comps if c > 0.0)
        neg_sum = sum(c for c in comps if c < 0.0)
        irange_lo[i] = int(np.floor(neg_sum))
        irange_hi[i] = int(np.floor(pos_sum)) + 1  # +1 because range() is exclusive

    print(' irange_lo: ',irange_lo)
    print(' irange_hi: ',irange_hi)
    expos = []
    symbols = nsys.get_symbols()
    print(' symbols :',symbols)
    exsymbols = []
    print(' Expanding the original system...')
    for n3 in range(irange_lo[2],irange_hi[2]):
        for n2 in range(irange_lo[1],irange_hi[1]):
            for n1 in range(irange_lo[0],irange_hi[0]):
                for ia in range(len(spos)):
                    sposi = copy.deepcopy(spos[ia])
                    sposi[0] += n1
                    sposi[1] += n2
                    sposi[2] += n3
                    posi = np.dot(hmat,sposi)
                    symbol = symbols[ia]
                    # print(ia,n1,n2,n3,symbol,sposi)
                    expos.append(posi)
                    exsymbols.append(symbol)

    print(' Extracting the atoms inside the new unit vectors...')
    hmat= psnew.get_hmat()
    hi = np.linalg.inv(hmat)
    for ia,posi in enumerate(expos):
        sposi = np.dot(hi,posi)
        if 0.0 <= sposi[0] < 1.0 and \
           0.0 <= sposi[1] < 1.0 and \
           0.0 <= sposi[2] < 1.0:
            symbol = exsymbols[ia]
            # print(' {0:5d} {1:s}'.format(ia,symbol)
            #       +' {0:12.5f} {1:12.5f} {2:12.5f}'.format(sposi[0],
            #                                                sposi[1],
            #                                                sposi[2]))
            
            # atom.set_symbol(symbol)
            # atom.set_pos(sposi[0],sposi[1],sposi[2])
            # psnew.add_atom(atom)
            psnew.add_atoms([symbol],[sposi])
            
    # tmp = None
    # #tmp = raw_input('Input periodic shift vector if you want: ')
    # tmp = ' 0.5, 0.0, 0.5'
    # if tmp:
    #     shift = [ float(x) for x in tmp.split(',')]
    #     for a in psnew.atoms:
    #         a.pos[0] += shift[0]
    #         a.pos[1] += shift[1]
    #         a.pos[2] += shift[2]
    #     psnew.assign_pbc()
    # psnew.write_POSCAR(infile+'.new')
    # print('Check '+infile+'.new')

    psnew.assign_pbc()

    #...Remove atoms that are too close each other
    psnew.make_pair_list(rcut=1.0)
    poss = psnew.get_scaled_positions()
    hmat = psnew.get_hmat()
    to_remove = []
    for i in range(len(psnew)):
        pi = poss[i]
        pinorm = np.dot(pi,pi)
        for j in psnew.neighbors_of(i):
            if i in to_remove or j in to_remove:
                continue
            pj = poss[j]
            rij = pj-pi -np.round(pj-pi)
            rij_real = np.dot(hmat, rij)
            dij2 = np.dot(rij_real, rij_real)
            #...dij2 < 0.1 Ang^2 (~0.316 Ang) is too close, remove one of them
            if dij2 < 0.1:
                pjnorm = np.dot(pj,pj)
                if pjnorm < pinorm:
                    to_remove.append(j)
                else:
                    to_remove.append(i)

    psnew.remove_atoms(*to_remove)
                
    return psnew


def convert_to_ortho(nsys):
    """Convert to a nearly orthogonal cell where a1->x, a2->y, a3~z.
Requires gamma=90 deg (a1 perp a2).

Steps:
  1. Build rotation R: a1_dir->x, a2_dir->y, (a1xa2)_dir->z
  2. Express a3 in the rotated frame (a3r); flip sign if a3r[z] < 0
  3. Find integers n1, n2 minimizing the off-z components of a3r
  4. Apply R to the cell and (equivalently) atom positions, then call
     convert_by_multiply with ax3=(n1, n2, sign3)
    """
    nsys.assign_pbc()
    nsys.a1 = nsys.a1 *nsys.alc
    nsys.a2 = nsys.a2 *nsys.alc
    nsys.a3 = nsys.a3 *nsys.alc
    nsys.alc = 1.0

    a1 = nsys.a1.copy()
    a2 = nsys.a2.copy()
    a3 = nsys.a3.copy()

    la = np.linalg.norm(a1)
    lb = np.linalg.norm(a2)
    cos_gamma = np.dot(a1, a2) / (la * lb)
    if abs(cos_gamma) > 0.035:
        raise ValueError(f'gamma is not close to 90 deg (cos_gamma={cos_gamma:.4f}).'
                         ' convert_to_ortho requires gamma=90 deg.')

    # Rotation R: rows are the new x, y, z unit vectors expressed in old Cartesian.
    # R @ v rotates any Cartesian vector v into the new frame.
    e1 = a1 / la
    e2 = a2 / lb
    e3 = np.cross(e1, e2)
    e3 /= np.linalg.norm(e3)
    R = np.array([e1, e2, e3])

    a1r = R @ a1   # ~ (la, 0, 0)
    a2r = R @ a2   # ~ (0, lb, 0)
    a3r = R @ a3   # rotated a3 (z-component may be negative)

    print(' Rotation matrix R (rows = new x,y,z in old Cartesian):')
    for row in R:
        print('  [{:10.6f} {:10.6f} {:10.6f}]'.format(*row))
    print(' a1 (rotated) =', a1r)
    print(' a2 (rotated) =', a2r)
    print(' a3 (rotated) =', a3r)

    # Ensure z-component of a3 is positive for a right-handed system.
    sign3 = 1
    if a3r[2] < 0.0:
        a3r = -a3r
        sign3 = -1
    print(' sign3 =', sign3)

    # Integers n1, n2 that minimise the off-z (x, y) components of a3r.
    n1 = int(np.round(-a3r[0] / a1r[0]))
    n2 = int(np.round(-a3r[1] / a2r[1]))
    a3_new = a3r + n1*a1r + n2*a2r
    print(f' n1={n1}, n2={n2}  ->  a3_new = {a3_new}')

    vol = np.dot(a1r, np.cross(a2r, a3_new))
    if vol < 0.0:
        raise ValueError('Resulting cell is left-handed.')
    print(f' Volume = {vol:.4f} Ang^3')

    # Build a rotated NAPSystem: fractional positions are invariant under a
    # rotation that acts consistently on both cell vectors and atom positions,
    # so we reuse the original scaled positions with the rotated cell.
    nsys_rot = NAPSystem(specorder=nsys.specorder)
    nsys_rot.set_lattice(1.0, a1r, a2r, R @ a3)  # a3 before sign flip
    spos = nsys.get_scaled_positions()
    symbols = nsys.get_symbols()
    nsys_rot.add_atoms(symbols, spos)
    nsys_rot.assign_pbc()

    # Integer cell conversion in the rotated frame: new_a3 = n1*a1 + n2*a2 + sign3*a3
    ax3_v = [float(n1), float(n2), float(sign3)]
    print(' ax3 for convert_by_multiply =', ax3_v)
    return convert_by_multiply(nsys_rot, [1.0,0.0,0.0], [0.0,1.0,0.0], ax3_v)


def main():
    args = docopt(__doc__,version=__version__)
    infile = args['INFILE']
    nsys = nappy.io.read(infile)
    if args['--ortho']:
        newsys = convert_to_ortho(nsys)
    else:
        ax1 = [ float(x) for x in args['--ax1'].split(',') ]
        ax2 = [ float(x) for x in args['--ax2'].split(',') ]
        ax3 = [ float(x) for x in args['--ax3'].split(',') ]
        if len(ax1) != 3:
            raise ValueError('len(ax1) != 3')
        if len(ax2) != 3:
            raise ValueError('len(ax2) != 3')
        if len(ax3) != 3:
            raise ValueError('len(ax3) != 3')
        newsys = convert_by_multiply(nsys,ax1,ax2,ax3)
    nappy.io.write(newsys,infile+'_new')
    return None

if __name__ == '__main__':
    
    main()
