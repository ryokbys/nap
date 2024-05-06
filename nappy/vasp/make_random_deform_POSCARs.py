#!/usr/bin/env python
"""
Make randomly deformed POSCAR files from the non-deformed POSCAR file.

Usage:
  make_random_deform_POSCARs.py [options] POSCAR

Options:
  -h, --help  Show this help message and exit.
  -n NUM      Number of output POSCARs to be created. [default: 100]
  --deform-range=DEFORM_RANGE
              Range of the lattice deformation ratio. [default: 0.01,0.2]
  --displace-range=DISPLACE_RANGE
              Range of the atom displacements in ratio to (atomic volume)^{-3}. [default: 0.1]
  --num-displace=NUM_DISPLACE
              Number of displacements per deformation. [default: 2]
"""
from docopt import docopt
import numpy as np
import random
from ase.io import read,write

#======================================== subroutines and functions
def rotate(vec,axis,angle):
    """
    Rotate the vector `vec` by `angle` around `axis`.
    `vec` and `axis` are 3-D vector, `angle` is in degree.
    """
    theta = angle/180. *np.pi
    n = np.array(axis)
    n = n/np.linalg.norm(n)
    newvec = np.zeros(3,dtype=float)
    newvec = n *np.dot(n,vec) \
             +(vec- n*np.dot(n,vec))*np.cos(theta) \
             +np.cross(vec,n)*np.sin(theta)
    return newvec

def deformed_vector(a,maxlen):
    if type(a) is not np.ndarray:
        a = np.array(a)
    
    r = np.array([0.5-random.random(), 0.5-random.random(), 0.5-random.random()])
    #...make random vector normal to the vector `a`
    axr = np.cross(a,r)
    #...n is the axis of rotation
    n = axr / np.linalg.norm(axr)
    #...rotation angle between 0 to 180
    angle = 180.0 *random.random()
    #...length of deformation vector
    length = maxlen *random.random()
    da = a /np.linalg.norm(a) *length
    da = rotate(da,n,angle)
    return a+da


############################################################ main

if __name__ == "__main__":

    args = docopt(__doc__)

    num_data = int(args['-n'])
    deform_range = [ float(x) for x in args['--deform-range'].split(',') ]
    displace_range = float(args['--displace-range'])
    num_displace = int(args['--num-displace'])
    infname= args['POSCAR']

    atoms0 = read(infname,format='vasp')

    dlatd = (deform_range[1] -deform_range[0])/(num_data-1)
    ddisp = displace_range/num_displace
    print('deformation range = ',deform_range[0],deform_range[1])
    print('displace range    = ',displace_range)
    print('num of displace   = ',num_displace)
    print('dlatd = ',dlatd)
    print('ddisp = ',ddisp)

    n = 0
    while True:
        latd = dlatd*n +deform_range[0]
        cell0 = atoms0.get_cell()
        a = np.array(cell0[0])
        b = np.array(cell0[1])
        c = np.array(cell0[2])
        anew = deformed_vector(a,np.linalg.norm(a)*latd)
        bnew = deformed_vector(b,np.linalg.norm(b)*latd)
        cnew = deformed_vector(c,np.linalg.norm(c)*latd)
        cell = [anew,bnew,cnew]
        for idisp in range(num_displace):
            n += 1
            if n > num_data:
                break
            disp = ddisp*(idisp+1)
            atoms = atoms0.copy()
            atoms.set_cell(cell,scale_atoms=True)
            avol = atoms.get_volume()/len(atoms)
            avol3= avol**(1.0/3)
            dmax = avol3*disp
            pos = atoms.get_positions()
            for ia in range(1,len(atoms)):
                for l in range(3):
                    pos[ia,l] = pos[ia,l] +dmax*(2.0*random.random()-1.0)
            atoms.set_positions(pos)
            outfname = 'POSCAR_{0:05d}'.format(n)
            write(outfname,images=atoms,format='vasp',
                  direct=True,vasp5=True,sort=True)
            print('fname,avol= {0:s} {1:7.2f}'.format(outfname,
                                                      atoms.get_volume()/len(atoms)))
        if n > num_data:
            break
