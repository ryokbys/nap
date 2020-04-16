#!/usr/bin/env python
"""
Convert cell from the original one.

Usage:
  cell_convertor.py [options] INFILE OUTFILE

Options:
  -h,--help  Show this message and exit.
  --specorder=SPECORDER
             Set species order. [default: None]
"""
from __future__ import print_function

from docopt import docopt
import numpy as np

from nappy.napsys import NAPSystem

__author__ = 'Ryo KOBAYASHI'
__version__ = '190524'

def monocli_to_ortho(nsys):
    """
    Convert monoclinic cell to an orthogonal cell.
    """
    #...Determine monoclinic axis
    a0 = np.zeros((3,3),dtype=float)
    a0[0] = nsys.a1 *nsys.alc
    a0[1] = nsys.a2 *nsys.alc
    a0[2] = nsys.a3 *nsys.alc
    la0 = np.zeros(3)
    la0[0] = np.linalg.norm(a0[0])
    la0[1] = np.linalg.norm(a0[1])
    la0[2] = np.linalg.norm(a0[2])
    aa0 = np.zeros(3)
    aa0[0] = np.arccos(np.dot(a0[1],a0[2])/la0[1]/la0[2])/np.pi*180.0
    aa0[1] = np.arccos(np.dot(a0[0],a0[2])/la0[0]/la0[2])/np.pi*180.0
    aa0[2] = np.arccos(np.dot(a0[0],a0[1])/la0[0]/la0[1])/np.pi*180.0
    axis = -1
    if abs(aa0[2]-90.0) > 5.0:
        axis = 2
    elif abs(aa0[1]-90.0) > 5.0:
        axis = 1
    elif abs(aa0[0]-90.0) > 5.0:
        axis = 0
    if axis < 0:
        raise ValueError('The system seems to be orthogonal...')

    #...Determine how many multiplication to each axes and get new cell vectors, a1
    a1 = np.zeros((3,3),dtype=float)
    a1[axis] = a0[axis]
    axis1 = (axis +1) % 3
    axis2 = (axis1+1) % 3
    a1[axis1] = a0[axis1]
    n2 = round(abs(la0[axis1])/abs(la0[axis2]*np.cos(aa0[axis]/180.0*np.pi)))
    a1[axis2] = a0[axis2]*n2 -a1[axis1]*np.sign(np.dot(a0[axis1],a0[axis2]))

    #...Wrap atoms that are outside the new cell
    nrep = np.zeros((3),dtype=int)
    nrep[axis] = 1
    nrep[axis1] = 1
    nrep[axis2] = n2
    nsys.repeat(*nrep)
    rpos = nsys.get_real_positions()
    hmat1 = np.zeros((3,3))
    hmat1[:,0] = a1[0]
    hmat1[:,1] = a1[1]
    hmat1[:,2] = a1[2]
    hmat1i = np.linalg.inv(hmat1)
    sposs = np.zeros((nsys.num_atoms(),3))
    for i,ri in enumerate(rpos):
        sposs[i] = np.dot(hmat1i,ri)
    nsys.set_scaled_positions(sposs)
    nsys.set_hmat(hmat1)
    return nsys


if __name__ == '__main__':
    
    args = docopt(__doc__,version=__version__)
    infile = args['INFILE']
    outfile = args['OUTFILE']
    specorder = args['--specorder'].split(',')

    if specorder[0] == 'None':
        nsys = NAPSystem(fname=infile)
    else:
        nsys = NAPSystem(fname=infile,specorder=specorder)

    newsys = monocli_to_ortho(nsys)
    newsys.write(outfile)
