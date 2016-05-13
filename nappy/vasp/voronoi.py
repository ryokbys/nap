#!/bin/env python
# -*- coding: utf-8 -*-
"""
Voronoi analysis based on Delauney tetrahedra.

USAGE:
    $ python ./voronoi.py [options] POSCAR

INPUT: (these files must be in the working directory)
    - POSCAR (for the cell information)
"""

import os,optparse
from atom_system import AtomSystem

def voronoi(ia,aSys):


#=======================================================================
if __name__ == '__main__':

    usage= '$ python %prog [options] POSCAR'

    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-r",dest="rcut",type="float", \
                      default=5.0, \
                      help="Cutoff radius in Angstrom. Default is 5.0.")
    (options,args)= parser.parse_args()

    rcut= options.rcut
    infname= args[0]
    
    aSys= AtomSystem()
    if not os.path.exists(infname):
        print '[Error] File does not exist !!!'
        print infname
        exit()
    aSys.read_POSCAR(infname)

    #.....make neighbor list
    aSys.make_pair_list(rcut)

    #.....make Voronoi polygons from Delauney tetrahedra
    vrns= []
    for ia in range(aSys.num_atoms):
        vrn.append(voronoi(ia,aSys))
        
