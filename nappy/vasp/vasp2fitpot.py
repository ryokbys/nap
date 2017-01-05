#!/usr/bin/env python
"""
Extract data necessary for *fitpot* program from vasprun.xml:
total energy, atom positions, forces on atoms, and cell information.
This uses *Vasprun* class in the *pymatgen* package.
Following files will be written:

* pos
* erg.ref
* frc.ref

If the vasprun.xml includes MD trajectory data, directories
with 5-digit name are created and above three files are
put in those directories.

Usage:
    vasp2fitpot.py [options] <vasprun.xml>

Options:
    -h, --help   show this help message and exit.
    -s, --skip=N
                 skip every N steps from the output of MD data. [default: 1]
    --idoffset=I
                 offset of species ID converted from vasprun.xml to pos. [default: 1]
"""

__author__    = "Ryo KOBAYASHI"
__email__     = "ryo.kbys@gmail.com"
__copyright__ = "Copyright 2015, Ryo KOBAYASHI"
__license__   = "MIT"
__version__   = "0.1"

_usage     = "vasp2fitpot.py [options] <vasprun.xml>"
_fname_pos = "pos"
_fname_erg = "erg.ref"
_fname_frc = "frc.ref"

import os,sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

import numpy as np
from docopt import docopt
from pymatgen.io.vaspio.vasp_output import Vasprun

from napsys import NAPSystem
from atom import Atom



def structure2aSys(structure,idoffset=1):
    """
    Converts Structure object of pymatgen to NAPSystem object in nap.

    Args:
        structure (Structure): pymatgen Structure object to be converted
        to NAPSystem object..

    Returns:
        aSys (NAPSystem): 
    """
    lattice= structure.lattice
    alc= 1.0
    a1= np.array(lattice.matrix[0])
    a2= np.array(lattice.matrix[1])
    a3= np.array(lattice.matrix[2])
    #... rescale a? vectors
    a1= a1/alc
    a2= a2/alc
    a3= a3/alc
    aSys= NAPSystem()
    aSys.set_lattice(alc,a1,a2,a3)
    for ia in range(structure.num_sites):
        ai= Atom()
        si= structure[ia]
        crd= si.frac_coords
        ai.set_pos(crd[0],crd[1],crd[2])
        sid= structure.symbol_set.index(si.species_string)+idoffset
        ai.set_sid(sid)
        ai.set_id(ia+1)
        aSys.add_atom(ai)
    return aSys

def write_ergref(fname,energy):
    f=open(fname,'w')
    f.write(' {0:20.8f}\n'.format(energy))
    f.close()
    
def write_frcref(fname,forces):
    f= open(fname,'w')
    f.write(' {0:10d}\n'.format(len(forces)))
    for fi in forces:
        f.write('  {0:12.7f}'.format(fi[0]))
        f.write('  {0:12.7f}'.format(fi[1]))
        f.write('  {0:12.7f}'.format(fi[2]))
        f.write('\n')
    f.close()

if __name__ == "__main__":

    args= docopt(__doc__)
    
    _infname= args['<vasprun.xml>']
    iskip= int(args['--skip'])
    idoffset= int(args['--idoffset'])

    #...parse vasprun.xml using Vasprun in pymatgen
    vasprun= Vasprun(_infname)
    
    ibrion= vasprun.incar['IBRION']

    if ibrion == 0: # MD
        nstps= len(vasprun.ionic_steps)
        for istp in range(0,nstps,iskip):
            dirname= '{0:05d}'.format(istp+1)
            os.mkdir(dirname)
            vaspstep= vasprun.ionic_steps[istp]
            aSys= structure2aSys(vaspstep['structure'],idoffset)
            energy= vaspstep['electronic_steps'][-1]['e_fr_energy']
            forces= vaspstep['forces']
            aSys.write_pmd(dirname+'/'+_fname_pos)
            write_ergref(dirname+'/'+_fname_erg,energy)
            write_frcref(dirname+'/'+_fname_frc,forces)
            print ' write ',dirname,_fname_pos,_fname_erg,_fname_frc
    else: # only final structure is needed
        aSys= structure2aSys(vasprun.final_structure,idoffset)
        energy= vasprun.ionic_steps[-1]['electronic_steps'][-1]['e_fr_energy']
        forces= vasprun.ionic_steps[-1]['forces']
        aSys.write_pmd(_fname_pos)
        write_ergref(_fname_erg,energy)
        write_frcref(_fname_frc,forces)
        print ' write ',_fname_pos,_fname_erg,_fname_frc
    
