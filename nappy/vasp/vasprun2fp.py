#!/usr/bin/env python
"""
Convert from vasprun.xml in DIR direcotry
to erg.ref, frc.ref, pos, and POSCAR files.
Note that reading a lot of vasprun.xml takes a lot of time.

Usage:
  vasprun2fp.py [options] DIR [DIR...]

Options:
  -h,--help  Show this message and exit.
  --specorder=SPECORDER
             Specify the order of species needed to convert POSCAR to pos. [default: Al,Mg,Si]
  --index=INDEX
             Convert a snapshot of INDEX. [default: -1]
  --remove-constraints
             Remove constraints originally set to the system. [default: False]
"""

import os,sys
from ase.io import read,write
from glob import glob
from docopt import docopt

__author__ = "Ryo KOBAYASHI"
__version__ = "160507"

_specorder = []
_kb2gpa = 160.2176487

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



if __name__ == "__main__":

    args=docopt(__doc__)
    dirs= args['DIR']
    specorder= args['--specorder']
    index= int(args['--index'])
    remove_const = args['--remove-constraints']

    _specorder = specorder.split(',')
    print 'specorder = ',_specorder
    print 'index   = ',index
    print 'remove_const   = ',remove_const

    ndirs= len(dirs)
    print 'number of directories = ',ndirs

    cwd=os.getcwd()
    for i,d in enumerate(dirs):
        os.chdir(cwd)
        print '{0:5d}/{1:d}: '.format(i+1,ndirs)+d
        os.chdir(d)
        if not os.path.exists('vasprun.xml'):
            print 'No vasprun so skip.'
            continue
        if os.path.exists('erg.ref') and \
           os.stat('erg.ref').st_mtime > os.stat('vasprun.xml').st_mtime:
            print 'Since there is newer erg.ref, skip it.'
            continue
        try:
            atoms= read('vasprun.xml',index=index,format='vasp-xml')
        except:
            print 'Failed to read vasprun.xml, so skip it.'
            continue
        if remove_const:
            del atoms.constraints
        write('POSCAR',images=atoms,format='vasp',direct=True,vasp5=True)
        with open('erg.ref','w') as f:
            f.write("{0:12.7f}\n".format(atoms.get_potential_energy()))
        with open('frc.ref','w') as f:
            f.write("{0:6d}\n".format(len(atoms)))
            frcs= atoms.get_forces()
            for frc in frcs:
                f.write("{0:12.7f} {1:12.7f} {2:12.7f}\n".format(frc[0],frc[1],frc[2]))
        write_pos(atoms,fname='pos')
        with open('strs.ref','w') as f:
            strs = atoms.get_stress()
            for s in strs:
                f.write(" {0:15.7f}".format(s*_kb2gpa)) # converting from kBar to GPa
            f.write('\n')
    os.chdir(cwd)

