#!/usr/bin/env python
"""
Convert from vasprun.xml in DIR direcotry
to erg.ref, frc.ref, and pos files.
Note that reading a lot of vasprun.xml takes a lot of time.

Usage:
  vasprun2fp.py [options] DIR [DIR...]

Options:
  -h,--help  Show this message and exit.
  --specorder=SPECORDER
             Specify the order of species needed to convert POSCAR to pos. [default: Al,Mg,Si]
  --index=INDEX
             Convert a snapshot of INDEX. [default: -1]
  --sequence
             Extract all the sequence of MD or relaxation stored in vasprun.xml.
  --keep-constraints
             Keep constraints originally set to the system. 
             Otherwise all the constratins are removed. [default: False]
"""

import os
from ase.io import read,write
from docopt import docopt

__author__ = "Ryo KOBAYASHI"
__version__ = "170620"

_kb2gpa = 160.2176487

def get_tag(symbol,atom_id,specorder):
    sid= specorder.index(symbol)+1
    tag= float(sid) +0.1 +atom_id*1e-14
    return '{0:16.14f}'.format(tag)

def write_pos(atoms,fname="pos",specorder=[]):
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
            f.write(' {0:s}'.format(get_tag(atom.symbol,i+1,specorder)))
            f.write(' {0:12.8f} {1:12.8f} {2:12.8f}'.format(pos[i,0],pos[i,1],pos[i,2]))
            f.write(' 0.0 0.0 0.0 ')
            f.write(' 0.0 0.0 ' 
                    +' 0.0 0.0 0.0 0.0 0.0 0.0\n')


def output_for_fitpot(atoms,keep_const,dirname='./',specorder=[]):
    if not keep_const:
        del atoms.constraints
    #write(dirname+'/POSCAR',images=atoms,format='vasp',direct=True,vasp5=True)
    try:
        epot = atoms.get_potential_energy()
    except:
        print ' Failed to get_potential_energy(), so skip it.'
        return None
    with open(dirname+'/erg.ref','w') as f:
        f.write("{0:12.7f}\n".format(epot))
    with open(dirname+'/frc.ref','w') as f:
        f.write("{0:6d}\n".format(len(atoms)))
        frcs= atoms.get_forces()
        for frc in frcs:
            f.write("{0:12.7f} {1:12.7f} {2:12.7f}\n".format(frc[0],frc[1],frc[2]))
    write_pos(atoms,fname=dirname+'/pos',specorder=specorder)
    if not os.path.exists(dirname+'/POSCAR'):
        write(dirname+'/POSCAR',images=atoms,format='vasp',
              direct=True,vasp5=True,sort=True)
    with open(dirname+'/strs.ref','w') as f:
        strs = atoms.get_stress()
        for s in strs:
            f.write(" {0:15.7f}".format(s*_kb2gpa)) # converting from kBar to GPa
        f.write('\n')
    

if __name__ == "__main__":

    args=docopt(__doc__)
    dirs= args['DIR']
    specorder= args['--specorder'].split(',')
    index= int(args['--index'])
    sequence = args['--sequence']
    keep_const = args['--keep-constraints']

    print ' specorder = ',specorder
    if sequence:
        print ' All the sequence are to be extracted.'
        index = ':'
    else:
        print ' index   = ',index
    print ' keep_const   = ',keep_const

    ndirs= len(dirs)
    print ' number of directories = ',ndirs

    cwd=os.getcwd()
    for i,d in enumerate(dirs):
        os.chdir(cwd)
        print '{0:5d}/{1:d}: '.format(i+1,ndirs)+d
        os.chdir(d)
        if not os.path.exists('vasprun.xml'):
            print ' No vasprun.xml, so skip.'
            continue
        if os.path.exists('erg.ref') and \
           os.stat('erg.ref').st_mtime > os.stat('vasprun.xml').st_mtime:
            print ' Since there is newer erg.ref, skip it.'
            continue
        try:
            #...Since there is a bug in vasp, species "r" needs to be replaced by "Zr"
            os.system("sed -i '' -e 's|<c>r </c>|<c>Zr</c>|g' vasprun.xml")
            atoms= read('vasprun.xml',index=index,format='vasp-xml')
        except Exception as e:
            print ' Failed to read vasprun.xml, so skip it.'
            print e
            continue

        if sequence:  # MD sequence
            print(' Extracting sequence...')
            for j,a in enumerate(atoms):
                dirname = '{0:05d}/'.format(j)
                print('  {0:s}'.format(dirname))
                os.system('mkdir -p {0:s}'.format(dirname))
                output_for_fitpot(a,keep_const,dirname=dirname,
                                  specorder=specorder)
            pass
        else:   # snapshopt
            dirname = './'
            output_for_fitpot(atoms,keep_const,dirname=dirname,
                              specorder=specorder)
    os.chdir(cwd)

