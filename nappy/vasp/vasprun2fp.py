#!/usr/bin/env python
"""
Convert from vasprun.xml in DIR direcotry
to erg.ref, frc.ref, and pos files.
Note that reading a lot of vasprun.xml takes a lot of time.

Usage:
  vasprun2fp.py [options]

Options:
  -h,--help  Show this message and exit.
  --specorder=SPECORDER
             Specify the order of species needed to convert POSCAR to pos. [default: None]
  --index=INDEX
             Convert a snapshot of INDEX. Comma separated indices can be specified. 
             If three digits are separated by colons like 0:1000:10, indices of slices 
             of initial(0), final(1000) and skip(10) are spcified. [default: -1]
  --sequence
             Extract all the sequence of MD or relaxation stored in vasprun.xml.
             If the index is specified as a list of indices, this option will be omitted.
  --velocity
             Compute velocities of atoms at time t from the diference of positions between t and t+dt. 
             This option is obsolete and it is not available in this version. [default: False]
"""
import os
from ase.io import read,write
from docopt import docopt
import numpy as np
import nappy

__author__ = "Ryo KOBAYASHI"
__version__ = "230917"

_kb2gpa = 160.2176487

def get_tag(symbol,atom_id,specorder):
    sid= specorder.index(symbol)+1
    tag= float(sid) +0.1 +atom_id*1e-14
    return '{0:17.14f}'.format(tag)

def output_for_fitpot(nsys, fname='smpl',
                      specorder=[]):
    try:
        epot = nsys.get_potential_energy()
    except:
        print(' Failed to get_potential_energy(), so skip it.')
        return None
    try:
        strs = nsys.get_stress_tensor()*_kb2gpa # Convert from kBar to GPa
        strs = strs.reshape(9)[[0, 4, 8, 5, 2, 1]]
    except:
        with open('WARNING','w') as f:
            f.write(' Since failed to get stress tensor, put 0.0s into strs.ref file.\n')
        strs = np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    nappy.io.write_pmd(nsys, fname=fname,
                       auxs=['fx','fy','fz'])
    return None

def main():
    args=docopt(__doc__)
    # dirs= args['DIR']
    specorder= args['--specorder'].split(',')
    sequence = args['--sequence']
    velocity = args['--velocity']
    if velocity:
        raise ValueError('velocity option is obsolete and not available in this version.')

    if specorder == 'None':
        raise ValueError('specorder must be specified.')
    print(' specorder = ',specorder)
    
    index= args['--index']
    if ',' in index:
        index = [ int(x) for x in index.split(',') ]
    elif ':' in index:
        index = [ int(x) for x in index.split(':') ]
        index = slice(*index)
    else:
        index = int(index)

    if type(index) is list:
        print(' The following steps are to be extracted: ',end='')
        for i in index:
            print(i,end='')
        print('')
        ase_index = ':'
    elif type(index) is slice:
        print(' The sliced steps are to be extracted, '+args['--index'])
        ase_index = index
    elif sequence:
        print(' All the sequence are to be extracted.')
        ase_index = ':'
    else:
        ase_index = index
        print(' index   = ',ase_index)

    if not os.path.exists('vasprun.xml'):
        raise FileNotFoundError(' No vasprun.xml found!')
    if os.path.exists('erg.ref') and \
       os.stat('erg.ref').st_mtime > os.stat('vasprun.xml').st_mtime:
        raise ValueError(' There is newer erg.ref than vasprun.xml.')
    try:
        #...Since there is a bug in vasp, species "r" needs to be replaced by "Zr"
        sysname, nodename, release, version, machine = os.uname()
        if 'Darwin' in sysname:
            os.system("sed -i '' -e 's|<c>r </c>|<c>Zr</c>|g' vasprun.xml")
        else:
            os.system("sed -i -e 's|<c>r </c>|<c>Zr</c>|g' vasprun.xml")
        #atoms= read('vasprun.xml',index=ase_index,format='vasp-xml')
        #...read_vasprun_xml always returns a list object
        nsyss = nappy.io.read_vasprun_xml(fname='vasprun.xml', velocity=velocity)
    except Exception as e:
        print(f' Failed to read vasprun.xml because of {e}.')
        raise

    if type(index) is list:
        print(' Extracting specified steps from ',len(nsyss),' steps in total')
        n = 0
        for j,nsys in enumerate(nsyss):
            if j not in index:
                continue
            fname = f'smpl_{j:05d}'
            output_for_fitpot(nsys,fname=fname,
                              specorder=specorder)
            n += 1
    elif sequence:  # Whole MD sequence
        print(' Extracting sequence of ',len(nsyss),' steps')
        indices = []
        # for j,nsys in enumerate(nsyss):
        for j,nsys in enumerate(nsyss):
            fname = f'smpl_{j:05d}'
            output_for_fitpot(nsys, fname=fname,
                              specorder=specorder)
    elif type(index) is slice:
        indices = index.indices(len(nsyss))
        print(' Extracting indiecs by the slice:', indices)
        inc = 0
        for nsys in nsyss[indices]:
            fname = f'smpl_{inc:05d}'
            output_for_fitpot(nsys, fname=fname,
                              specorder=specorder)
            inc += 1
    else:   # snapshot
        if index < 0:
            fname = f'smpl_{len(nsyss)+index}'
        else:
            fname = f'pmd_{index}'
        output_for_fitpot(nsyss[index],
                          fname=fname,
                          specorder=specorder)
    return None

if __name__ == "__main__":

    main()
