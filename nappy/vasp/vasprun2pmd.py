#!/usr/bin/env python
"""
Convert from vasprun.xml to pmdini file that contains the following information,
  - potential energy (eV)
  - kinetic energy (eV)
  - stress (GPa)
  - forces (eV/Ang.)

Usage:
  vasprun2pmd.py [options]

Options:
  -h,--help  Show this message and exit.
  --index=INDEX
             Convert a snapshot of INDEX. Comma separated indices can be specified. 
             If three digits are separated by colons like 0:1000:10, indices of slices 
             of initial(0), final(1000) and skip(10) are spcified. [default: -1]
  --sequence
             Extract all the sequence of MD or relaxation stored in vasprun.xml.
             If the index is specified as a list of indices, this option will be omitted.
  --velocity
             Compute velocities of atoms at time t from the diference of positions between t and t+dt. [default: False]
"""
import os
from docopt import docopt
import numpy as np
import nappy

__author__ = "Ryo KOBAYASHI"
__version__ = "231031"

_kb2gpa = 160.2176487

def output_pmdini(nsys, postfix='_',
                  velocity=False, force=False):
    try:
        epot = nsys.get_potential_energy()
    except:
        print(' Failed to get_potential_energy(), so skip it.')
        return None
    try:
        strs = nsys.get_stress_tensor()
        strs = strs.reshape(9)[[0, 4, 8, 5, 2, 1]]
    except:
        with open('WARNING','w') as f:
            f.write(' Since failed to get stress tensor,'
                    +' put 0.0s into strs.ref file.\n')
        strs = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #write_pmdini(nsys,fname=dirname+'/pmdini')
    nappy.io.write_pmd(nsys, fname='pmdini'+postfix,
                       potential_energy=epot,
                       stress=strs,
                       forces=force)
    return None

def main():
    args=docopt(__doc__)
    # dirs= args['DIR']
    sequence = args['--sequence']
    velocity = args['--velocity']

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
        raise FileNotFoundError(' No vasprun.xml found in the working directory!')
    try:
        #...Since there is a bug in vasp, species "r" needs to be replaced by "Zr"
        sysname, nodename, release, version, machine = os.uname()
        if 'Darwin' in sysname:
            os.system("sed -i '' -e 's|<c>r </c>|<c>Zr</c>|g' vasprun.xml")
        else:
            os.system("sed -i -e 's|<c>r </c>|<c>Zr</c>|g' vasprun.xml")
        #atoms= read('vasprun.xml',index=ase_index,format='vasp-xml')
        #...read_vasprun_xml always returns a list object
        nsyss = nappy.io.read_vasprun_xml(fname='vasprun.xml',
                                          velocity=velocity)
    except Exception as e:
        print(f' Failed to read vasprun.xml because of {e}.')
        raise

    if type(index) is list:
        print(' Extracting specified steps from ',len(nsyss),
              ' steps in total')
        n = 0
        for j,nsys in enumerate(nsyss):
            if j not in index:
                continue
            postfix = f'_{n:05d}'
            print(f'  {postfix}')
            # dirname = '{0:05d}/'.format(n)
            # print('  {0:s}'.format(dirname))
            # os.system('mkdir -p {0:s}'.format(dirname))
            output_pmdini(nsys, postfix=postfix,
                          velocity=velocity,
                          force=True)
            n += 1
    elif sequence:
        if type(index) is slice:
            nsyss_slc = nsyss[index]
            print(' Extracting sequence of ',len(nsyss_slc),' steps')
            for j,nsys in enumerate(nsyss_slc):
                postfix = f'_{j:05d}'
                output_pmdini(nsys, postfix=postfix,
                              velocity=velocity,
                              force=True)
        else: # Whole MD sequence
            print(' Extracting sequence of ',len(nsyss),' steps')
            for j,nsys in enumerate(nsyss):
                postfix = f'_{j:05d}'
                output_pmdini(nsys, postfix=postfix,
                              velocity=velocity,
                              force=True)
    else:   # snapshopt
        # dirname = './'
        postfix = ''
        output_pmdini(nsyss[0], postfix=postfix,
                      velocity=velocity,
                      force=True)
    return None

if __name__ == "__main__":

    main()
