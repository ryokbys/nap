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
             Compute velocities of atoms at time t from the diference of positions between t and t+dt. [default: False]
  --force
             Write forces in pmdini files. [default: False]
"""
from __future__ import print_function

import os
from ase.io import read,write
from docopt import docopt
import numpy as np
import nappy

__author__ = "Ryo KOBAYASHI"
__version__ = "230917"

_kb2gpa = 160.2176487

def read_vasprun_xml(fname='vasprun.xml', velocity=False):
    """
    Based on read_vasp_xml in ase.io.vasp.py.
    """
    import xml.etree.ElementTree as ET
    
    tree = ET.iterparse(fname, events=['start', 'end'])
    calcs = []
    dt = -1.0
    specorder = []
    try:
        for event, elem in tree:
            if event == 'end':
                if elem.tag == 'atominfo':
                    species = []
                    for entry in elem.find("array[@name='atoms']/set"):
                        species.append(entry[0].text.strip())
                    natoms = len(species)
                    specorder = sorted(list(set(species)),key=species.index)
                elif elem.tag == 'incar':
                    for e in elem.iter():
                        if 'name' in e.attrib.keys() \
                          and e.attrib['name'] == 'POTIM':
                            dt = float(e.text)
            elif event == 'start' and elem.tag == 'calculation':
                calcs.append(elem)
    except ET.ParseError as e:
        raise e
    
    if len(calcs)==0:
        raise ValueError(f'There is no calculation in {fname}')
    
    nsyss = []
    for calc in calcs:
        nsys = nappy.napsys.NAPSystem(specorder=specorder)

        lastscf = calc.findall('scstep/energy')[-1]
        de = (float(lastscf.find('i[@name="e_0_energy"]').text) -
              float(lastscf.find('i[@name="e_fr_energy"]').text))
        e_free = float(calc.find('energy/i[@name="e_fr_energy"]').text)
        epot = e_free + de
        nsys.set_potential_energy(epot)
        
        cell = np.zeros((3, 3), dtype=float)
        for i, vector in enumerate(
                calc.find('structure/crystal/varray[@name="basis"]')):
            cell[i] = np.array([float(val) for val in vector.text.split()])

        sposs = np.zeros((natoms, 3), dtype=float)
        for i, vector in enumerate(
                calc.find('structure/varray[@name="positions"]')):
            sposs[i] = np.array([float(val) for val in vector.text.split()])

        frcs = None
        fblocks = calc.find('varray[@name="forces"]')
        if fblocks is not None:
            frcs = np.zeros((natoms, 3), dtype=float)
            for i, vector in enumerate(fblocks):
                frcs[i] = np.array(
                    [float(val) for val in vector.text.split()])

        strs = None
        sblocks = calc.find('varray[@name="stress"]')
        if sblocks is not None:
            strs = np.zeros((3, 3), dtype=float)
            for i, vector in enumerate(sblocks):
                strs[i] = np.array(
                    [float(val) for val in vector.text.split()])
            #strs *= -0.1 * GPa
            #strs = strss.reshape(9)[[0, 4, 8, 5, 2, 1]]
        
        nsys.set_hmat(cell.T) # hmat and cell are in transpose relation
        nsys.add_atoms(species, sposs, frcs=frcs)
        nsys.set_stress_tensor(strs)
        nsyss.append(nsys)
    
    if velocity:
        # Assuming that the change/velocity of the cell does not contribute to atom velocities,
        # which is a common assumption in Andersen and Parrinello-Rahman methods.
        # But this assumption may not be apropriate when the precise measurement of velocity is important.
        nsys0 = nsyss[0]
        vels = np.zeros((natoms,3), dtype=float)
        for n in range(1,len(nsyss)):
            nsys1 = nsyss[n]
            hmat0 = nsys0.get_hmat()
            sposs0 = nsys0.get_scaled_positions()
            sposs1 = nsys1.get_scaled_positions()
            for i in range(len(sposs0)):
                dpos = sposs1[i] -sposs0[i]
                dpos = dpos -np.round(dpos)
                vels[i] = np.dot(hmat0, dpos) /dt
            nsyss[n-1].set_real_velocities(vels)
            nsys0 = nsys1
        # Since velocities of the last step cannot be computed in this definition,
        # remove the last one for data consistency...
        del nsyss[-1]
    return nsyss

def get_tag(symbol,atom_id,specorder):
    sid= specorder.index(symbol)+1
    tag= float(sid) +0.1 +atom_id*1e-14
    return '{0:17.14f}'.format(tag)

def write_pmdini(nsys,fname="pmdini",velocity=False):
    hmat = nsys.get_hmat()
    a1,a2,a3 = nsys.get_lattice_vectors()
    specorder = nsys.specorder
    spos= nsys.get_scaled_positions()
    if velocity:
        vels = nsys.get_velocities()
    else:
        vels = np.zeros((len(nsys),3), dtype=float)

    with open(fname,'w') as f:
        f.write('!\n')
        f.write('!  specorder: ')
        for s in specorder:
            f.write(' {0:<3s}'.format(s))
        f.write('\n')
        f.write('!\n')
        f.write('   1.000  \n')
        # f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(hmat[0,0],hmat[0,1],hmat[0,2]))
        # f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(hmat[1,0],hmat[1,1],hmat[1,2]))
        # f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(hmat[2,0],hmat[2,1],hmat[2,2]))
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(*a1))
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(*a2))
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e}\n'.format(*a3))
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' 0.00000000 0.00000000 0.00000000\n')
        f.write(' {0:10d}\n'.format(len(nsys)))
        symbols = nsys.get_symbols()
        for i in range(len(nsys)):
            f.write(' {0:s}'.format(get_tag(symbols[i],i+1,specorder))
                    +' {0:12.8f} {1:12.8f} {2:12.8f}'.format(*spos[i])
                    +' {0:12.4e} {1:12.4e} {2:12.4e}'.format(*vels[i])
                    +'\n')
    return None

def output_for_fitpot(nsys, dirname='./', specorder=[],
                      velocity=False, force=False):
    try:
        epot = nsys.get_potential_energy()
    except:
        print(' Failed to get_potential_energy(), so skip it.')
        return None
    with open(dirname+'/erg.ref','w') as f:
        f.write("{0:12.7f}\n".format(epot))
    with open(dirname+'/frc.ref','w') as f:
        f.write("{0:6d}\n".format(len(nsys)))
        frcs= nsys.get_forces()
        for frc in frcs:
            f.write("{0:12.7f} {1:12.7f} {2:12.7f}\n".format(frc[0],frc[1],frc[2]))
    if not os.path.exists(dirname+'/POSCAR'):
        nappy.io.write_POSCAR(nsys, fname=dirname+'/POSCAR')
    try:
        strs = nsys.get_stress_tensor()
        strs = strs.reshape(9)[[0, 4, 8, 5, 2, 1]]
    except:
        with open(dirname+'/WARNING','w') as f:
            f.write(' Since failed to get stress tensor, put 0.0s into strs.ref file.\n')
        strs = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    with open(dirname+'/strs.ref','w') as f:
        for s in strs:
            # Convert from kBar to GPa
            # The minus sign comes from the difference of definition betw pmd.
            # f.write(" {0:15.7f}".format(s*(-_kb2gpa)))
            #...It seems no need to multiply (-1).
            f.write(" {0:15.7f}".format(s*_kb2gpa))
        f.write('\n')
    #write_pmdini(nsys,fname=dirname+'/pmdini')
    nappy.io.write_pmd(nsys, fname=dirname+'/pmdini',
                       potential_energy=epot,
                       stress=strs,
                       forces=force)
    return None

def main():
    args=docopt(__doc__)
    # dirs= args['DIR']
    specorder= args['--specorder'].split(',')
    sequence = args['--sequence']
    velocity = args['--velocity']
    force = args['--force']

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
        print(' The sliced steps are to be extracted: ')
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
        nsyss = read_vasprun_xml(fname='vasprun.xml', velocity=velocity)[index]
    except Exception as e:
        raise Exception(' Failed to read vasprun.xml because of {0}.'.format(e))

    if type(index) is list:
        print(' Extracting specified steps from ',len(nsyss),' steps in total')
        n = 0
        for j,nsys in enumerate(nsyss):
            if j not in index:
                continue
            dirname = '{0:05d}/'.format(n)
            print('  {0:s}'.format(dirname))
            os.system('mkdir -p {0:s}'.format(dirname))
            output_for_fitpot(nsys,dirname=dirname,
                              specorder=specorder,
                              velocity=velocity,
                              force=force)
            n += 1
    elif sequence or type(index) is slice:  # Whole MD sequence
        print(' Extracting sequence of ',len(nsyss),' steps')
        indices = []
        for j,nsys in enumerate(nsyss):
            dirname = '{0:05d}/'.format(j)
            print('  {0:s}'.format(dirname))
            os.system('mkdir -p {0:s}'.format(dirname))
            output_for_fitpot(nsys,dirname=dirname,
                              specorder=specorder,
                              velocity=velocity,
                              force=force)
        pass
    else:   # snapshopt
        dirname = './'
        output_for_fitpot(nsyss[0],dirname=dirname,
                          specorder=specorder,
                          velocity=velocity,
                          force=force)
    return None

if __name__ == "__main__":

    main()
