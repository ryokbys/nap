#!/usr/bin/env python
"""
Convert POSCAR to QuantumESPRESSO input.

Usage:
  poscar2qein.py [options] POSCAR

Options:
  -h, --help  Show this message and exit.
  --calc=CALC
              Calculation type. [default: scf]
  --constraints=CONSTRAINTS
              Constraint string added to each entry of atom.
              '1' indicates the atom moving, while '0' being fixed.
              [default: 1,1,1]
  -o OUTFNAME
              Output file name. [default: in.pw]
  -p,--pitch=PITCH
              Pitch of k-points. [default: 0.0968]
  -s SCALEFACTOR
              Scaling factor for the original lattice constant.
              [default: 1.0]
  --pseudo-dir PSUEDODIR
              Path to the directory where pseudopotential data exist,
              which is appended to $HOME environmental variable.
              [default: /local/espresso/SSSP_acc_PBESOL/]
"""
from __future__ import print_function

import os,sys
import glob
from docopt import docopt
from ase.io import read,write
import ase.data
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "170314"

qe_params = {
    '&CONTROL':{
        'calculation': 'vc-relax',
        'max_seconds': 85400,
        'outdir': './out/',
        'disk_io': 'none',
        'prefix': 'pw',
        'pseudo_dir': '/Users/kobayashi/local/espresso/SSSP_acc_PBESOL/',
        'restart_mode': 'from_scratch',
        'verbosity': 'high',
        'wf_collect': False,
        'tprnfor': True,
        'tstress': True,
    },
    '&SYSTEM':{
        'degauss': 0.05,
        'ecutrho': 280,
        'ecutwfc': 35,
        'occupations': 'smearing',
        'smearing': 'gauss',
        'nat': 0,
        'ntyp': 3,
        'ibrav': 0,
        'nbnd': 0,
    },
    '&ELECTRONS':{
        'conv_thr': 1.0e-04,
        'mixing_beta': 7.0e-01,
        'mixing_mode': 'plain',
        'diagonalization': 'cg',
    },
    '&IONS':{
        'ion_dynamics': 'bfgs',
    },
    '&CELL':{
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'cell_dofree': 'all',
    },
}

namelist_order = (
    '&CONTROL',
    '&SYSTEM',
    '&ELECTRONS',
    '&IONS',
    '&CELL',
)

nelectrons = {
    'Al': 3,
    'Mg': 10,
    'Si': 4,
}

atomic_species_txt = """ATOMIC_SPECIES
Al  26.981538 Al.pbe-n-rrkjus_psl.1.0.0.UPF
Mg  24.305    Mg.pbe-v1.4-uspp.F.UPF
Si  28.0855   Si.pbe-n-rrkjus_psl.1.0.0.UPF
"""


def override_qe_params(**kwargs):
    global qe_params
    
    for namelist in qe_params.keys():
        for k,v in qe_params[namelist].items():
            if k in kwargs.keys():
                qe_params[namelist][k] = kwargs[k]
    

def write_espresso_in(atoms,outfname,k_shift=(0.,0.,0.),pitch=0.0968,
                      symbols=[],ppfiles=[],valences=[],
                      constraints=[],
                      **kwargs):
    """
    Write input file for PWscf of QuantumESPRESSO.
    """
    global qe_params

    fo = open(outfname,'w')

    override_qe_params(**kwargs)
    psdir = kwargs['pseudo_dir']

    #...Some necessary paramters.
    elms = atoms.get_chemical_symbols()
    for namelist in namelist_order:
        ctx = qe_params[namelist]
        fo.write(namelist+'\n')
        for k,v in ctx.items():
            if k == 'nat':
                fo.write("   {0:s} = {1:d}\n".format(k,len(atoms)))
                continue
            elif k == 'nbnd':
                nbnd = 0
                for e in elms:
                    nbnd += valences[e]
                fo.write("   {0:s} = {1:d}\n".format(k,int(nbnd)))
                continue
            if isinstance(v,bool):
                if v:
                    fo.write("   {0:s} = .true.\n".format(k,v))
                else:
                    fo.write("   {0:s} = .false.\n".format(k,v))
            elif isinstance(v,str):
                fo.write("   {0:s} = '{1:s}'\n".format(k,v))
            elif isinstance(v,int):
                fo.write("   {0:s} = {1:d}\n".format(k,v))
            elif isinstance(v,float):
                fo.write("   {0:s} = {1:f}\n".format(k,v))
        fo.write('/\n')

    # if 'relax' in qe_params['&CONTROL']['calculation']:
    #     fo.write('&IONS\n')
    #     fo.write("   ion_dynamics = 'bfgs'\n")
    #     fo.write('/\n')
    # if 'vc' in qe_params['&CONTROL']['calculation']:
    #     fo.write('&CELL\n')
    #     fo.write("   cell_dynamics = 'bfgs'\n")
    #     fo.write('/\n')

    ####### INPUT_CARDS hereafter #######

    #...ATOMIC_SPECIES
    fo.write('ATOMIC_SPECIES\n')
    for i,s in enumerate(symbols):
        anum = ase.data.atomic_numbers[s]
        mass = ase.data.atomic_masses[anum]
        fo.write('{0:s}  {1:8.3f}  {2:s}\n'.format(s,mass,ppfiles[i]))
    #fo.write(atomic_species_txt)
    
    #...CELL_PARAMETERS
    fo.write('CELL_PARAMETERS angstrom\n')
    cell = atoms.get_cell()
    fo.write('{0:15.7f} {1:15.7f} {2:15.7f}\n'.format(cell[0,0],cell[0,1],cell[0,2]))
    fo.write('{0:15.7f} {1:15.7f} {2:15.7f}\n'.format(cell[1,0],cell[1,1],cell[1,2]))
    fo.write('{0:15.7f} {1:15.7f} {2:15.7f}\n'.format(cell[2,0],cell[2,1],cell[2,2]))

    #...ATOMIC_POSITIONS
    fo.write('ATOMIC_POSITIONS crystal\n')
    spos = atoms.get_scaled_positions()
    for e,sp in zip(elms,spos):
        fo.write('{0:s} '.format(e))
        fo.write('{0:15.7f} {1:15.7f} {2:15.7f} '.format(sp[0],sp[1],sp[2]))
        if constraints:
            fo.write(' {:s} '.format(constraints[0])+
                     ' {:s} '.format(constraints[1])+
                     ' {:s} '.format(constraints[2]))
        fo.write('\n')

    #...K_POINTS
    fo.write('K_POINTS automatic\n')
    k1,k2,k3 = get_kpoints(cell,pitch)
    fo.write('{0:4d} {1:4d} {2:4d} '.format(k1,k2,k3))
    fo.write('{0:d} {1:d} {2:d}\n'.format(0,0,0))  # k_shift, not implemented
    fo.close()


def get_kpoints(cell,pitch):
    a1 = cell[0]
    a2 = cell[1]
    a3 = cell[2]
    vol= np.dot(a1,np.cross(a2,a3))
    b1 = np.linalg.norm(np.cross(a2,a3)/vol)*2.0*np.pi 
    b2 = np.linalg.norm(np.cross(a3,a1)/vol)*2.0*np.pi
    b3 = np.linalg.norm(np.cross(a1,a2)/vol)*2.0*np.pi
    k1 = max(int(b1/pitch),1)
    k2 = max(int(b2/pitch),1)
    k3 = max(int(b3/pitch),1)
    return k1,k2,k3

def uniq(arr):
    uniqarr = []
    for a in arr:
        if not a in uniqarr:
            uniqarr.append(a)
    return uniqarr
    
def get_PP_files(symbols,psdir):
    files = [ os.path.basename(f) for f in glob.glob(psdir+'/*[Uu][Pp][Ff]') ]
    ppfiles = []
    for s in symbols:
        for f in files:
            if s.lower() in f[:4].lower():
                ppfiles.append(f)
                break
    if len(symbols) != len(ppfiles):
        print('symbols = ',symbols)
        print('ppfiles = ',ppfiles)
        print('files = ',files)
        raise ValueError('Something wrong...')
    return ppfiles

def get_valence(upf):
    valence = 0.0
    with open(upf,'r') as f:
        for line in f.readlines():
            if 'z_val' in line:
                data = line.split('=')[1]
                valence = float(data.replace('"',''))
                break
            elif 'Z val' in line:
                valence = float(line.split()[0])
                break
    return valence

if __name__ == "__main__":

    args = docopt(__doc__)
    poscar = args['POSCAR']
    outfname = args['-o']
    sfac = float(args['-s'])
    calc = args['--calc']
    constraints = args['--constraints'].split(',')
    pitch = float(args['--pitch'])
    psdir = args['--pseudo-dir']

    #...PP directory
    psdir = os.environ['HOME']+'/'+psdir
    
    atoms = read(poscar,format="vasp")
    print('system: {}'.format(atoms.get_chemical_formula()))

    #...some scaling to the system
    cell = atoms.get_cell()
    cell *= sfac
    atoms.set_cell(cell,scale_atoms=True)
    print('cell:')
    print(cell)

    #...get PP file
    symbols = uniq(atoms.get_chemical_symbols())
    ppfiles = get_PP_files(symbols,psdir)
    print('PP files:')
    for ppf in ppfiles:
        print(ppf)

    #...get_valence electrons
    valences = {}
    for s,ppf in zip(symbols,ppfiles):
        valences[s] = get_valence(psdir+'/'+ppf)
    
    #...write it out in QE format
    write_espresso_in(atoms,outfname,
                      symbols=symbols,
                      ppfiles=ppfiles,
                      valences=valences,
                      pitch=pitch,
                      constraints=constraints,
                      calculation=calc,
                      pseudo_dir=psdir,
                      ntyp=len(symbols))
    

