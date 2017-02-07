#!/bin/env python
"""
Make appropriate INCAR and KPOINTS files for VASP calculation
from POSCAR and POTCAR files. Number of division in each k-space direction
must be set, otherwise default NORMAL accuracy setting will be chosen.

Usage:
  prepare.py [options] POSCAR

Options:
  -h, --help  Show this help message and exit.
  -e, --even  Set even number to the k-points in a direction.
  -p PITCH    PITCH of k in a direction. [default: 0.2]
  --spin-polarize
              Set spin polarization true.
  --break-symmetry
              Set to allow symmetry breakage.
  --metal     Set metal flag True, and set ISMEAR=2,
              otherwise ISMEAR=-5.
  --potcar-dir POTCARDIR
              Specify the directory that contains POTCAR files.
              [default: ~/local/vasp/potpaw_PBE]
  --potcar-postfix POTCAR_POSTFIX
              Postfix of POTCAR directory name. If the species does not have this
              postfix, use directory without it. [default: ]
"""
from __future__ import print_function

import os
import math
from docopt import docopt

import nappy.vasp.poscar
import nappy.vasp.potcar

__author__ = "Ryo KOBAYASHI"
__version__ = "170206"


_SYSTEM='system made by prepare.py '+ __version__
_metal= False
_spin_polarized= False
_break_symmetry= False
_INCAR_name= 'INCAR'
_KPOINTS_name= 'KPOINTS'
_KPOINTS_type= 'Monkhorst-Pack' # or 'Gamma'

_IBRION= -1  # -1:no update, 0:MD, 1:q-Newton, 2:CG, 3:damped MD
_ISIF= 2     # 2: relax ions only, 3:shell-shape too, 4:shell volume too
_NSW= 0      # number of ion relaxation steps

_NCORE= 4

def determine_num_kpoint(b_length,pitch,leven):
    # maximum 11
    # minimum 1
    nk= int(2.0 *math.pi /b_length /pitch)
    if nk > 11:
        return 11
    elif nk < 1:
        return 1
    if leven:
        if nk % 2 == 1:
            nk= nk +1
    else:
        if nk % 2 == 0:
            nk= nk +1
    return nk

def write_KPOINTS(fname,type,ndiv):
    with open(fname,'w') as f:
        f.write('{0:d}x{1:d}x{1:d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
        f.write('0\n')
        f.write(type+'\n')
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(0,0,0))
        f.close()

def write_INCAR(fname,encut,nbands,break_symmetry,spin_polarized,metal):
    
    with open(fname,'w') as f:
        f.write("SYSTEM ="+_SYSTEM+"\n")
        f.write("\n")
        f.write("ISTART = 1\n")
        f.write("ICHARG = 1\n")
        f.write("INIWAV = 1\n")
        if spin_polarized:
            f.write("ISPIN  = 2\n")
            f.write("IMIX     = 4\n")
            f.write("AMIX     = 0.05\n")
            f.write("BMIX     = 0.0001\n")
            f.write("AMIX_MAG = 0.2\n")
            f.write("BMIX_MAG = 0.0001\n")
            f.write("MAXMIX   = 40\n")
        else:
            f.write("ISPIN  = 1\n")
    
        if break_symmetry:
            f.write("ISYM   = 0\n")
        else:
            f.write("ISYM   = 2\n")
    
        f.write("\n")
        f.write("ENCUT  = {0:7.3f}\n".format(encut))
        f.write("LREAL  = Auto\n")
        f.write("EDIFF  = 1.0e-6\n")
        f.write("ALGO   = Fast\n")
        f.write("PREC   = Normal\n")
        f.write("\n")
        f.write("NELMIN = 4\n")
        f.write("NELM   = 100\n")
        f.write("NBANDS = {0:4d}\n".format(nbands))
        f.write("\n")
        if metal:
            f.write("ISMEAR = 2\n")
            f.write("SIGMA  = 0.2\n")
        else:
            f.write("ISMEAR = -5\n")
            f.write("SIGMA  = 0.00001\n")
    
        f.write("\n")
        f.write("ISIF   = {0:2d}\n".format(_ISIF))
        f.write("IBRION = {0:2d}\n".format(_IBRION))
        f.write("POTIM  = 0.5\n") 
        f.write("SMASS  = 0.4\n") 
        f.write("NSW    = {0:4d}\n".format(_NSW))
        f.write("\n")
        
        f.write("NCORE   = {0:4d}\n".format(_NCORE)) 
        f.write("\n")
        f.close()

def prepare_potcar(poscar,potcar_dir,potcar_postfix):
    """
    Create a POTCAR if there is not in the directory.
    The directory that contains pseudo-potential files should be specified.
    Potentials without `_h`, `_s`, or `_GW` are to be used.
    """
    if not potcar_dir:
        return None
    if not os.path.exists(potcar_dir):
        raise RuntimeError(potcar_dir+' does not exist.')
    if len(poscar.species) == 0:
        return None
    if os.path.exists('POTCAR'):
        os.system('rm POTCAR')
    for sp in poscar.species:
        spdir = potcar_dir+'/'+sp+potcar_postfix
        if not os.path.exists(spdir):
            raise RuntimeError(spdir+' does not exist.')
        os.system('cat '+spdir+'/POTCAR >> ./POTCAR')
    return None

def prepare_vasp(poscar_fname,pitch,even,spin_polarized,break_symmetry,
                 metal,potcar_dir,potcar_postfix):
    
    print(' Pitch of k points = {0:5.1f}'.format(pitch))

    poscar= nappy.vasp.poscar.POSCAR(poscar_fname)
    poscar.read(poscar_fname)

    if os.path.exists('./POTCAR'):
        potcar = nappy.vasp.potcar.read_POTCAR()
    else:
        prepare_potcar(poscar,potcar_dir,potcar_postfix)
        potcar = nappy.vasp.potcar.read_POTCAR()
    species= potcar['species']
    encut= max(potcar['encut'])
    valences= potcar['valence']
    a1= poscar.h[:,0]
    a2= poscar.h[:,1]
    a3= poscar.h[:,2]
    al= poscar.afac
    natms= poscar.num_atoms

    print(" species:",species)
    print(" encut:",encut)
    print(" valences:",valences)
    print(" natms:",natms)
    ntot= 0
    nele= 0
    for i in range(len(natms)):
        ntot= ntot +natms[i]
        nele= nele +natms[i]*int(valences[i])
    
    if spin_polarized:
        nbands= int(nele/2 *1.8)
    else:
        nbands= int(nele/2 *1.4)
    
    if nbands < 50:
        nbands= nele

    l1= al *math.sqrt(a1[0]**2 +a1[1]**2 +a1[2]**2)
    l2= al *math.sqrt(a2[0]**2 +a2[1]**2 +a2[2]**2)
    l3= al *math.sqrt(a3[0]**2 +a3[1]**2 +a3[2]**2)
    print(' Length of each axes:')
    print('   l1 = {0:10.3f}'.format(l1))
    print('   l2 = {0:10.3f}'.format(l2))
    print('   l3 = {0:10.3f}'.format(l3))
    k1= determine_num_kpoint(l1,pitch,even)
    k2= determine_num_kpoint(l2,pitch,even)
    k3= determine_num_kpoint(l3,pitch,even)
    print(' Number of k-points: {0:2d} {1:2d} {2:2d}'.format(k1,k2,k3))
    ndiv= [k1,k2,k3]
    
    write_KPOINTS(_KPOINTS_name,_KPOINTS_type,ndiv)
    write_INCAR(_INCAR_name,encut,nbands,break_symmetry,
                spin_polarized,metal)

#=======================================================================

if __name__ == '__main__':

    args= docopt(__doc__)

    pitch= float(args['-p'])
    leven= args['--even']
    _spin_polarized= args['--spin-polarize']
    _break_symmetry= args['--break-symmetry']
    _metal= args['--metal']
    poscar_fname= args['POSCAR']
    potcar_dir = os.path.expanduser(args['--potcar-dir'])
    potcar_postfix = args['--potcar-postfix']

    prepare_vasp(poscar_fname,pitch,leven,_spin_polarized,_break_symmetry,
                 _metal,potcar_dir,potcar_postfix)
