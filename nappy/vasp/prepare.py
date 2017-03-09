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
  --encut ENCUT
              Cutoff energy. [default: None]
  --ediff EDIFF
              Convergence criteria for the energy difference. [default: 1.0e-6]
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
  --relax     Enable structure relaxation. This sets IBRION=1 and NSW=100.
  --relax-cell RELAX_CELL
              Type of relaxation of cell, either (ion_only, cell, volume)
              [default: ion_only]
  --isif ISIF  Directly specify ISIF value. [default: 2]
"""
from __future__ import print_function

import os
import math
from docopt import docopt

import nappy.vasp.poscar
import nappy.vasp.potcar

__author__ = "Ryo KOBAYASHI"
__version__ = "170206"

_magnetic_elements = ('Cr','Mn','Fe','Co','Ni')

_INCAR_name= 'INCAR'
_KPOINTS_name= 'KPOINTS'
_KPOINTS_type= 'Monkhorst-Pack' # or 'Gamma'

_IBRION= -1  # -1:no update, 0:MD, 1:q-Newton, 2:CG, 3:damped MD
_ISIF= 2     # 2: relax ions only, 3:shell-shape too, 4:shell volume too
_NSW= 0      # number of ion relaxation steps

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

def estimate_ncore(nbands):
    """
    Estimate NCORE value from NBANDS info.
    NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE). 
    VASP master recommends that 

        NCORE = 4 - approx. SQRT( # of cores )

    And according to the site below,
    https://www.nsc.liu.se/~pla/blog/2015/01/12/vasp-how-many-cores/
    NCORE ~ NBANDS/8 is also a good approximation.
    And we limit NCORE as a multiplier of 2.
    """
    nb8 = int(nbands/8)
    if nb8 < 4:
        ncore = 4
        return ncore
    else:
        ncore = 4
        for i in range(20):
            if ncore*2 > nb8:
                return ncore
            ncore *= 2
    raise RuntimeError('Something is wrong.')

def estimate_nbands(nel):
    """
    Estimate appropriate number of bands to be computed.

    The policy is that the NBANDS:
    - should be multiples of 4 considering the efficient
      parallelization (which cannot be taken account here, though),
    - should be not be less than number of electrons, NEL, 
      if it is less than 50,
    """
    nbands = nel
    if nbands % 4 != 0:
        nbands += nbands % 4
    return nbands

def write_KPOINTS(fname,type,ndiv):
    with open(fname,'w') as f:
        f.write('{0:d}x{1:d}x{1:d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
        f.write('0\n')
        f.write(type+'\n')
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(0,0,0))
        f.close()

def write_INCAR(fname,encut,nbands,break_symmetry,spin_polarized,metal,
                ediff,
                relax=None,relax_cell=None,isif=2):
    SYSTEM = 'system made by prepare.py '+ __version__
    
    with open(fname,'w') as f:
        f.write("SYSTEM = "+SYSTEM+"\n")
        f.write("\n")
        f.write("ISTART = 1\n")
        f.write("ICHARG = 1\n")
        f.write("INIWAV = 1\n")
        if spin_polarized:
            f.write("ISPIN    = 2\n")
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
        f.write("EDIFF  = {0:7.1e}\n".format(ediff))
        f.write("ALGO   = Normal\n")
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
        if relax_cell == 'cell':
            f.write("ISIF   = {0:2d}\n".format(3))
        elif relax_cell == 'volume':
            f.write("ISIF   = {0:2d}\n".format(4))
        else:
            f.write("ISIF   = {0:2d}\n".format(isif))
        if relax:
            f.write("IBRION = {0:2d}\n".format(1))
            f.write("NSW    = {0:4d}\n".format(100))
        else:
            f.write("IBRION = {0:2d}\n".format(_IBRION))
            f.write("NSW    = {0:4d}\n".format(_NSW))
        f.write("POTIM  = 0.5\n") 
        f.write("SMASS  = 0.4\n") 
        f.write("\n")

        ncore = estimate_ncore(nbands)
        f.write("NCORE  = {0:4d}\n".format(ncore)) 
        f.write("\n")
        f.close()

def prepare_potcar(poscar,potcar_dir,potcar_postfix=''):
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

    postfixes = [potcar_postfix, '', '_sv', '_pv', '_s', '_h']
    for sp in poscar.species:
        found = False
        for pf in postfixes:
            spdir = potcar_dir+'/'+sp +pf
            if os.path.exists(spdir):
                found = True
            else:
                continue
            os.system('cat '+spdir+'/POTCAR >> ./POTCAR')
            print('POTCAR for {0:3s}: {1:s}'.format(sp,spdir))
            break
        if not found:
            raise RuntimeError('POTCAR was not found for '+sp)
    return None

def prepare_vasp(poscar_fname,pitch,even,spin_polarized,break_symmetry,
                 metal,potcar_dir,potcar_postfix,
                 encut=None,ediff=None,
                 relax=None,relax_cell=None,isif=None):
    
    print(' Pitch of k points = {0:5.1f}'.format(pitch))

    poscar= nappy.vasp.poscar.POSCAR(poscar_fname)
    poscar.read(poscar_fname)

    if os.path.exists('./POTCAR'):
        potcar = nappy.vasp.potcar.read_POTCAR()
    else:
        prepare_potcar(poscar,potcar_dir,potcar_postfix)
        potcar = nappy.vasp.potcar.read_POTCAR()
    species= potcar['species']
    valences= potcar['valence']
    a1= poscar.h[:,0]
    a2= poscar.h[:,1]
    a3= poscar.h[:,2]
    al= poscar.afac
    natms= poscar.num_atoms
    if not encut:
        encut= max(potcar['encut'])

    print(" species:",species)
    print(" encut:",encut)
    print(" valences:",valences)
    print(" natms:",natms)
    ntot= 0
    nele= 0
    for i in range(len(natms)):
        ntot= ntot +natms[i]
        nele= nele +natms[i]*int(valences[i])

    for e in _magnetic_elements:
        if e in species:
            spin_polarized = True
    
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

    nbands = estimate_nbands(nele)

    write_KPOINTS(_KPOINTS_name,_KPOINTS_type,ndiv)
    write_INCAR(_INCAR_name,encut,nbands,break_symmetry,
                spin_polarized,metal,
                ediff,relax=relax,relax_cell=relax_cell,isif=isif)

#=======================================================================

if __name__ == '__main__':

    args= docopt(__doc__)

    pitch= float(args['-p'])
    leven= args['--even']
    spin_polarized= args['--spin-polarize']
    break_symmetry= args['--break-symmetry']
    metal= args['--metal']
    poscar_fname= args['POSCAR']
    potcar_dir = os.path.expanduser(args['--potcar-dir'])
    potcar_postfix = args['--potcar-postfix']
    encut = args['--encut']
    ediff = float(args['--ediff'])
    relax = args['--relax']
    relax_cell = args['--relax-cell']
    isif = int(args['--isif'])

    if encut[0].isdigit():
        encut = float(encut)
    else:
        encut = None

    prepare_vasp(poscar_fname,pitch,leven,spin_polarized,break_symmetry,
                 metal,potcar_dir,potcar_postfix,
                 encut=encut,ediff=ediff,
                 relax=relax,relax_cell=relax_cell,isif=isif)
