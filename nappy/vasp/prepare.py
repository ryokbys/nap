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
"""

import math,optparse
from docopt import docopt

import poscar, potcar

__author__ = "Ryo KOBAYASHI"
__version__ = "0.1a"


_SYSTEM='system made by prepare-vasp.py '+ __version__
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
    f=open(fname,'w')
    f.write('{0:d}x{1:d}x{1:d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
    f.write('0\n')
    f.write(type+'\n')
    f.write(' {0:2d} {1:2d} {2:2d}\n'.format(ndiv[0],ndiv[1],ndiv[2]))
    f.write(' {0:2d} {1:2d} {2:2d}\n'.format(0,0,0))
    f.close()

def write_INCAR(fname,encut,nbands):
    file=open(fname,'w')
    file.write("SYSTEM ="+_SYSTEM+"\n")
    file.write("\n")
    file.write("ISTART = 1\n")
    file.write("ICHARG = 1\n")
    file.write("INIWAV = 1\n")
    if _spin_polarized:
        file.write("ISPIN  = 2\n")
        file.write("IMIX     = 4\n")
        file.write("AMIX     = 0.05\n")
        file.write("BMIX     = 0.0001\n")
        file.write("AMIX_MAG = 0.2\n")
        file.write("BMIX_MAG = 0.0001\n")
        file.write("MAXMIX   = 40\n")
    else:
        file.write("ISPIN  = 1\n")

    if _break_symmetry:
        file.write("ISYM   = 0\n")
    else:
        file.write("ISYM   = 2\n")

    file.write("\n")
    file.write("ENCUT  = {0:7.3f}\n".format(encut))
    file.write("LREAL  = Auto\n")
    file.write("EDIFF  = 1.0e-6\n")
    file.write("ALGO   = Fast\n")
    file.write("PREC   = Normal\n")
    file.write("\n")
    file.write("NELMIN = 4\n")
    file.write("NELM   = 100\n")
    file.write("NBANDS = {0:4d}\n".format(nbands))
    file.write("\n")
    if _metal:
        file.write("ISMEAR = 2\n")
        file.write("SIGMA  = 0.2\n")
    else:
        file.write("ISMEAR = -5\n")
        file.write("SIGMA  = 0.00001\n")

    file.write("\n")
    file.write("ISIF   = {0:2d}\n".format(_ISIF))
    file.write("IBRION = {0:2d}\n".format(_IBRION))
    file.write("POTIM  = 0.5\n") 
    file.write("SMASS  = 0.4\n") 
    file.write("NSW    = {0:4d}\n".format(_NSW))
    file.write("\n")
    
    file.write("NCORE   = {0:4d}\n".format(_NCORE)) 
    file.write("\n")
    file.close()

def check_POTCAR():
    

#=======================================================================

if __name__ == '__main__':

    args= docopt(__doc__)

    pitch= float(args['-p'])
    leven= args['--even']
    _spin_polarized= args['--spin-polarize']
    _break_symmetry= args['--break-symmetry']
    _metal= args['--metal']
    poscar_fname= args['POSCAR']

    print ' Pitch of k points = {0:5.1f}'.format(pitch)

    poscar= poscar.POSCAR()
    poscar.read(poscar_fname)

    
    potcar= potcar.read_POTCAR()
    species= potcar['species']
    encut= max(potcar['encut'])
    valences= potcar['valence']
    a1= poscar.h[:,0]
    a2= poscar.h[:,1]
    a3= poscar.h[:,2]
    al= poscar.afac
    natms= poscar.num_atoms

    print " species:",species
    print " encut:",encut
    print " valences:",valences
    print " natms:",natms
    ntot= 0
    nele= 0
    for i in range(len(natms)):
        ntot= ntot +natms[i]
        nele= nele +natms[i]*int(valences[i])
    
    if _spin_polarized:
        nbands= int(nele/2 *1.8)
    else:
        nbands= int(nele/2 *1.4)
    
    if nbands < 50:
        nbands= nele

    l1= al *math.sqrt(a1[0]**2 +a1[1]**2 +a1[2]**2)
    l2= al *math.sqrt(a2[0]**2 +a2[1]**2 +a2[2]**2)
    l3= al *math.sqrt(a3[0]**2 +a3[1]**2 +a3[2]**2)
    print ' Length of each axes:'
    print '   l1 = {0:10.3f}'.format(l1)
    print '   l2 = {0:10.3f}'.format(l2)
    print '   l3 = {0:10.3f}'.format(l3)
    k1= determine_num_kpoint(l1,pitch,leven)
    k2= determine_num_kpoint(l2,pitch,leven)
    k3= determine_num_kpoint(l3,pitch,leven)
    print ' Number of k-points: {0:2d} {1:2d} {2:2d}'.format(k1,k2,k3)
    ndiv= [k1,k2,k3]
    
    write_KPOINTS(_KPOINTS_name,_KPOINTS_type,ndiv)
    write_INCAR(_INCAR_name,encut,nbands)
