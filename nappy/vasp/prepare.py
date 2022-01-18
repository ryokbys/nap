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
  --ediffg EDIFFG
              Convergence criteria for ionic relaxation. 
              Negative value for force criterion. [default: -0.05]
  --spin-polarize
              Set spin polarization true. If '--high-spin' is set, this is also set.
  --high-spin 
              Set initial spin state as high. Otherwise set it low.
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
  --mode MODE
              Mode of VASP calculation, either one of the followings:
              scf, relax-ion, relax-cell, relax-shape, md-ion, md-cell, md-shape.
              [default: scf]
  --nsw NSW   Number of MD/relaxation steps. [default: 0]
  --isif ISIF  Directly specify ISIF value. [default: 2]
  --ismear ISMEAR
              Specify smearing type. [default: 0]
  --sigma SIGMA
              Sigma (temperature) of smearing. [default: 0.01]
  --extra-nbands REXTNBANDS
              Ratio of extra number of bands multiplied to the estimated number of bands.
              [default: 1.0]
"""
from __future__ import print_function

import os
import math
from docopt import docopt

import nappy.vasp.poscar
import nappy.vasp.potcar

__author__ = "Ryo KOBAYASHI"
__version__ = "190425"

_magnetic_elements = ('Cr','Mn','Fe','Co','Ni')

_INCAR_name= 'INCAR'
_KPOINTS_name= 'KPOINTS'
_KPOINTS_type= 'Monkhorst-Pack'  # or 'Gamma'

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
    But we cannot know how many cores will be used when preparing inputs.
    Here we use the following original criterion, 
       NCORE ~ SQRT(NBANDS)/2
    """
    nbo = int(math.sqrt(nbands)/2)
    if nbo < 4:
        ncore = 4
        return ncore
    else:
        ncore = 4
        while True:
            if ncore >= nbo:
                return ncore
            ncore += 2
            if ncore >= nbands:
                raise RuntimeError('ncore >= nbands')
    raise RuntimeError('Could not estimate NCORE.')

def estimate_kpar(kpnts):
    """
    Estimte KPAR which is the number of parallel process among k-points.
    """
    return max(kpnts)

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

def get_magmom_str(high_spin,species=[],natms=[],valences=[]):
    magmom = ""
    if high_spin:
        for i,spcs in enumerate(species):
            n = natms[i]
            spin = valences[i]
            magmom += ' {0:d}*{1:d}'.format(n,int(spin))
    else:
        for i,spcs in enumerate(species):
            n = natms[i]
            magmom += ' {0:d}*{1:4.2f}'.format(n,0.1)
    return magmom

def write_KPOINTS(fname,ktype,kpnts):
    with open(fname,'w') as f:
        f.write('{0:d}x{1:d}x{1:d}\n'.format(*kpnts))
        f.write('0\n')
        f.write(ktype+'\n')
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(*kpnts))
        f.write(' {0:2d} {1:2d} {2:2d}\n'.format(0,0,0))
        f.close()
    return None

def write_INCAR(fname,encut,nbands,break_symmetry,
                spin_polarized,ediff,ediffg,
                high_spin,species,natms,valences,kpnts,
                mode=None,nsw=0,isif=2,ismear=0,sigma=0.01):
    from datetime import datetime as dt
    tdate = dt.now()
    dstr = tdate.strftime('%Y-%m-%d')
    SYSTEM = 'System made by nappy/vasp/prepare.py '+ dstr
    
    if mode == 'scf' and nsw != 0:
        print('NSW is meaningless for mode==scf, so reset NSW=0.')
        nsw = 0
    elif mode != 'scf' and nsw == 0:
        print(' Since NSW==0 is meaningless for mode==(relax|md), change NSW=1000.')
        nsw = 1000
        
    with open(fname,'w') as f:
        f.write("SYSTEM = "+SYSTEM+"\n")
        f.write("\n")
        f.write("ISTART = 1  # startjob:  0) new,  1) cont,  2) samecut \n")
        f.write("ICHARG = 1  # charge:  1) file,  2) atom,  10) const \n")
        f.write("INIWAV = 1  # initial wavefunc: 0) lowe,  1) random \n")
        if spin_polarized:
            f.write("\n")
            f.write("# Spin poralization and mixing information \n")
            f.write("ISPIN    = 2  # 2) spin-polarized,  1) not \n")
            f.write("IMIX     = 4 \n")
            f.write("AMIX     = 0.05 \n")
            f.write("BMIX     = 0.0001 \n")
            f.write("AMIX_MAG = 0.2 \n")
            f.write("BMIX_MAG = 0.0001 \n")
            f.write("MAXMIX   = 40 \n")
            magmom = get_magmom_str(high_spin,species,natms,valences)
            f.write("MAGMOM   = {0:s} \n".format(magmom))
        else:
            f.write("ISPIN  = 1  # 1) no-spin,  2) spin-polarized \n")
    
        if break_symmetry:
            f.write("ISYM   = 0  # 0) symmetry OFF,  1-3) symmetry ON, \n")
        else:
            f.write("ISYM   = 2  # 0) symmetry OFF,  1-3) symmetry ON, \n")
    
        f.write("\n")
        f.write("ENCUT  =  {0:7.3f}\n".format(encut))
        f.write("LREAL  =  Auto  # non-local projectors in real space \n")
        f.write("EDIFF  =  {0:9.1e}\n".format(ediff))
        if 'relax' in mode:
            f.write("EDIFFG =  {0:10.2e}\n".format(ediffg))
        f.write("ALGO   =  Very Fast\n")
        f.write("PREC   =  High\n")
        if 'md' in mode:
            f.write("LWAVE  =  F  # Not to write wave function data into file. \n")
            f.write("LCHARG =  F  # Not to write charge distribution into file. \n")
            f.write("NWRITE =  0  # Reduce verbosity \n")
        else:
            f.write("NWRITE =  2  # Default verbosity \n")
        f.write("\n")
        f.write("NELMIN =  4  # minimum electronic steps \n")
        f.write("NELM   =  100  # maximum electronic steps \n")
        f.write("NBANDS = {0:4d}  # number of bands \n".format(nbands))
        f.write("\n")
        # if metal:
        #     f.write("ISMEAR =   2\n")
        #     f.write("SIGMA  =   0.2\n")
        # else:
        #     f.write("ISMEAR =   0\n")
        #     f.write("SIGMA  =   0.01\n")
        f.write("ISMEAR =   {0:d}  # method for partial occupancy: -5) Blochl, -1) Fermi, 0) Gaussian, 1-5) Methfessel-Paxton order-N \n".format(ismear))
        f.write("SIGMA  =   {0:8.4f}  # width of Gaussian\n".format(sigma))
    
        f.write("\n")
        f.write("# Ionic updates IBRION: -1) no update, 0) MD, 1) quasi-Newton, 2) CG \n")
        if 'scf' in mode:
            f.write("IBRION = {0:4d}\n".format(-1))
            f.write("NSW    = {0:4d}\n".format(0))
        elif 'relax' in mode:
            f.write("IBRION = {0:4d}\n".format(1))
            f.write("NSW    = {0:4d}\n".format(nsw))
            f.write("POTIM  = 0.5\n")
        elif 'md' in mode:
            f.write("IBRION = {0:4d}\n".format(0))
            f.write("NSW    = {0:4d}\n".format(nsw))
            f.write("POTIM  = 1.0\n")
            f.write("SMASS  = 0.4\n")
            f.write("MDALGO = 3   # 0)NVE, 1)Andersen thermostat, 2)Nose-Hoover, 3)Langevin \n")
            gmms = [ 10.0 for i in range(len(species)) ]
            f.write("LANGEVIN_GAMMA = ")
            for g in gmms:
                f.write("{0:4.1f} ".format(g))
            f.write("\n")
            if 'cell' in mode or 'shape' in mode:
                f.write("LANGEVIN_GAMMA_L = 10.0\n")
            f.write("TEBEG  = 1000.0\n")
            f.write("# TEEND  = 1000.0\n")
        else:
            raise ValueError('mode is wrong, mode = ',mode)

        f.write("# What to relax, ISIF: 2) ion-only,  3) ion & cell,  4) ion & cell-shape (vol conserving) \n")
        if 'ion' in mode:
            f.write("ISIF   = {0:4d}\n".format(2))
        elif 'cell' in mode:
            f.write("ISIF   = {0:4d}\n".format(3))
        elif 'shape' in mode:
            f.write("ISIF   = {0:4d}\n".format(4))
        else:
            f.write("ISIF   = {0:4d}\n".format(2))

        f.write("\n")

        ncore = estimate_ncore(nbands)
        f.write("NCORE  = {0:4d}\n".format(ncore))
        kpar = estimate_kpar(kpnts)
        f.write("KPAR   = {0:4d}\n".format(kpar))
        f.write("\n")
        f.close()

    return None


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

    if potcar_postfix == '_sv':
        postfixes = [potcar_postfix, '_pv', '', '_s', '_h']
    elif potcar_postfix == '_pv':
        postfixes = [potcar_postfix, '', '_sv', '_s', '_h']
    else:
        postfixes = [potcar_postfix, '', '_pv', '_sv', '_s', '_h']
    for sp in poscar.species:
        found = False
        for pf in postfixes:
            spdir = potcar_dir+'/'+sp +pf
            if os.path.exists(spdir):
                found = True
            else:
                continue
            os.system('cat '+spdir+'/POTCAR >> ./POTCAR')
            print(' POTCAR for {0:3s}: {1:s}'.format(sp,spdir))
            break
        if not found:
            raise RuntimeError('POTCAR was not found for '+sp)
    return None

def prepare_vasp(poscar_fname,pitch,even,spin_polarized,break_symmetry,
                 potcar_dir,potcar_postfix,
                 encut=None,ediff=None,ediffg=-0.05,
                 mode=None,nsw=0,isif=None,high_spin=False,
                 ismear=0,sigma=0.01,extra_nbands=1.0):
    print(' Pitch of k points = {0:5.1f}'.format(pitch))

    poscar= nappy.vasp.poscar.POSCAR(poscar_fname)
    poscar.read(poscar_fname)

    if os.path.exists('./POTCAR'):
        print(' WARNING:')
        print('   Use ./POTCAR already exists.')
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
    kpnts= [k1,k2,k3]

    nbands = estimate_nbands(nele)
    nbands = int(nbands*extra_nbands)

    write_KPOINTS(_KPOINTS_name,_KPOINTS_type,kpnts)
    write_INCAR(_INCAR_name,encut,nbands,break_symmetry,
                spin_polarized,ediff,ediffg,
                high_spin,species,natms,valences,kpnts,
                mode=mode,nsw=nsw,isif=isif,ismear=ismear,sigma=sigma)
    return None

def main():
    args= docopt(__doc__)

    pitch= float(args['-p'])
    leven= args['--even']
    spin_polarized= args['--spin-polarize']
    high_spin = args['--high-spin']
    break_symmetry= args['--break-symmetry']
    metal= args['--metal']
    poscar_fname= args['POSCAR']
    potcar_dir = os.path.expanduser(args['--potcar-dir'])
    potcar_postfix = args['--potcar-postfix']
    encut = args['--encut']
    ediff = float(args['--ediff'])
    ediffg = float(args['--ediffg'])
    mode = args['--mode']
    nsw  = int(args['--nsw'])
    isif = int(args['--isif'])
    ismear = int(args['--ismear'])
    sigma = float(args['--sigma'])
    rextnb = float(args['--extra-nbands'])

    if encut[0].isdigit():
        encut = float(encut)
    else:
        encut = None

    if high_spin:
        spin_polarized = True

    if metal:  # if metal is specified, ismear and sigma are overwritten
        ismear = 2
        sigma = 0.2
        print('Since metal==True, ISMEAR and SIGMA are overwritten.')

    prepare_vasp(poscar_fname,pitch,leven,spin_polarized,break_symmetry,
                 potcar_dir,potcar_postfix,
                 encut=encut,ediff=ediff,ediffg=ediffg,
                 mode=mode,nsw=nsw,isif=isif,
                 high_spin=high_spin,
                 ismear=ismear,sigma=sigma,extra_nbands=rextnb)
    return None

if __name__ == '__main__':
    main()
    
