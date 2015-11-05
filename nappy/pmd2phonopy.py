#!/bin/env python
"""
Compute forces of the supercell system created by phonopy and
create FORCE_SETS file, which includes forces on atoms.

Usage:
    pmd2phonopy.py [options] PMDFILE

Options:
    -h, --help  Show this help message and exit.
    -c, --cutoff=CUTOFF
                Cutoff radius of the interatomic potential used for
                the calculation of forces. [default: 4.0]
    --pmddir=PMDDIR
                Path to where ``pmd`` executable exists.
                [default: ~/src/nap/pmd/]
    --exec=EXEC
                Name of the executable. [default: pmd]
    -p, --plot  Show band graph with phonopy.
"""

import os,sys,math,copy,re
from docopt import docopt
import numpy as np
import glob,yaml

from atom import Atom
from pmdsys import PMDSystem

__author__="Ryo KOBAYASHI"
__version__="0.1a"

_default_band_conf="""ATOM_NAME = Si
DIM = 3 3 3
BAND = 0 0 0  1/2 0 1/2,  1/2 1/2 1  0 0 0  1/2 1/2 1/2
"""

################################################# Functions ############

def calc_extention_ratio(hmat,rcut):
    a1= np.zeros(3)
    a2= np.zeros(3)
    a3= np.zeros(3)
    # a1[:]= hmat[:,0]
    # a2[:]= hmat[:,1]
    # a3[:]= hmat[:,2]
    a1[:]= hmat[0,:]
    a2[:]= hmat[1,:]
    a3[:]= hmat[2,:]
    vol= abs(np.dot(a1,np.cross(a2,a3)))
    print(' vol of unit cell=',vol)
    a23= np.cross(a2,a3)
    a31= np.cross(a3,a1)
    a12= np.cross(a1,a2)
    area23= abs(np.sqrt(np.dot(a23,a23)))
    area31= abs(np.sqrt(np.dot(a31,a31)))
    area12= abs(np.sqrt(np.dot(a12,a12)))
    n1= int(2.0*rcut/vol*area23)+1
    n2= int(2.0*rcut/vol*area31)+1
    n3= int(2.0*rcut/vol*area12)+1
    return n1,n2,n3


################################################ Main routine ##########
if __name__ == "__main__":

    args= docopt(__doc__)
    rcut= float(args['--cutoff'])
    pmddir= args['--pmddir']
    execname= args['--exec']
    infname= args['PMDFILE']
    plot= args['--plot']
    #print(args)
    
    sys0= PMDSystem()
    sys0.read_pmd(infname)
    sys0.write_POSCAR()
    print ' POSCAR was written.'
    
    natm0= sys0.num_atoms()
    
    h0=np.zeros((3,3))
    h0[0]= sys0.a1 *sys0.alc
    h0[1]= sys0.a2 *sys0.alc
    h0[2]= sys0.a3 *sys0.alc
    n1,n2,n3= calc_extention_ratio(h0,rcut)
    print("n1,n2,n3 = {0:d} {1:d} {2:d}".format(n1,n2,n3))

    # Make displaced POSCARS via phonopy
    os.system("phonopy -d --dim=\"{0:d} {1:d} {2:d}\"".format(n1,n2,n3))
    with open('disp.yaml','r') as f:
        disp= yaml.load(f)
    
    # Make directories for the calculation of those POSCARS
    poscars= glob.glob("POSCAR-[0-9]??")
    n= 0
    for poscar in poscars:
        n+=1
        dname="disp-{0:03d}".format(n)
        os.system("mkdir -p {0}".format(dname))
        os.system("cp {0} {1}/POSCAR".format(poscar,dname))
        print("processing {0}...".format(dname))
        # Prepare pmd calculation
        psystmp= PMDSystem()
        psystmp.read_POSCAR(poscar)
        os.system("cp in.* {0}/".format(dname))
        if execname == "pmd":
            os.system("mkdir -p {0}/0000".format(dname))
            psystmp.write_pmd("{0}/0000/pmd00000".format(dname))
        elif execname == "smd":
            psystmp.write_pmd("{0}/smd0000".format(dname))
        os.system("cd {0}/; {1}/{2} > out.{2}; cd ..".format(dname,pmddir,execname))

    # Make FORCE_SETS file from the pmd results
    natm= int(disp['natom'])
    ndisps= len(disp['displacements'])
    with open('FORCE_SETS','w') as f:
        f.write("{0:d}\n".format(natm))
        f.write("{0:d}\n".format(ndisps))
        n=0
        for ds in disp['displacements']:
            n+=1
            f.write("\n")
            f.write("{0:d}\n".format(n))
            d= ds['displacement']
            f.write("{0:22.15f} {1:22.15f} {2:22.15f}\n".format(d[0],d[1],d[2]))
            # forces
            dname="disp-{0:03d}".format(n)
            g= open(dname+'/frc0000','r')
            g.readline() # skip 1st line, it should be same as natm
            for ia in range(natm):
                buff= [ float(a) for a in g.readline().split()]
                f.write(" {0:13.7f} {1:13.7f} {2:13.7f}\n".format(buff[0],buff[1],buff[2]))
            g.close()
    
    # Perform phonopy using FORCE_SETS
    if not os.path.exists("./band.conf"):
        print("There is no band.conf, which is configuration file for phonopy band calculation.")
        print("")
        print("Default band.conf is like:")
        print(__default_band_conf)
        sys.exit()

    # Change DIM parameters in band.conf
    dim= re.compile('DIM')
    with open('band.conf','r') as f:
        g= open('band.conf.tmp','w')
        for line in f.readlines():
            if dim.search(line):
                g.write('DIM = {0:2d} {1:2d} {2:2d}\n'.format(n1,n2,n3))
            else:
                g.write(line)
        g.close()
    os.system('mv band.conf.tmp band.conf')

    if plot:
        os.system("phonopy -p band.conf")
    else:
        os.system("phonopy band.conf")
        
    print("bandplot --gnuplot band.yaml > out_band")
    os.system("bandplot --gnuplot band.yaml > out_band")
