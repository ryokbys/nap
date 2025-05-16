#!/bin/env python
"""
Compute forces of the supercell system created by phonopy and
create FORCE_SETS file, which includes forces on atoms.

Usage:
    {0:s} [options] PMDINI

Options:
    -h, --help  Show this help message and exit.
    -c, --cutoff=CUTOFF
                Cutoff radius of the interatomic potential used for
                the calculation of forces. [default: 4.0]
    --exec=EXEC
                Name of the executable. [default: ~/src/nap/pmd/pmd]
    -p, --plot  Show band graph with phonopy.
"""
import os,sys,math,copy,re
from docopt import docopt
import numpy as np
import glob,yaml

import nappy

__author__="Ryo KOBAYASHI"
__version__="250515"

_default_band_conf="""ATOM_NAME = Si
DIM = 3 3 3
BAND = 0 0 0  1/2 0 1/2,  1/2 1/2 1  0 0 0  1/2 1/2 1/2
"""

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

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    rcut= float(args['--cutoff'])
    execname= args['--exec']
    infname= args['PMDINI']
    plot= args['--plot']
    #print(args)
    
    sys0= nappy.io.read(infname)
    nappy.io.write_POSCAR(sys0, fname='POSCAR')
    print(' POSCAR was written.')
    
    natm0= sys0.num_atoms()
    
    h0=np.zeros((3,3))
    h0[0]= sys0.a1 *sys0.alc
    h0[1]= sys0.a2 *sys0.alc
    h0[2]= sys0.a3 *sys0.alc
    n1,n2,n3= calc_extention_ratio(h0,rcut)
    print("n1,n2,n3 = {0:d} {1:d} {2:d}".format(n1,n2,n3))

    # Make displaced POSCARS via phonopy
    os.system(f"phonopy -d --dim=\"{n1:d} {n2:d} {n3:d}\"")
    with open('phonopy_disp.yaml','r') as f:
        disp= yaml.safe_load(f)
    
    # Make directories for the calculation of those POSCARS
    poscars= glob.glob("POSCAR-[0-9]??")
    n= 0
    for poscar in poscars:
        n+=1
        dname=f"disp-{n:03d}"
        os.system(f"mkdir -p {dname}")
        os.system(f"cp {poscar} {dname}/POSCAR")
        print(f"processing {dname}...")
        # Prepare pmd calculation
        psystmp = nappy.io.read_POSCAR(poscar)
        os.system(f"cp in.* {dname}/")
        nappy.io.write_pmd(psystmp, fname=f"{dname}/pmdini")
        basename=execname.split('/')[-1]
        os.system(f"cd {dname}/; {execname} > out.{basename}; cd ..")

    # Make FORCE_SETS file from the pmd results
    supercell = disp['supercell']
    #natm= int(disp['natom'])
    natm = len(supercell['points'])
    ndisps= len(disp['displacements'])
    with open('FORCE_SETS','w') as f:
        f.write(f"{natm:d}\n")
        f.write(f"{ndisps:d}\n")
        n=0
        for ds in disp['displacements']:
            n+=1
            f.write("\n")
            f.write(f"{n:d}\n")
            d= ds['displacement']
            f.write(f"{d[0]:22.15f} {d[1]:22.15f} {d[2]:22.15f}\n")
            # forces
            dname=f"disp-{n:03d}"
            nsystmp = nappy.io.read(f"{dname}/pmdfin", format='pmd')
            frcs = nsystmp.get_real_forces()
            for fi in frcs:
                f.write(f" {fi[0]:13.7f}  {fi[1]:13.7f}  {fi[2]:13.7f}\n")
            # g= open(dname+'/frc.pmd','r')
            # g.readline() # skip 1st line, it should be same as natm
            # for ia in range(natm):
            #     buff= [ float(a) for a in g.readline().split()]
            #     f.write(" {0:13.7f} {1:13.7f} {2:13.7f}\n".format(buff[0],buff[1],buff[2]))
            # g.close()
    
    # Perform phonopy using FORCE_SETS
    if not os.path.exists("./band.conf"):
        print("There is no band.conf, which is configuration file for phonopy band calculation.")
        print("")
        print("Default band.conf is like:")
        print(_default_band_conf)
        sys.exit()

    # Change DIM parameters in band.conf
    dim= re.compile('DIM')
    with open('band.conf','r') as f:
        g= open('band.conf.tmp','w')
        for line in f.readlines():
            if dim.search(line):
                g.write(f'DIM = {n1:2d} {n2:2d} {n3:2d}\n')
            else:
                g.write(line)
        g.close()
    os.system('mv band.conf.tmp band.conf')

    if plot:
        os.system("phonopy -p band.conf")
    else:
        os.system("phonopy band.conf")
        
    print("phonopy-bandplot --gnuplot band.yaml > out_band")
    os.system("phonopy-bandplot --gnuplot band.yaml > out_band")
    return None


if __name__ == "__main__":

    main()
