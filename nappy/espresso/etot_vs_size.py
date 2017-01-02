#!/opt/local/bin/python
"""
Calculate the total energy as a function of lattice constant,
by altering the lattice constant in in.pw file.
And if possible, calculate equilibrium lattice size and
bulk modulus, too.
Assume using Quantum ESPRESSO (QE).

PREPARE mode creates directories and in.pw files to be computed by QE.
ANALYZE mode reads the data obtained via QE and computes some values.

Usage:
  etot_vs_size.py prepare [options] POSCAR
  etot_vs_size.py analyze

Options:
  -h, --help  Shows this message and exit.
  -n NITER    Number of points to be calculated. [default: 11]
              In case of even number, original size is not used.
  --no-LS     Do not perform least square fitting. [default: False]
  -s STRAIN   Maximum strain value in a direction. [default: 0.05]
  -x          Change in x-direction. [default: False]
  -y          Change in y-direction. [default: False]
  -z          Change in z-direction. [default: False]
  --calc=CALC
              Calculation type. [default: scf]
  -o OUTFNAME
              Output file name. [default: in.pw]
  -p,--pitch=PITCH
              Pitch of k-points. [default: 0.0968]
"""

import sys,os,commands,copy
import numpy as np
from docopt import docopt
from ase.io import read
from scipy.optimize import leastsq
import json

sys.path.append(os.path.dirname(__file__))
from poscar2qein import write_espresso_in

__author__  = 'Ryo KOBAYASHI'
__version__ = '160922'
__licence__ = 'MIT'

_confname = 'conf.etot_vs_size.json'

def read_POSCAR(fname='POSCAR'):
    f=open(fname,'r')
    #...1st line: comment
    cmmt= f.readline()
    #...read 1st line and get current lattice size
    al= float(f.readline().split()[0])
    hmat= np.zeros((3,3))
    hmat[0]= [ float(x) for x in f.readline().split() ]
    hmat[1]= [ float(x) for x in f.readline().split() ]
    hmat[2]= [ float(x) for x in f.readline().split() ]
    buffer= f.readline().split()
    if buffer[0].isdigit():
        natm= 0
        for b in buffer:
            natm += int(b)
    else:
        natm= 0
        for b in f.readline().split():
            natm += int(b)
    f.close()
    return (al,hmat,natm)

def get_vol(al,hmat):
    a1= hmat[0:3,0] *al
    a2= hmat[0:3,1] *al
    a3= hmat[0:3,2] *al
    return np.dot(a1,np.cross(a2,a3))

def replace_lattice_constant(x,infname='POSCAR',outfname='POSCAR'):
    f=open(infname,'r')
    ini= f.readlines()
    f.close()
    g=open(outfname,'w')
    for l in range(len(ini)):
        if l == 1:
            g.write(' {0:10.4f}\n'.format(x))
        else:
            g.write(ini[l])
    g.close()

def replace_hmat(hmat,infname='POSCAR',outfname='POSCAR'):
    f=open(infname,'r')
    ini= f.readlines()
    f.close()
    g=open(outfname,'w')
    for l in range(len(ini)):
        if 2 <= l <= 4:
            g.write(' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(hmat[l-2,0],hmat[l-2,1],hmat[l-2,2]))
        else:
            g.write(ini[l])
    g.close()

def residuals(p,y,x):
    b,bp,v0,ev0= p
    err= y -( b*x/(bp*(bp-1.0)) *(bp*(1.0-v0/x) +(v0/x)**bp -1.0) +ev0 )
    return err

def peval(x,p):
    b,bp,v0,ev0= p
    return b*x/(bp*(bp-1.0)) *(bp*(1.0-v0/x) +(v0/x)**bp -1.0) +ev0

def prepare(args):
    niter= int(args['-n'])
    no_LS= args['--no-LS']
    mvx= args['-x']
    mvy= args['-y']
    mvz= args['-z']
    strain= float(args['-s'])
    sfac = float(args['-s'])
    calc = args['--calc']
    pitch = float(args['--pitch'])
    outfname = args['-o']
    poscar = args['POSCAR']

    if niter < 2:
        raise RuntimeError('NITER must be larger than 1')

    atoms0 = read(poscar,format='vasp')
    atoms = atoms0.copy()

    al_orig,hmat_orig,natm= read_POSCAR()
    hmat= copy.copy(hmat_orig)
    hmat_min= copy.copy(hmat_orig)
    hmat_max= copy.copy(hmat_orig)
    dhmat= np.zeros((3,3),dtype=float)
    if not mvx and not mvy and not mvz:
        al_min= al_orig*(1.0-strain)
        al_max= al_orig*(1.0+strain)
        dl= (al_max-al_min)/(niter-1)
    else:
        if mvx:
            hmat_min[0]= hmat_orig[0]*(1.0-strain)
            hmat_max[0]= hmat_orig[0]*(1.0+strain)
        if mvy:
            hmat_min[1]= hmat_orig[1]*(1.0-strain)
            hmat_max[1]= hmat_orig[1]*(1.0+strain)
        if mvz:
            hmat_min[2]= hmat_orig[2]*(1.0-strain)
            hmat_max[2]= hmat_orig[2]*(1.0+strain)
        dhmat= (hmat_max -hmat_min)/(niter-1)

    for iter in range(niter):
        dname= "etot_vs_size_{0:05d}".format(iter)
        print dname
        os.system("mkdir -p "+dname)
        os.system('cp INCAR KPOINTS POTCAR '+dname+'/')
        if not mvx and not mvy and not mvz:
            al= al_min +dl*iter
            hmat= hmat_orig
            replace_lattice_constant(al,infname='POSCAR',outfname=dname+'/POSCAR')
        else:
            al= al_orig
            hmat= hmat_min +dhmat*iter
            replace_hmat(hmat,infname='POSCAR',outfname=dname+'/POSCAR')
    print 'prepare done.'
    print ''
    print 'Perform VASP calculations in those directories and then run the following command,'
    print 'python etot_vs_size.py analyze'
    print ''

def analyze(config):
    """
    Analyze etot_vs_size using VASP results in etot_vs_size_##### directories.
    Configuration are read from config.etot_vs_size.json.
    """

    niter= int(config['-n'])
    no_LS= config['--no-LS']
    mvx= config['-x']
    mvy= config['-y']
    mvz= config['-z']
    strain= float(config['-s'])

    al_orig,hmat_orig,natm= read_POSCAR()
    hmat= copy.copy(hmat_orig)
    hmat_min= copy.copy(hmat_orig)
    hmat_max= copy.copy(hmat_orig)
    dhmat= np.zeros((3,3),dtype=float)
    if not mvx and not mvy and not mvz:
        al_min= al_orig*(1.0-strain)
        al_max= al_orig*(1.0+strain)
        dl= (al_max-al_min)/(niter-1)
    else:
        if mvx:
            hmat_min[0]= hmat_orig[0]*(1.0-strain)
            hmat_max[0]= hmat_orig[0]*(1.0+strain)
        if mvy:
            hmat_min[1]= hmat_orig[1]*(1.0-strain)
            hmat_max[1]= hmat_orig[1]*(1.0+strain)
        if mvz:
            hmat_min[2]= hmat_orig[2]*(1.0-strain)
            hmat_max[2]= hmat_orig[2]*(1.0+strain)
        dhmat= (hmat_max -hmat_min)/(niter-1)

    logfile= open('log.etot_vs_size','w')
    outfile1= open('out.etot_vs_size','w')
    for iter in range(niter):
        dname= "etot_vs_size_{0:05d}".format(iter)
        if not mvx and not mvy and not mvz:
            al= al_min +dl*iter
            hmat= hmat_orig
            #replace_lattice_constant(al)
        else:
            al= al_orig
            hmat= hmat_min +dhmat*iter
            #replace_hmat(hmat)
        erg= float(commands.getoutput("tail -n1 {0:s}/OSZICAR".format(dname) \
                                      +" | awk '{print $5}'"))
        #os.system("mkdir -p "+dname)
        #os.system("cp vasprun.xml {0}/".format(dname))
        vol= get_vol(al,hmat)
        print ' {0:10.4f} {1:10.4f} {2:15.7f}'.format(al,vol,erg)
        outfile1.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
        logfile.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
    outfile1.close()

    #...prepare for Murnaghan fitting
    f= open('out.etot_vs_size','r')
    lines= f.readlines()
    xarr= np.zeros((len(lines)))
    yarr= np.zeros((len(lines)))
    for l in range(len(lines)):
        dat= lines[l].split()
        xarr[l]= float(dat[1])
        yarr[l]= float(dat[2])
    f.close()
    #...set initial values
    b= 1.0
    bp= 2.0
    ev0= min(yarr)
    v0= xarr[len(xarr)/2]
    p0= np.array([b,bp,v0,ev0])
    #...least square fitting
    plsq= leastsq(residuals,p0,args=(yarr,xarr))

    #...output results
    print ' plsq=',plsq[0]
    print '{0:=^72}'.format(' RESULTS ')
    logfile.write('{0:=^72}\n'.format(' RESULTS '))
    a1= hmat_orig[0:3,0]
    a2= hmat_orig[0:3,1]
    a3= hmat_orig[0:3,2]
    uvol= np.dot(a1,np.cross(a2,a3))
    lc= (plsq[0][2]/uvol)**(1.0/3)
    print ' Lattice constant = {0:10.4f} Ang.'.format(lc)
    print ' Cohesive energy  = {0:10.3f} eV'.format(plsq[0][3]/natm)
    print ' Bulk modulus     = {0:10.2f} GPa'.format(plsq[0][0]*1.602e+2)
    logfile.write(' Lattice constant = {0:10.4f} Ang.\n'.format(lc))
    logfile.write(' Cohesive energy  = {0:10.3f} eV\n'.format(plsq[0][3]/natm))
    logfile.write(' Bulk modulus     = {0:10.2f} GPa\n'.format(plsq[0][0]*1.602e+2))

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * out.etot_vs_size'
    print ' * log.etot_vs_size'


if __name__ == '__main__':

    args= docopt(__doc__,version=__version__)

    if args['prepare']:
        with open(_confname,'w') as f:
            f.write(json.dumps(args,sort_keys=True,indent=4))
        prepare(args)
    elif args['analyze']:
        try:
            with open(_confname,'r') as f:
                config= json.load(f)
        except:
            raise RuntimeError('Cannot read '+_confname)
        analyze(config)

    
