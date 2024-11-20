#!/opt/local/bin/python
"""
Calculate the total energy as a function of lattice constant,
by altering the lattice constant in POSCAR file.
And if possible, calculate equilibrium lattice size and
bulk modulus, too.

Usage:
  energy_vs_size.py [options]

Options:
  -h, --help  Show this help message and exit.
  -n NITER    Number of points to be calculated. [default: 10]
  -p          Show a graph on the screen.
  -s STRAIN   STRAIN (%) applied to the lattice constant. [default: 10.0]
  -x          Not isotropic deformation, *x* direction is changed.
  -y          Not isotropic deformation, *y* direction is changed.
  -z          Not isotropic deformation, *z* direction is changed.
  --LS        Perform least square fitting to obtain lattice constant
              and bulk modulus.
  --cmd=CMD   VASP execution command. [default: \"vasp > out.vasp\"]
"""
import sys,os,copy
import subprocess
from docopt import docopt
import numpy as np

from scipy.optimize import leastsq
_no_pyplot=False
try:
    import matplotlib.pyplot as plt
except:
    _no_pyplot=True

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

def replace_1st_line(x,fname='POSCAR'):
    f=open(fname,'r')
    ini= f.readlines()
    f.close()
    g=open(fname,'w')
    for l in range(len(ini)):
        if l == 1:
            g.write(' {0:10.4f}\n'.format(x))
        else:
            g.write(ini[l])
    g.close()

def replace_hmat(hmat,fname='POSCAR'):
    f=open(fname,'r')
    ini= f.readlines()
    f.close()
    g=open(fname,'w')
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

if __name__ == '__main__':
    
    args= docopt(__doc__)

    niter= int(args['-n'])
    show_graph= args['-p']
    strain= float(args['-s'])
    perform_ls= args['--LS']
    cmd= args['--cmd']
    mvx= args['-x']
    mvy= args['-y']
    mvz= args['-z']

    if show_graph and _no_pyplot:
        print("matplotlib.pyplot is not available in this sysytem.")
        print("Run this script without -p option.")
        sys.exit()

    strain= strain/100
    al_orig,hmat_orig,natm= read_POSCAR()
    hmat= copy.copy(hmat_orig)
    hmat_min= copy.copy(hmat_orig)
    hmat_max= copy.copy(hmat_orig)
    dhmat= np.zeros((3,3),dtype=float)
    if not mvx and not mvy and not mvz:
        al_min= al_orig*(1.0-strain)
        al_max= al_orig*(1.0+strain)
        dl= (al_max-al_min)/niter
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
        dhmat= (hmat_max -hmat_min)/niter

    logfile= open('log.energy_vs_size','w')
    outfile1= open('out.energy_vs_size','w')
    for iter in range(niter+1):
        dname= "energy-{0:05d}".format(iter)
        if not mvx and not mvy and not mvz:
            al= al_min +dl*iter
            hmat= hmat_orig
            replace_1st_line(al)
        else:
            al= al_orig
            hmat= hmat_min +dhmat*iter
            replace_hmat(hmat)
        #os.system('vasp > out.vasp')
        os.system(cmd)
        erg= float(subprocess.getoutput("tail -n1 OSZICAR | awk '{print $5}'"))
        os.system("mkdir -p "+dname)
        os.system("cp INCAR POSCAR OSZICAR OUTCAR vasprun.xml {0}/".format(dname))
        vol= get_vol(al,hmat)
        print(' {0:10.4f} {1:10.4f} {2:15.7f}'.format(al,vol,erg))
        outfile1.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
        logfile.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
    outfile1.close()

    if not mvx and not mvy and not mvz:
        replace_1st_line(al_orig)
    else:
        replace_hmat(hmat_orig)
        print(' energy_vs_size finished because mvx,y,z=False.')
        sys.exit()

    if not perform_ls:
        print(' energy_vs_size finished without performing least square fitting...')
        sys.exit()

    #...prepare for Murnaghan fitting
    f= open('out.energy_vs_size','r')
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
    print(' plsq=',plsq[0])
    print('{0:=^72}'.format(' RESULTS '))
    logfile.write('{0:=^72}\n'.format(' RESULTS '))
    a1= hmat_orig[0:3,0]
    a2= hmat_orig[0:3,1]
    a3= hmat_orig[0:3,2]
    uvol= np.dot(a1,np.cross(a2,a3))
    lc= (plsq[0][2]/uvol)**(1.0/3)
    print(' Lattice constant = {0:10.4f} Ang.'.format(lc))
    print(' Cohesive energy  = {0:10.3f} eV'.format(plsq[0][3]/natm))
    print(' Bulk modulus     = {0:10.2f} GPa'.format(plsq[0][0]*1.602e+2))
    logfile.write(' Lattice constant = {0:10.4f} Ang.\n'.format(lc))
    logfile.write(' Cohesive energy  = {0:10.3f} eV\n'.format(plsq[0][3]/natm))
    logfile.write(' Bulk modulus     = {0:10.2f} GPa\n'.format(plsq[0][0]*1.602e+2))
    if show_graph:
        plt.plot(xarr,peval(xarr,plsq[0]),xarr,yarr,'o')
        plt.title('Data fitted with Murnaghan eq.')
        plt.legend(['fitted','data'])
        plt.xlabel('Volume (Ang.^3)')
        plt.ylabel('Energy (eV)')
        plt.savefig('graph.energy_vs_size.eps',dpi=150)
        plt.show()

    print('{0:=^72}'.format(' OUTPUT '))
    print(' * out.energy_vs_size')
    print(' * log.energy_vs_size')
    print(' * graph.energy_vs_size.eps')
