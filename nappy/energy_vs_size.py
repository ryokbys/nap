#!/opt/local/bin/python
"""
Calculate the cohesive energy as a function of lattice constant,
by altering the lattice constant in pmdini file.

And if possible, calculate equilibrium lattice size and
bulk modulus, too.

MIN and MAX values correspond to the original value in *pmdini* file.
So check the value at the 2nd line in *pmdini* file and usually MIN and MAX
values should smaller and larger than the original value, respectively.

Usage:
  energy_vs_size.py [options] MIN MAX

Options:
  -h, --help  Show this message and exit.
  -n ITER     Num of points to be calculated. [default: 10]
  --mdexec=MDEXEC
              Path to *pmd*. [default: ~/src/nap/pmd/pmd]
"""

import sys,os,commands
import numpy as np
from docopt import docopt
from scipy.optimize import leastsq


def read_pmd(fname="pmdini"):
    f=open(fname,'r')
    #...read 1st line and get current lattice size
    al= float(f.readline().split()[0])
    hmat= np.zeros((3,3))
    hmat[0]= [ float(x) for x in f.readline().split() ]
    hmat[1]= [ float(x) for x in f.readline().split() ]
    hmat[2]= [ float(x) for x in f.readline().split() ]
    f.readline()
    f.readline()
    f.readline()
    natm= int(f.readline().split()[0])
    f.close()
    return (al,hmat,natm)

def get_vol(al,hmat):
    a1= hmat[0:3,0] *al
    a2= hmat[0:3,1] *al
    a3= hmat[0:3,2] *al
    return np.dot(a1,np.cross(a2,a3))

def replace_1st_line(x,fname="pmdini"):
    f=open(fname,'r')
    ini= f.readlines()
    f.close()
    g=open(fname,'w')
    for l in range(len(ini)):
        if l == 0:
            g.write(' {0:10.4f}\n'.format(x))
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

    args = docopt(__doc__)

    niter = int(args['-n'])
    mdexec = args['--mdexec']

    fname = 'pmdini'
    tmpfname = fname+'.tmp'
    os.system('cp '+fname+' '+tmpfname)
    al_orig,hmat,natm= read_pmd(fname)
    al_min = float(args['MIN'])
    al_max = float(args['MAX'])

    if al_orig < al_min or al_orig > al_max:
        print ' [Warning] min and max maybe wrong.'
        print '   I hope you know what you are doing.'
        print '   al_min, al_orig, al_max=',al_min, al_orig, al_max
        #sys.exit()

    logfile= open('log.energy_vs_size','w')
    outfile1= open('out.energy_vs_size','w')
    dl= (al_max -al_min)/niter
    for iter in range(niter+1):
        al= al_min +dl*iter
        replace_1st_line(al,fname)
        os.system(mdexec +' > out.pmd')
        erg= float(commands.getoutput("grep 'Potential energy' out.pmd | tail -n1 | awk '{print $3}'"))
        vol= get_vol(al,hmat)
        print ' {0:10.4f} {1:10.4f} {2:15.7f}'.format(al,vol,erg)
        outfile1.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
        logfile.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
    outfile1.close()

    #...revert pmdini
    os.system('cp '+tmpfname+' '+fname)
    os.system('rm '+tmpfname)
    #replace_1st_line(al_orig,fname)

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
    print ' plsq=',plsq[0]
    print '{0:=^72}'.format(' RESULTS ')
    logfile.write('{0:=^72}\n'.format(' RESULTS '))
    a1= hmat[0:3,0]
    a2= hmat[0:3,1]
    a3= hmat[0:3,2]
    uvol= np.dot(a1,np.cross(a2,a3))
    lc= (plsq[0][2]/uvol)**(1.0/3)
    print ' Lattice constant = {0:10.4f} Ang.'.format(lc)
    print ' Cohesive energy  = {0:10.3f} eV'.format(plsq[0][3]/natm)
    print ' Bulk modulus     = {0:10.2f} GPa'.format(plsq[0][0]*1.602e+2)
    logfile.write(' Lattice constant = {0:10.4f} Ang.\n'.format(lc))
    logfile.write(' Cohesive energy  = {0:10.3f} eV\n'.format(plsq[0][3]/natm))
    logfile.write(' Bulk modulus     = {0:10.2f} GPa\n'.format(plsq[0][0]*1.602e+2))

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * out.energy_vs_size'
    print ' * log.energy_vs_size'
    #print ' * graph.Ecoh-vs-size.eps'
