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
              Path to *pmd*. [default: pmd]
"""
from __future__ import print_function

import sys
import os
import subprocess
import numpy as np
from docopt import docopt
from scipy.optimize import leastsq

_infname = 'pmdini'

def read_pmd(fname="pmdini"):
    with open(fname,'r') as f:
        lines = f.readlines()
    nl = 0
    hmat= np.zeros((3,3))
    for line in lines:
        if line[0] in ('#','!'):
            continue
        nl += 1
        #...read 1st line and get current lattice size
        if( nl==1 ):
            al= float(line.split()[0])
        elif( nl==2 ):
            hmat[0]= [ float(x) for x in line.split() ]
        elif( nl==3 ):
            hmat[1]= [ float(x) for x in line.split() ]
        elif( nl==4 ):
            hmat[2]= [ float(x) for x in line.split() ]
        elif( nl in (5,6,7) ):
            continue
        elif( nl==8 ):
            natm= int(line.split()[0])
        else:
            break

    return (al,hmat,natm)

def get_vol(al,hmat):
    a1= hmat[0:3,0] *al
    a2= hmat[0:3,1] *al
    a3= hmat[0:3,2] *al
    return np.dot(a1,np.cross(a2,a3))

def replace_1st_line(x,fname="pmdini"):
    with open(fname,'r') as f:
        orig_txt = f.readlines()

    with open(fname,'w') as g:
        nl = 0
        for line in orig_txt:
            if line[0] in ('#','!'):
                g.write(line)
            else:
                nl += 1
                if nl == 1:
                    g.write(' {0:10.4f}\n'.format(x))
                else:
                    g.write(line)
    return None

def residuals(p,y,x):
    b,bp,v0,ev0= p
    err= y -( b*x/(bp*(bp-1.0)) *(bp*(1.0-v0/x) +(v0/x)**bp -1.0) +ev0 )
    return err

def peval(x,p):
    b,bp,v0,ev0= p
    return b*x/(bp*(bp-1.0)) *(bp*(1.0-v0/x) +(v0/x)**bp -1.0) +ev0

def erg_vs_size(fname,al_min,al_max,niter,mdexec):
    tmpfname = fname+'.tmp'
    os.system('cp '+fname+' '+tmpfname)
    al_orig,hmat,natm= read_pmd(fname)

    if al_orig < al_min or al_orig > al_max:
        print(' [Warning] min and max maybe wrong.')
        print('   I hope you know what you are doing.')
        print('   al_min, al_orig, al_max=',al_min, al_orig, al_max)
        #sys.exit()

    dl= (al_max -al_min)/(niter-1)

    als = np.zeros(niter)
    vols = np.zeros(niter)
    prss = np.zeros(niter)
    ergs = np.zeros(niter)

    print('     al        vol        erg             prss')
    for i in range(niter):
        al= al_min +dl*i
        replace_1st_line(al,fname)
        #os.system('rm -f out.pmd')
        #os.system(mdexec +' 2>&1 > out.pmd')
        pmdout = subprocess.check_output(mdexec)
        with open('out.pmd','w') as f:
            f.write(pmdout.decode('utf-8'))
        cmdstr = "grep 'Potential energy' out.pmd | tail -n1 | awk '{print $3}'"
        erg = float(subprocess.check_output(cmdstr,shell=True))
        cmdstr = "grep 'Pressure' out.pmd | tail -n1 | awk '{print $3}'"
        prs = float(subprocess.check_output(cmdstr,shell=True))
        vol= get_vol(al,hmat)
        als[i] = al
        vols[i] = vol
        ergs[i] = erg
        prss[i] = prs
        text = ' {0:10.4f} {1:10.4f} {2:15.7f} {3:12.4f}'.format(al,vol,erg,prs)
        print(text)
        outfname = 'out.pmd.{0:03d}'.format(i)
        os.system('rm -f {0:s} && mv out.pmd {0:s}'.format(outfname))
        dumpfname = 'dump_0.{0:03d}'.format(i)
        os.system('rm -f {0:s} && mv dump_0 {0:s}'.format(dumpfname))
    #...revert pmdini
    os.system('rm -f '+fname)
    os.system('cp '+tmpfname+' '+fname)
    os.system('rm -f '+tmpfname)
    return als,vols,ergs,prss,al_orig,hmat,natm

def main():
    args = docopt(__doc__)

    niter = int(args['-n'])
    mdexec = args['--mdexec']
    al_min = float(args['MIN'])
    al_max = float(args['MAX'])

    als,vols,ergs,prss,al_orig,hmat,natm = erg_vs_size(_infname,al_min,al_max,niter,mdexec)

    logfile= open('log.energy_vs_size','w')
    outfile1= open('out.energy_vs_size','w')
    nprss = []
    text = '#    al,       volume,     energy,   press (pmd),   press (numerical)'
    outfile1.write(text+'\n')
    logfile.write(text+'\n')
    for i in range(len(als)):
        al = als[i]
        vol = vols[i]
        erg = ergs[i]
        prs = prss[i]
        if i == 0:
            nprs = -(ergs[i+1]-ergs[i])/(vols[i+1]-vols[i])
        elif i == len(vols)-1:
            nprs = -(ergs[i]-ergs[i-1])/(vols[i]-vols[i-1])
        else:
            nprs = -(ergs[i+1]-ergs[i-1])/(vols[i+1]-vols[i-1])
        # Convert eV/A^3 to GPa
        nprs *= 1.602e-19/(1.0e-10)**3 *1.0e-9
        nprss.append(nprs)
        text = ' {0:10.4f} {1:10.2f} {2:11.3f} {3:11.3f} {4:11.3f}'.format(al,vol,erg,prs,nprs)
        # print text
        outfile1.write(text+'\n')
        logfile.write(text+'\n')
    outfile1.close()

    #...set initial values
    b= 1.0
    bp= 2.0
    ev0= min(ergs)
    v0= vols[int(len(vols)/2)]
    p0= np.array([b,bp,v0,ev0])
    #...least square fitting
    plsq= leastsq(residuals,p0,args=(ergs,vols))

    #...output results
    print(' plsq=',plsq[0])
    print('{0:=^72}'.format(' RESULTS '))
    logfile.write('{0:=^72}\n'.format(' RESULTS '))
    a1= hmat[0:3,0]
    a2= hmat[0:3,1]
    a3= hmat[0:3,2]
    uvol= np.dot(a1,np.cross(a2,a3))
    lc= (plsq[0][2]/uvol)**(1.0/3)
    print(' Lattice constant = {0:10.4f} Ang.'.format(lc))
    print(' Cohesive energy  = {0:10.3f} eV'.format(plsq[0][3]/natm))
    print(' Bulk modulus     = {0:10.2f} GPa'.format(plsq[0][0]*1.602e+2))
    logfile.write(' Lattice constant = {0:10.4f} Ang.\n'.format(lc))
    logfile.write(' Cohesive energy  = {0:10.3f} eV\n'.format(plsq[0][3]/natm))
    logfile.write(' Bulk modulus     = {0:10.2f} GPa\n'.format(plsq[0][0]*1.602e+2))

    print('{0:=^72}'.format(' OUTPUT '))
    print(' * out.energy_vs_size')
    print(' * log.energy_vs_size')
    #print ' * graph.Ecoh-vs-size.eps'
    return None
    

if __name__ == '__main__':
    main()
