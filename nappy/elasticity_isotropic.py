#!/opt/local/bin/python
"""
Calculate elastic constants, C11, C12, C44,
Young's modulus, poison's ratio, and shear modulus,
by static method which measures energy differences
w.r.t. given strains.

Usage:
  elastic_constants.py [options]

Options:
  -h, --help  Show this message and exit.
  -n NITER    Num of points to be calculated. 
              Even number is same as an odd number NITER+1. [default: 5]
  -d DLTMAX   Max deviation of finite difference. [default: 0.01]
  --mdexec=MDEXEC
              Path to *pmd*. [default: ~/src/nap/pmd/pmd]
  --relax     Relax internal coordinates. [default: False]
"""

import sys,os,commands
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit

#...constants
outfname='out.elastic_constants'
logfname='log.elastic_constants'

def read_pmd():
    f=open('pmdini','r')
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

def replace_hmat(hmat):
    f=open('pmdini','r')
    ini= f.readlines()
    f.close()
    g=open('pmdini','w')
    for l in range(len(ini)):
        if l in (1,2,3): #...hmat lines
            g.write(' {0:15.7f}'.format(hmat[l-1,0]))
            g.write(' {0:15.7f}'.format(hmat[l-1,1]))
            g.write(' {0:15.7f}'.format(hmat[l-1,2]))
            g.write('\n')
        else:
            g.write(ini[l])
    g.close()

def quad_func(x,a,b):
    return a *x**2 +b

if __name__ == '__main__':
    
    args = docopt(__doc__)
    niter = int(args['-n'])
    dltmax = float(args['-d'])
    mdexec = args['--mdexec']
    relax = args['--relax']

    if niter < 2:
        raise ValueError('NITER {0:d} should be larger than 1.'.format(niter))

    if niter % 2 == 0:
        niter += 1

    al,hmat0,natm= read_pmd()
    hmax= np.max(hmat0)

    logfile= open(logfname,'w')
    outfile1= open(outfname,'w')
    #...get reference energy
    os.system(mdexec+' > out.pmd')
    erg0= float(commands.getoutput("grep 'potential energy' out.pmd | tail -n1 | awk '{print $3}'"))
    # print ' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(0.0,erg0,erg0,erg0)
    # outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(0.0,erg0,erg0,erg0))
    # logfile.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(0.0,erg0,erg0,erg0))
    dltmin = -dltmax
    ddlt= (dltmax-dltmin)/(niter-1)
    dlts = np.zeros(niter,dtype=float)
    e11s = np.zeros(niter,dtype=float)
    e12s = np.zeros(niter,dtype=float)
    e44s = np.zeros(niter,dtype=float)
    
    #for iter in range(-niter/2,niter/2+1):
    for it in range(niter):
        #dlt= (ddlt*(iter+1))
        if it == niter/2:
            dlts[it] = 0.0
            e11s[it] = erg0
            e12s[it] = erg0
            e44s[it] = erg0
        else:
            dlt= dltmin +ddlt*it
            dh= hmax*dlt
            #...uniaxial strain for calc C11
            hmat= np.copy(hmat0)
            hmat[0,0]= hmat[0,0] +dh
            replace_hmat(hmat)
            os.system(mdexec+' > out.pmd')
            erg11= float(commands.getoutput("grep 'potential energy' out.pmd | tail -n1 | awk '{print $3}'"))
            dlts[it] = dlt
            e11s[it] = erg11
    
            #...orthorhombic volume-conserving strain for (C11-C12)
            hmat= np.copy(hmat0)
            hmat[0,0]= hmat[0,0] +dh
            hmat[1,1]= hmat[1,1] -dh
            hmat[2,2]= hmat[2,2] +dh**2/(1.0-dh**2)
            replace_hmat(hmat)
            os.system(mdexec+' > out.pmd')
            erg12= float(commands.getoutput("grep 'potential energy' out.pmd | tail -n1 | awk '{print $3}'"))
            e12s[it] = erg12
    
            #...monoclinic volume-conserving strain for C44
            hmat= np.copy(hmat0)
            hmat[0,1]= hmat[0,1] +dh/2
            hmat[1,0]= hmat[1,0] +dh/2
            hmat[2,2]= hmat[2,2] +dh**2/(4.0-dh**2)
            replace_hmat(hmat)
            os.system(mdexec+' > out.pmd')
            erg44= float(commands.getoutput("grep 'potential energy' out.pmd | tail -n1 | awk '{print $3}'"))        
            e44s[it] = erg44

        print ' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(dlts[it],e11s[it],e12s[it],e44s[it])
        logfile.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(dlts[it],e11s[it],e12s[it],e44s[it]))
        outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(dlts[it],e11s[it],e12s[it],e44s[it]))
    outfile1.close()

    #...revert pmdini
    replace_hmat(hmat0)

    #...prepare for Murnaghan fitting
    # f= open(outfname,'r')
    # lines= f.readlines()
    # dlts= np.zeros((len(lines)))
    # e11s= np.zeros((len(lines)))
    # e12s= np.zeros((len(lines)))
    # e44s= np.zeros((len(lines)))
    # for l in range(len(lines)):
    #     dat= lines[l].split()
    #     dlts[l]= float(dat[0])
    #     e11s[l]= float(dat[1])
    #     e12s[l]= float(dat[2])
    #     e44s[l]= float(dat[3])
    # f.close()

    #...set initial values
    a= 1.0
    b= erg0
    p0= np.array([a,b])
    vol= get_vol(al,hmat0)
    #...least square fitting
    popt11,pcov11= curve_fit(quad_func,dlts,e11s,p0=p0)
    popt12,pcov12= curve_fit(quad_func,dlts,e12s,p0=p0)
    popt44,pcov44= curve_fit(quad_func,dlts,e44s,p0=p0)

    c11= popt11[0]/vol*2 *160.218
    c11_c12= popt12[0]/vol *160.218
    c12= c11 -c11_c12
    c44= popt44[0]/vol*2 *160.218

    #...output results
    print '{0:=^72}'.format(' RESULTS ')
    logfile.write('{0:=^72}\n'.format(' RESULTS '))
    print ' C11     = {0:10.3f} GPa'.format(c11)
    print ' C11-C12 = {0:10.3f} GPa'.format(c11_c12)
    print ' C12     = {0:10.3f} GPa'.format(c12)
    print ' C44     = {0:10.3f} GPa'.format(c44)
    logfile.write(' C11     = {0:10.3f} GPa\n'.format(c11))
    logfile.write(' C11-C12 = {0:10.3f} GPa\n'.format(c11_c12))
    logfile.write(' C12     = {0:10.3f} GPa\n'.format(c12))
    logfile.write(' C44     = {0:10.3f} GPa\n'.format(c44))
    ymod= c44*(2.0*c44+3.0*c12)/(c11+c44)
    prto= c12/2.0/(c11+c44)
    smod= ymod/2.0/(1.0+prto)
    print ' Following values maybe only valid for isotropic materials...'
    print ' Young\'s modulus = {0:10.3f} GPa'.format(ymod)
    print ' shear modulus   = {0:10.3f} GPa'.format(smod)
    print ' Poisson\'s ratio = {0:10.3f}'.format(prto)
    logfile.write(' Following values maybe only valid for isotropic materials...\n')
    logfile.write(' Young\'s modulus = {0:10.3f} GPa\n'.format(ymod))
    logfile.write(' shear modulus   = {0:10.3f} GPa\n'.format(smod))
    logfile.write(' Poisson\'s ratio = {0:10.3f}\n'.format(prto))
    logfile.close()

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * '+outfname
