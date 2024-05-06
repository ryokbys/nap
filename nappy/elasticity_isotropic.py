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
"""
import sys,os
import subprocess
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit
import copy
#from nappy.napsys import NAPSystem
from nappy.io import read

#...constants
outfname='out.elast_iso'

def quad_func(x,a,b):
    return a *x**2 +b


if __name__ == '__main__':
    
    args = docopt(__doc__)
    niter = int(args['-n'])
    dltmax = float(args['-d'])
    mdexec = args['--mdexec']

    if niter < 2:
        raise ValueError('NITER {0:d} should be larger than 1.'.format(niter))

    if niter % 2 == 0:
        niter += 1

    #nsys0 = NAPSystem(fname='pmdini')
    nsys0 = read('pmdini')
    os.system('cp pmdini pmdini.bak')
    hmat0 = nsys0.get_hmat()
    print('Original H-matrix:')
    print(hmat0)
    # al,hmat0,natm= read_pmd()
    hmax= np.max(hmat0)
    print('Maximum in H-matrix elements = ',hmax)

    outfile1= open(outfname,'w')
    #...get reference energy
    os.system(mdexec+' > out.pmd')
    erg0 = float(subprocess.check_output("cat erg.pmd",shell=True))
    #erg0= float(commands.getoutput("grep 'potential energy' out.pmd | tail -n1 | awk '{print $3}'"))
    # print ' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(0.0,erg0,erg0,erg0)
    # outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(0.0,erg0,erg0,erg0))
    dltmin = -dltmax
    ddlt= (dltmax-dltmin)/(niter-1)
    dlts = np.zeros(niter,dtype=float)
    e11s = np.zeros(niter,dtype=float)
    e12s = np.zeros(niter,dtype=float)
    e44s = np.zeros(niter,dtype=float)
    
    #for iter in range(-niter/2,niter/2+1):
    for it in range(niter):
        nsys = copy.deepcopy(nsys0)
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
            nsys.set_hmat(hmat)
            write(nsys,fname='pmdini')
            os.system(mdexec+' > out.pmd')
            erg11 = float(subprocess.check_output('cat erg.pmd',shell=True))
            dlts[it] = dlt
            e11s[it] = erg11
    
            #...orthorhombic volume-conserving strain for (C11-C12)
            hmat= np.copy(hmat0)
            hmat[0,0]= hmat[0,0] +dh
            hmat[1,1]= hmat[1,1] -dh
            hmat[2,2]= hmat[2,2] +dh**2/(1.0-dh**2)
            nsys.set_hmat(hmat)
            nsys.write('pmdini')
            os.system(mdexec+' > out.pmd')
            erg12 = float(subprocess.check_output('cat erg.pmd',shell=True))
            e12s[it] = erg12
    
            #...monoclinic volume-conserving strain for C44
            hmat= np.copy(hmat0)
            hmat[0,1]= hmat[0,1] +dh/2
            hmat[1,0]= hmat[1,0] +dh/2
            hmat[2,2]= hmat[2,2] +dh**2/(4.0-dh**2)
            nsys.set_hmat(hmat)
            nsys.write('pmdini')
            os.system(mdexec+' > out.pmd')
            erg44 = float(subprocess.check_output('cat erg.pmd',shell=True))
            e44s[it] = erg44

        print(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(dlts[it],e11s[it],e12s[it],e44s[it]))
        outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(dlts[it],e11s[it],e12s[it],e44s[it]))
    outfile1.close()

    #...revert pmdini
    os.system('mv pmdini.bak pmdini')

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
    # vol= get_vol(al,hmat0)
    vol = nsys0.get_volume()
    #...least square fitting
    popt11,pcov11= curve_fit(quad_func,dlts,e11s,p0=p0)
    popt12,pcov12= curve_fit(quad_func,dlts,e12s,p0=p0)
    popt44,pcov44= curve_fit(quad_func,dlts,e44s,p0=p0)

    c11= popt11[0]/vol*2 *160.218
    c11_c12= popt12[0]/vol *160.218
    c12= c11 -c11_c12
    c44= popt44[0]/vol*2 *160.218

    cij = np.zeros((6,6),dtype=float)
    cij[0,0] = cij[1,1] = cij[2,2] = c11
    cij[0,1] = cij[0,2] = cij[1,2] = cij[2,1] = cij[2,0] = cij[1,0] = c12
    cij[3,3] = cij[4,4] = cij[5,5] = c44
    sij = np.linalg.inv(cij)

    #...output results
    print('{0:=^72}'.format(' RESULTS '))

    print(' Cij [GPa]:')
    for i in range(6):
        for j in range(6):
            print(' {0:8.2f}'.format(cij[i,j]),end='')
        print('')

    c112233 = cij[0,0]+cij[1,1]+cij[2,2]
    c122331 = cij[0,1]+cij[0,2]+cij[1,2]
    c445566 = cij[3,3]+cij[4,4]+cij[5,5]
    s112233 = sij[0,0]+sij[1,1]+sij[2,2]
    s122331 = sij[0,1]+sij[0,2]+sij[1,2]
    s445566 = sij[3,3]+sij[4,4]+sij[5,5]
    kv = (c112233 +2.0*c122331)/9
    kr = 1.0/(s112233 +2.0*(s122331))
    gv = (c112233 -c122331 +3.0*c445566)/15
    gr = 15.0 /(4.0*s112233 -4.0*s122331 +3.0*s445566)
    kvrh = (kv+kr)/2
    gvrh = (gv+gr)/2
    prto2 = (3.0*kvrh -2.0*gvrh)/(6.0*kvrh +2.0*gvrh)

    print('')
    # print ' Definition of the following values, see ' \
    #     +'https://materialsproject.org/wiki/index.php/Elasticity_calculations'
    # print ' K_V   = {0:10.3f} GPa'.format(kv)
    # print ' K_R   = {0:10.3f} GPa'.format(kr)
    # print ' G_V   = {0:10.3f} GPa'.format(gr)
    # print ' G_R   = {0:10.3f} GPa'.format(gv)
    print(' Bulk modulus    = {0:10.3f} GPa'.format(kvrh))
    print(' shear modulus   = {0:10.3f} GPa'.format(gvrh))
    print(' Poisson\'s ratio = {0:10.3f}'.format(prto2))
    

    print('{0:=^72}'.format(' OUTPUT '))
    print(' * '+outfname)
