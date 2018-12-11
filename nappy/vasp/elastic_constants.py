#!/opt/local/bin/python
"""
Calculate elastic constants, C11, C12, C44,
Young's modulus, poison's ratio, and shear modulus,
by static method which measures energy differences
w.r.t. given strains.

Usage:
  elastic_constants.py [options]

Options:
  -h, --help  Show this help message and exit.
  -n NITER    Number of points to be calculated. [default: 10]
  -d DLTMAX   Max deviation of finite difference. [default: 0.01]
  -p          Plot a graph on the screen.
  --cmd=CMD   VASP execution command. [default: \"~/bin/vasp > out.vasp\"]
"""
from __future__ import print_function

import sys,os,commands
import numpy as np
from docopt import docopt
from scipy.optimize import curve_fit
_no_pyplot=False
try:
    import matplotlib.pyplot as plt
except:
    _no_pyplot=True

#...constants
outfname='out.elastic-constants'
logfname='log.elastic-constants'
graphname='graph.elastic-constants.eps'

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

def replace_hmat(hmat,fname='POSCAR'):
    f=open(fname,'r')
    ini= f.readlines()
    f.close()
    g=open(fname,'w')
    for l in range(len(ini)):
        if l in (2,3,4): #...hmat lines
            g.write(' {0:15.7f}'.format(hmat[l-2,0]))
            g.write(' {0:15.7f}'.format(hmat[l-2,1]))
            g.write(' {0:15.7f}'.format(hmat[l-2,2]))
            g.write('\n')
        else:
            g.write(ini[l])
    g.close()

def quad_func(x,a,b):
    return a *x**2 +b

if __name__ == '__main__':
    
    # usage= '%prog [options]'
    # 
    # parser= optparse.OptionParser(usage=usage)
    # parser.add_option("-n",dest="niter",type="int",default=10,
    #                   help="Number of points to be calculated.")
    # parser.add_option("-d",dest="dltmax",type="float",default=0.01,
    #                   help="Max deviation of finite difference..")
    # parser.add_option("-p",action="store_true",
    #                   dest="plot",default=False,
    #                   help="Plot a graph on the screen.")
    # parser.add_option("--cmd",dest="cmd",type="string",
    #                   default='~/bin/vasp > out.vasp',
    #                   help="vasp execution command")
    # (options,args)= parser.parse_args()
    # 
    # if len(args) != 0:
    #     print ' [Error] number of arguments wrong !!!'
    #     print '  Usage: $ {0}'.format(args[0])
    #     sys.exit()
    # 
    # niter= options.niter
    # shows_graph= options.plot
    # cmd= options.cmd
    # dltmax= options.dltmax

    args= docopt(__doc__)

    niter= int(args['-n'])
    dltmax= float['-d']
    shows_graph= args['-p']
    cmd= args['--cmd']

    al,hmat0,natm= read_POSCAR()
    hmax= np.max(hmat0)

    outfile1= open(outfname,'w')
    logfile= open(logfname,'w')
    #...get reference energy
    os.system(cmd)
    #os.system('mpirun -np 4 vasp > out.vasp')
    erg0= float(commands.getoutput("tail -n1 OSZICAR | awk '{print $5}'"))
    ncalc = 0
    dname="elastic-{0:05d}".format(ncalc)
    os.system("mkdir -p "+dname)
    os.system("cp INCAR OSZICAR OUTCAR vasprun.xml {0}/".format(dname))
    #erg0= float(commands.getoutput("grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'"))
    print(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(0.0,erg0,erg0,erg0))
    outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(0.0,erg0,erg0,erg0))
    logfile.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(0.0,erg0,erg0,erg0))
    ddlt= dltmax/niter
    #for iter in range(-niter/2,niter/2+1):
    for iter in range(1,niter+1):
        #dlt= (ddlt*(iter+1))
        dlt= ddlt*iter
        dh= hmax*dlt
        #...uniaxial strain for calc C11
        hmat= np.copy(hmat0)
        hmat[0,0]= hmat[0,0] +dh
        replace_hmat(hmat)
        os.system(cmd)
        #erg11= float(commands.getoutput("grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'"))
        erg11= float(commands.getoutput("tail -n1 OSZICAR | awk '{print $5}'"))
        ncalc += 1
        dname="elastic-{0:05d}".format(ncalc)
        os.system("mkdir -p "+dname)
        os.system("cp INCAR OSZICAR OUTCAR vasprun.xml {0}/".format(dname))

        #...orthorhombic volume-conserving strain for (C11-C12)
        hmat= np.copy(hmat0)
        hmat[0,0]= hmat[0,0] +dh
        hmat[1,1]= hmat[1,1] -dh
        hmat[2,2]= hmat[2,2] +dh**2/(1.0-dh**2)
        replace_hmat(hmat)
        os.system(cmd)
        #erg12= float(commands.getoutput("grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'"))
        erg12= float(commands.getoutput("tail -n1 OSZICAR | awk '{print $5}'"))
        ncalc += 1
        dname="elastic-{0:05d}".format(ncalc)
        os.system("mkdir -p "+dname)
        os.system("cp INCAR OSZICAR OUTCAR vasprun.xml {0}/".format(dname))

        #...monoclinic volume-conserving strain for C44
        hmat= np.copy(hmat0)
        hmat[0,1]= hmat[0,1] +dh/2
        hmat[1,0]= hmat[1,0] +dh/2
        hmat[2,2]= hmat[2,2] +dh**2/(4.0-dh**2)
        replace_hmat(hmat)
        os.system(cmd)
        #erg44= float(commands.getoutput("grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'"))
        erg44= float(commands.getoutput("tail -n1 OSZICAR | awk '{print $5}'"))
        ncalc += 1
        dname="elastic-{0:05d}".format(ncalc)
        os.system("mkdir -p "+dname)
        os.system("cp INCAR OSZICAR OUTCAR vasprun.xml {0}/".format(dname))

        print(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}'.format(dlt,erg11,erg12,erg44))
        outfile1.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(dlt,erg11,erg12,erg44))
        logfile.write(' {0:10.4f} {1:15.7f} {2:15.7f} {3:15.7f}\n'.format(dlt,erg11,erg12,erg44))
    outfile1.close()

    #...revert POSCAR
    replace_hmat(hmat0)

    #...prepare for Murnaghan fitting
    f= open(outfname,'r')
    lines= f.readlines()
    dlts= np.zeros((len(lines)))
    e11s= np.zeros((len(lines)))
    e12s= np.zeros((len(lines)))
    e44s= np.zeros((len(lines)))
    for l in range(len(lines)):
        dat= lines[l].split()
        dlts[l]= float(dat[0])
        e11s[l]= float(dat[1])
        e12s[l]= float(dat[2])
        e44s[l]= float(dat[3])
    f.close()

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
    ymod= c44*(2.0*c44+3.0*c12)/(c11+c44)
    prto= c12/2.0/(c11+c44)
    smod= ymod/2.0/(1.0+prto)
    str= '{0:=^72}\n'.format(' RESULTS ') \
         +' C11     = {0:10.3f} GPa\n'.format(c11) \
         +' C11-C12 = {0:10.3f} GPa\n'.format(c11_c12) \
         +' C12     = {0:10.3f} GPa\n'.format(c12) \
         +' C44     = {0:10.3f} GPa\n'.format(c44) \
         +' Young\'s modulus = {0:10.3f} GPa\n'.format(ymod) \
         +' shear modulus   = {0:10.3f} GPa\n'.format(smod) \
         +' Poisson\'s ratio = {0:10.3f}\n'.format(prto)
    print(str)
    logfile.write(str)
    logfile.close()

    if shows_graph:
        plt.plot(dlts,quad_func(dlts,*popt11),dlts,e11s,'o')
        plt.plot(dlts,quad_func(dlts,*popt12),dlts,e12s,'o')
        plt.plot(dlts,quad_func(dlts,*popt44),dlts,e44s,'o')
        plt.title('Energy vs. strain')
        plt.legend(['C11 fitted','C11 data'
                    ,'C12 fitted','C12 data'
                    ,'C44 fitted','C44 data'],loc=2)
        plt.xlabel('Strain')
        plt.ylabel('Energy (eV)')
        plt.savefig(graphname,dpi=150)
        plt.show()

    print('{0:=^72}'.format(' OUTPUT '))
    print(' * '+outfname)
    print(' * '+logfname)
    print(' * '+graphname)
