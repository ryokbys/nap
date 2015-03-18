#!/opt/local/bin/python
u"""
Calculate the cohesive energy as a function of lattice constant,
by altering the lattice constant in pmd00000 file.

And if possible, calculate equilibrium lattice size and
bulk modulus, too.
"""

import sys,os,commands
import numpy as np
import optparse
from scipy.optimize import leastsq
import matplotlib.pyplot as plt


def read_pmd():
    f=open('0000/pmd00000','r')
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

def replace_1st_line(x):
    f=open('0000/pmd00000','r')
    ini= f.readlines()
    f.close()
    g=open('0000/pmd00000','w')
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
    
    usage= '%prog [options] <min> <max>'

    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-n",dest="niter",type="int",default=10,
                      help="Number of points to be calculated.")
    parser.add_option("-p",action="store_true",
                      dest="plot",default=False,
                      help="Plot a graph on the screen.")
    parser.add_option("--pmdexec",dest="pmdexec",type="string",
                      default='~/src/nap/pmd/pmd',
                      help="path to the pmd executable.")
    (options,args)= parser.parse_args()

    if len(args) != 2:
        print ' [Error] number of arguments wrong !!!'
        print usage
        sys.exit()

    niter= options.niter
    shows_graph= options.plot
    pmdexec= options.pmdexec

    al_orig,hmat,natm= read_pmd()
    al_min = float(sys.argv[1])
    al_max = float(sys.argv[2])

    if al_orig < al_min or al_orig > al_max:
        print ' [Warning] min and max maybe wrong.'
        print '   hoping you know what you are doing.'
        print '   al_min, al_orig, al_max=',al_min, al_orig, al_max
        #sys.exit()

    logfile= open('log.Ecoh-vs-size','w')
    outfile1= open('out.Ecoh-vs-size','w')
    dl= (al_max -al_min)/niter
    for iter in range(niter+1):
        al= al_min +dl*iter
        replace_1st_line(al)
        os.system(pmdexec +' > out.pmd')
        erg= float(commands.getoutput("grep 'potential energy' out.pmd | head -n1 | awk '{print $3}'"))
        vol= get_vol(al,hmat)
        print ' {0:10.4f} {1:10.4f} {2:15.7f}'.format(al,vol,erg)
        outfile1.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
        logfile.write(' {0:10.4f} {1:10.4f} {2:15.7f}\n'.format(al,vol,erg))
    outfile1.close()

    #...revert 0000/pmd00000
    replace_1st_line(al_orig)

    #...prepare for Murnaghan fitting
    f= open('out.Ecoh-vs-size','r')
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
    if shows_graph:
        plt.plot(xarr,peval(xarr,plsq[0]),xarr,yarr,'o')
        plt.title('Data fitted with Murnaghan eq.')
        plt.legend(['fitted','data'])
        plt.xlabel('Volume (Ang.^3)')
        plt.ylabel('Energy (eV)')
        plt.savefig('graph.Ecoh-vs-size.eps',dpi=150)
        plt.show()

    print '{0:=^72}'.format(' OUTPUT ')
    print ' * out.Ecoh-vs-size'
    print ' * log.Ecoh-vs-size'
    print ' * graph.Ecoh-vs-size.eps'
