#!/bin/env python
"""
Smear 1D data of equidistance intervals with Gaussians.

Usage:
    gaussian_smear.py [options] INFILE

Options:
    -h, --help  Show this help messange and exit.
    -s, --sigma=SIGMA
                Width of Gaussian (sigma) in the unit of dx. [default: 2.0]
    -x=<x>      Position of x data (start from 1) [default: 1]
    -y=<y>      Position of y data (start from 1) [default: 1]
"""
from __future__ import print_function

import os,sys,math
import numpy as np
from docopt import docopt

def gsmear(xd,yd,sigma,ngwidth=3):
    """
    Compute Gaussian-smearing of given data with smearing width *sigma*.
    """
    dx= xd[1]-xd[0]
    sgm= sigma*dx
    pref= dx/(np.sqrt(2.0*np.pi)*sgm)
    ndat= len(xd)
    nwidth = int(sigma*ngwidth)
    #nwidth = ndat
    gdat= np.zeros(ndat,dtype=np.float)
    expfact = 1.0/sgm**2/2
    for ix in range(ndat):
        #gdat[ix]= yd[ix]
        for jx in range(-nwidth+1,nwidth-1):
            kx= ix+jx
            if kx < 0:
                kx = -kx
            elif kx >= ndat:
                kx = ndat -(kx-(ndat-1))
            gdat[ix] += yd[kx]*pref*np.exp(-(dx*(jx))**2*expfact)
    gdat /= 2
    return gdat

def gsmear_file(infname,sigma,xp,yp):
    """
    Read data from the given file and pass data to `gsmear` method.
    """
    infile= open(infname,'r')
    nline= 0
    for line in infile.readlines():
        if line[0] == "#":
            continue
        nline += 1
    infile.seek(0)
    xd= np.zeros(nline,dtype=np.float)
    yd= np.zeros(nline,dtype=np.float)
    il= 0
    for line in infile.readlines():
        if line[0] == "#":
            continue
        sline= line.split()
        xd[il]= float(sline[xp])
        yd[il]= float(sline[yp])
        il += 1
    infile.close()
    xd,gdat= gsmear(xd,yd,sigma)
    return xd,gdat
    

if __name__ == "__main__":

    args= docopt(__doc__)
    xp= int(args['-x']) -1
    yp= int(args['-y']) -1
    sigma= float(args['--sigma'])
    infname= args['INFILE']

    xd,gdat= gsmear_file(infname,sigma,xp,yp)
    
    outfile= open(infname+'.smeared','w')
    for il in range(len(gdat)):
        outfile.write(' {0:15.2f} {1:15.7f}\n'.format(xd[il],gdat[il]))
    outfile.close()
    
    print(' write '+infname+'.smeared')
    print(' GAUSSIAN_SMEAR done ')
