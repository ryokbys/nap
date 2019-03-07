#!/usr/bin/env python
"""
Make graphs from force_cspline data points.

Usage:
  cspline2graph.py [options] FILE

Options:
  -h, --help  Show this message and exit.
  --figsize FIGSIZE
              Size of figure, pass to matplotlib. [default: 15,10]
  --ylim YLIM
              Limits of y range. [default: -4.0,4.0]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
from scipy.interpolate import CubicSpline

__author__ = "RYO KOBAYASHI"
__version__ = ""

class CSpline():

    def __init__(self):
        self.cstype = None
        self.comb = None
        self.rcut = 0.0
        self.npnt = 0
        self.xdata = []
        self.ydata = []
        self.cspl = None
        
    def make_cspl(self):
        """
        Make coefficients from data that should be loaded beforehand.
        """
        if 'angular' in self.cstype: # 3body
            self.cspl = CubicSpline(self.xdata,self.ydata,)
        else: # 2body, left edge: natural, right edge: clamped
            self.cspl = CubicSpline(self.xdata,self.ydata,bc_type=((2,0.0),(1,0.0)))

    def curve(self,x):
        """
        Return y at given x points.
        """
        if self.cspl is None:
            self.make_cspl()
        y = self.cspl(x)
        return y

def read_params_cspline(fname):
    """
    Read in.params.cspline and return list of CSplines.
    """
    ncsp = 0
    csps = []
    with open(fname,'r') as f:
        lines = f.readlines()
    for il,line in enumerate(lines):
        data = line.split()
        if not len(data) > 0: # blank line
            continue
        if data[0][0] in ('#', '!'):  # comment line
            continue
        if ncsp == 0:
            ncsp = int(data[0])
            continue
        if 'radial' in data[0] or 'angular' in data[0]:
            csp = CSpline()
            csp.cstype = data[0]
            if data[0] == 'radial':
                isp = int(data[1])
                jsp = int(data[2])
                comb = (isp,jsp)
                rcut = float(data[3])
                npnt = int(data[4])
            else:  # angulars
                isp = int(data[1])
                jsp = int(data[2])
                ksp = int(data[3])
                comb = (isp,jsp,ksp)
                rcut = float(data[4])
                npnt = int(data[5])
            csp.comb = comb
            csp.rcut = rcut
            csp.npnt = npnt
        else:
            xi,yi = [ float(x) for x in data ]
            csp.xdata.append(xi)
            csp.ydata.append(yi)
            if len(csp.xdata) == csp.npnt:
                csps.append(csp)
                csp = None
    if len(csps) != ncsp:
        raise ValueError('ncsp != len(csps)')
    return csps

def get_num_species(csps):
    max_spcs = 0
    for c in csps:
        comb = c.comb
        max_spcs = max(max_spcs,*comb)
    return max_spcs

def get_rcut2b(csps):
    rcut = 0.0
    for c in csps:
        if c.cstype == 'radial':
            rcut = max(c.rcut,rcut)
    return rcut

def get_rmin2b(csps):
    rmin = 100000.0
    for c in csps:
        if c.cstype == 'radial':
            rmin = min(c.xdata[0],rmin)
    return rmin

if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['FILE']
    figsize = [ int(x) for x in args['--figsize'].split(',') ]
    ylim = [ float(x) for x in args['--ylim'].split(',') ]

    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
    except:
        raise ImportError('Failed to import matplotlib and seaborn.')

    csps = read_params_cspline(fname)
    nsp = get_num_species(csps)
    rcut = get_rcut2b(csps)
    rmin = get_rmin2b(csps)
    print(' nsp = ',nsp)
    print(' rmin,rcut = ',rmin,rcut)

    #...Plot 2-body potentials
    fig,axes = plt.subplots(nsp,nsp,figsize=figsize,sharex=True,sharey=True)
    for i in range(nsp):
        isp = i +1
        for j in range(nsp):
            jsp = j +1
            if jsp < isp:
                axes[i,j].axis('off')
                continue
            for c in csps:
                if c.cstype == 'radial' and (isp,jsp) == c.comb:
                    xs = np.linspace(c.xdata[0],rcut,200)
                    ys = c.curve(xs)
                    xd = c.xdata
                    yd = c.ydata
                    break
            axes[i,j].plot(xs,ys,'r-')
            axes[i,j].plot(xd,yd,'ro')
            axes[i,j].set_title('{0:d}-{1:d}'.format(isp,jsp))
            axes[i,j].set_ylim(*ylim)
            axes[i,j].set_xlim(rmin,rcut)
    plt.xlabel('Distance (Ang.)')
    plt.ylabel('Energy (eV)')
    plt.savefig("graph_cspline2b.png", format='png', dpi=300, bbox_inches='tight')
    print('Wrote graph_cspline2b.png')

    #...Plot 3-body potentials
    xs = np.linspace(-1.0, 1.0, 200)
    for c in csps:
        if not 'angular' in c.cstype:
            continue
        isp,jsp,ksp = c.comb
        ys = c.curve(xs)
        plt.clf()
        plt.figure(figsize=(8,6))
        plt.plot(xs,ys,'r-')
        plt.plot(c.xdata,c.ydata,'ro')
        comb = '{0:d}-{1:d}-{2:d}'.format(jsp,isp,ksp)
        plt.title(comb)
        plt.xlabel('cos({0:s})'.format(comb))
        plt.ylabel('Energy (eV)')
        plt.savefig("graph_cspline3b_{0:s}.png".format(comb), format='png', dpi=300, bbox_inches='tight')
        print('Wrote graph_cspline3b_{0:s}.png'.format(comb))
    
