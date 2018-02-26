#!/usr/bin/env python
"""
Make graphs from fitpot results.

Usage:
  evaluate_result.py [options]

Options:
  -h,--help  Show this message and exit.
  --histogram
             Make histogram of 
  -w,--width=WIDTH
             Width of the bin in histogram in eV. [default: 0.0005]
"""
from __future__ import print_function

import os
from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt

__author__ = "RYO KOBAYASHI"
__version__ = "180224"

def make_fit_graph(key='erg'):

    try:
        import seaborn as sns
        sns.set(context='poster',style='ticks')
    except:
        pass

    fnametrn = 'out.{0:s}.trn.fin'.format(key)
    fnametst = 'out.{0:s}.tst.fin'.format(key)
    if not os.path.exists(fnametrn):
        raise IOError('File does not exist: '+fnametrn)
    if not os.path.exists(fnametst):
        raise IOError('File does not exist: '+fnametst)

    epottrn = []
    ereftrn = []
    with open(fnametrn,'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0] == '#':
            continue
        data = line.split()
        ereftrn.append(float(data[0]))
        epottrn.append(float(data[1]))

    epottst = []
    ereftst = []
    with open(fnametst,'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0] == '#':
            continue
        data = line.split()
        ereftst.append(float(data[0]))
        epottst.append(float(data[1]))

    emax = -1.0e+30
    emax = max(max(ereftrn),emax)
    emax = max(max(epottrn),emax)
    emin = 1.0e+30
    emin = min(min(ereftrn),emin)
    emin = min(min(epottrn),emin)
    if ereftst:
        emax = max(max(ereftst),emax)
        emin = min(min(ereftst),emin)
    if epottst:
        emax = max(max(epottst),emax)
        emin = min(max(epottst),emin)

    erange = [emin,emax]
    cmap = plt.get_cmap('tab10')
    makersize = 5
    size = 8
    plt.clf()
    plt.figure(figsize=(size,size))
    plt.plot(erange,erange,'--',color='black',linewidth=1.0)
    plt.plot(ereftrn,epottrn,'o',color=cmap(0),mec='white',mew=0.5,
             ms=makersize,label='training data')
    plt.plot(ereftst,epottst,'o',color=cmap(1),mec='white',mew=0.5,
             ms=makersize,label='test data')
    if key == 'erg':
        plt.xlabel('DFT energy (eV/atom)')
        plt.ylabel('Model-potential energy (eV/atom)')
    elif key == 'frc':
        plt.xlabel('DFT force (eV/Ang)')
        plt.ylabel('Model-potential force (eV/Ang)')
    elif key == 'strs':
        plt.xlabel('DFT stress (GPa)')
        plt.ylabel('Model-potential stress (GPa)')
    
    plt.legend(loc='best')
    fname = 'graph.{0:s}.png'.format(key)
    plt.savefig(fname,format='png',dpi=300,bbox_inches='tight')
    print('- {0:s}'.format(fname))
    return

def make_iter_graph():
    if not os.path.exists('out.fitpot'):
        raise RuntimeError('File not exist: out.fitpot')
    if not os.path.exists('out.iter') or \
       os.path.getmtime('out.iter') < os.path.getmtime('out.fitpot'):
        os.system('grep "iter,ftrn" out.fitpot > out.iter')

    try:
        import seaborn as sns
        sns.set(context='poster',style='ticks')
    except:
        pass

    with open('out.iter','r') as f:
        lines = f.readlines()

    iters = []
    ftrns = []
    ftsts = []
    for line in lines:
        if line[0] == '#':
            continue
        data = line.split()
        iters.append(int(data[1]))
        ftrns.append(float(data[2]))
        ftsts.append(float(data[3]))

    plot_test = True
    for f in ftsts:
        if f < 1.0e-8:
            plot_test = False

    #...Plotting
    plt.clf()
    cmap = plt.get_cmap('tab10')
    plt.plot(iters,ftrns,'-',color=cmap(0),label='traning data')
    if plot_test:
        plt.plot(iters,ftsts,'-',color=cmap(1),label='test data')
    plt.xlabel('Iteration')
    plt.ylabel('Loss function value')
    plt.yscale('log')
    plt.legend(loc='best')
    fname = 'graph.iter.png'
    plt.savefig(fname,format='png',dpi=300,bbox_inches='tight')
    print('- {0:s}'.format(fname))
    return


def histogram(fname,width=0.0005):

    if not os.path.exists(fname):
        raise IOError('File does not exist: '+fname)
    
    with open(fname,'r') as f:
        lines = f.readlines()
        data = np.zeros((len(lines),2),dtype=float)
        for il,l in enumerate(lines):
            d = l.split()
            data[il,:] = [ float(d[0]),float(d[1]) ]

    diffmax = 0.0
    for i in range(len(data)):
        diff = abs(data[i,0]-data[i,1])
        diffmax = max(diff,diffmax)
    ndim = int(diffmax/width)+1
    hg = np.zeros((ndim),dtype=int)
    for i in range(len(data)):
        diff = abs(data[i,0]-data[i,1])
        ihg = int(diff/width)
        hg[ihg] += 1
    print("# histogram of "+fname
          +", width={0:10.5f}".format(width))
    ld = len(data)
    for i in range(len(hg)):
        print("{0:10.5f} {1:6d} {2:8.2f}".format(i*width+0.5*width,
                                                 hg[i],
                                                 float(hg[i])/ld))

if __name__ == "__main__":

    args = docopt(__doc__)
    # dfile = args['DATAFILE']
    width = float(args['--width'])
    hist = args['--histogram']

    make_fit_graph(key='erg')
    make_fit_graph(key='frc')
    make_fit_graph(key='strs')

    make_iter_graph()
    
    # if hist:
    #     histogram(dfile,width)

