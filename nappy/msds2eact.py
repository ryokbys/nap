#!/usr/bin/env python
"""
Compute activation energy from temperature vs diffusion coefficient.
If PATHS are directories,
PATHS should be like '300K' and the digits before 'K' is used as temperature.
else if PATHS are files,
PATHS should be like 'out.msd.300' and the digit at the end is used as temperature.

NOTICE that the errors in Ds are not exactly the statistical errors, but just regression errors,
which should not be confused and not used as errors in D vs 1/T plot.

Usage:
  msds2eact.py [options] PATHS [PATHS...]

Options:
  -h, --help  Show this message and exit.
  -o OUT      Output file name of D vs T. [default: out.D-T]
  --dim DIM   Spatial dimension of diffusion. [default: 3]
  --offset OFFSET
              Offset of the data. [default: 0]
  --plot
              Plot E_act vs 1/T graph (optional) or write a gnuplot script. [default: False]
  --specorder SPECORDER
              Species order used in the MD simulation, separated by comma. 
              e.g.) Li,Zr,P,O, [default: None]
  --spc SPC
              Migrating species. [default: None]
"""
import os,sys
from docopt import docopt
import numpy as np
from scipy import stats
from functools import cmp_to_key
from nappy.common import get_key

__author__ = "RYO KOBAYASHI"
__version__ = "250509"

_kB = 8.6173303e-5

def make_gnuplot_file(outDT,Eact,D0):
    txt = """set format y '10^{{%L}}'
set xl '1000/T (1/K)'
set yl 'D [cm^2/sec]'
set log y
f(x) = {0:.3e} *exp(-{1:.3f} /8.617e-5 *(x/1000))
plot '{2:s}' us (1000.0/$1):2:3 w yerr lc 'blue' pt 7 t 'data', f(x) t 'fitted' lc 'blue'
""".format(D0,Eact,outDT)
    with open('plot_D-T.gp','w') as f:
        f.write(txt)
    print(' Wrote plot_D-T.gp')
    return None

def cmp(a,b):
    if a == b: return 0
    return -1 if a < b else 1

def cmpstr_dir(a,b):
    return cmp(int(a.replace('K','')),int(b.replace('K','')))

def cmpstr_file(a,b):
    return cmp(int(a.replace('out.msd.','')),int(b.replace('out.msd.','')))

def paths2files(paths=[]):
    assert len(paths) > 0, 'len(paths) == 0'
    pathtype = None
    if os.path.isfile(paths[0]):
        pathtype = 'file'
    elif os.path.isdir(paths[0]):
        pathtype = 'dir'
    assert pathtype != None, 'Unknow pathtype...'

    files = []
    temps = []
    if pathtype == 'dir':
        paths.sort(key=cmp_to_key(cmpstr_dir),reverse=True)
        for d in paths:
            T = d.replace('K','')
            files.append(d+'/out.msd')
            temps.append(float(T))
    elif pathtype == 'file':
        paths.sort(key=cmp_to_key(cmpstr_file),reverse=True)
        for f in paths:
            T = f.replace('out.msd.','')
            files.append(f)
            temps.append(float(T))
    return files,temps

def msds2Ds(files,temps,dim=3,offset=0,specorder=[],spc=None):
    from nappy.msd2diff import read_out_msd, msd2D
    Ts = []
    Ds = []
    Dstds = []
    fac = 1.0e-16 /1.0e-15  # A^2/fs to cm^2/s
    for i,f in enumerate(files):
        T = temps[i]
        # ts,msds = read_out_msd(d+'/out.msd',offset=offset,column=specorder.index(spc)+1)
        ts,msds = read_out_msd(f,offset=offset,column=specorder.index(spc)+2)
        D,b,Dstd = msd2D(ts,msds,fac,dim=dim)
        if D < 0.0:
            continue
        print(' T,D = {0:5d}K, {1:0.4e} +/- {2:0.3e} [cm^2/s]'.format(int(T),D,Dstd))
        Ts.append(float(T))
        Ds.append(D)
        Dstds.append(Dstd)
    return np.array(Ts),np.array(Ds),np.array(Dstds)
    
    
if __name__ == "__main__":

    args = docopt(__doc__)
    paths = args['PATHS']
    outfname = args['-o']
    offset = int(args['--offset'])
    plot = args['--plot']
    dim = int(args['--dim'])
    specorder = args['--specorder'].split(',')
    spc = args['--spc']

    if 'None' in specorder or len(specorder) < 1:
        raise ValueError('SPECORDER must be set.')
    if spc == 'None':
        raise ValueError('SPC must be set.')
    if spc not in specorder:
        raise ValueError('SPC must be in SPECORDER')

    files, temps = paths2files(paths)
    Ts,Ds,Dstds = msds2Ds(files,temps,dim=dim,offset=offset,specorder=specorder,spc=spc)

    with open(outfname,'w') as f:
        f.write('# T [K],       D [cm^2/sec],    sgm(D) [cm^2/sec]\n')
        for T,D,Dstd in zip(Ts,Ds,Dstds):
            f.write('    {0:8.2f}  {1:12.4e}  {2:12.4e}\n'.format(T,D,Dstd))
    print(' Wrote {0:s}'.format(outfname))

    Tinvs = []
    logDs = []
    for i,D in enumerate(Ds):
        if D <= 0.0:
            continue
        Tinvs.append( 1.0/Ts[i] )
        logDs.append( np.log(D) )
    Tinvs = np.array(Tinvs)
    logDs = np.array(logDs)
    # Tinvs = np.array([ 1.0/T for T in Ts ])
    # logDs = np.array([ np.log(D) for D in Ds ])
    a,b,r,p,stderr = stats.linregress(Tinvs,logDs)
    Eact = a *(-_kB) /np.log(np.exp(1.0))
    Eaerr = stderr *_kB /np.log(np.exp(1.0))
    D0 = np.exp(b)
    print(' Ea = {0:.3e} +/- {1:.3e} [eV]'.format(Eact,Eaerr))
    print(' D0 = {0:.4e} [cm^2/s]'.format(D0))

    make_gnuplot_file(outfname,Eact,D0)
    
    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
        plt.xlabel('1000/T (1/K)')
        plt.ylabel('D (cm$^2$/s)')
        #plt.plot(Tinvs*1000,log10Ds,'bo',label='data')
        log10Dstds = [ np.log10(Dstd) for Dstd in Dstds ]
        plt.errorbar(Tinvs*1000,Ds,yerr=Dstds, capsize=5, fmt='o',
                     markersize=8,ecolor='k', markeredgecolor='k',
                     color='w',label='data')
        fvals = np.array([ D0*np.exp(-Eact/_kB/T) for T in Ts])
        plt.plot(Tinvs*1000,fvals,color='black',linestyle='dashed',label='fitted')
        plt.yscale('log')
        plt.savefig("graph_msds2eact.png", format='png',
                    dpi=300, bbox_inches='tight')
        print(' Wrote graph_msds2eact.png')
