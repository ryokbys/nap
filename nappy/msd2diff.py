#!/usr/bin/env python
"""
Compute diffusion coefficient from MSD data.
Time interval, DT, is obtained from in.pmd in the same directory.

Usage:
  msd2diff.py [options] MSD_FILE

Options:
  -h, --help  Show this message and exit.
  --offset OFFSET
              Offset of the given data. [default: 0]
  --plot      Plot a fitted graph. [default: False]
  --out4fp    Output data file for fp.py in any-target mode.
  --out4fp-name OUTNAME
              File name for --out4fp. [default: data.pmd.D] 
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "191212"


def read_out_msd(fname='out.msd',offset=0,specorder=[],spc=None):
    if specorder == [] or spc not in specorder:
        index = 1
    else:
        index = specorder.index(spc) +1
        
    with open(fname,'r') as f:
        lines = f.readlines()
    try:
        dname = os.path.dirname(fname)
        if len(dname) == 0:
            dname = '.'
        dt = dt_from_inpmd(fname='/'.join([dname,'in.pmd']))
    except Exception as e:
        raise RuntimeError('Failed to read in.pmd.')
    ts = []
    msds = []
    n0 = 0
    msd0 = 0.0
    for il,line in enumerate(lines):
        if line[0] == '#':
            continue
        data = line.split()
        if il < offset:
            n0 = int(data[0])
            msd0 = float(data[index])
            continue
        n = int(data[0])
        msd = float(data[index])
        ts.append((n-n0)*dt)
        msds.append(msd-msd0)
    return np.array(ts),np.array(msds)

def dt_from_inpmd(fname='in.pmd'):
    with open(fname,'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'time_interval' in line:
            time_interval = abs(float(line.split()[1]))
        elif 'num_iteration' in line:
            num_iteration = int(line.split()[1])
        elif 'num_out_pos' in line or 'num_out_pmd' in line:
            num_out_pos = int(line.split()[1])
    
    return time_interval*num_iteration/num_out_pos

def msd2D(ts,msds,fac,dim=3):
    """
    Compute diffusion coefficient from time [fs] vs MSD [Ang^2] data 
    by solving least square problem using numpy.
    Return diffusion coefficient multiplied by FAC.
    """
    A= np.array([ts, np.ones(len(ts))])
    A = A.T
    xvar = np.var(A[:,0])
    p,res,_,_ = np.linalg.lstsq(A,msds,rcond=None)
    a = p[0]
    b = p[1]
    # fac = 1.0e-16 /1.e-15
    a = a *fac /(2.0*dim)
    b = b *fac
    # print(res[0],xvar,np.mean(A[:,0]),len(ts))
    std = np.sqrt(res[0]/len(ts)/xvar) *fac /(2.0*dim)
    return a,b,std


if __name__ == "__main__":

    args = docopt(__doc__)
    fname = args['MSD_FILE']
    offset = int(args['--offset'])
    plot = args['--plot']
    out4fp = args['--out4fp']
    out4fpname = args['--out4fp-name']

    ts,msds = read_out_msd(fname,offset)
    #...Assuming input MSD unit in A^2/fs and output in cm^2/s
    fac = 1.0e-16 /1.0e-15
    #...Least square
    a,b,std = msd2D(ts,msds,fac)
    print(' Diffusion coefficient = {0:12.4e}'.format(a)+
          ' +/- {0:12.4e} [cm^2/s]'.format(std))

    if out4fp:
        from datetime import datetime
        with open(out4fpname,'w') as f:
            cmd = ' '.join(s for s in sys.argv)
            f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('#  {0:s}\n'.format(cmd))
            f.write('#\n')
            f.write('    1    1.0\n')
            f.write('   {0:12.4}\n'.format(a))

    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
        #...Original time unit == fs
        unit = 'fs'
        tfac = 1.0
        if ts[-1] > 1.0e+5: #...if max t > 100ps, time unit in ps
            unit = 'ps'
            tfac = 1.0e-3
        plt.xlabel('Time ({0:s})'.format(unit))
        plt.ylabel('MSD (A^2/{0:s})'.format(unit))
        fvals = np.array([ (t*a+b)/fac for t in ts ])
        plt.plot(ts*tfac,msds/tfac,'b-',label='MSD data')
        plt.plot(ts*tfac,fvals/tfac,'r-',label='Fitted curve')
        plt.savefig("graph_msd2D.png", format='png',
                    dpi=300, bbox_inches='tight')
        print(' Wrote graph_msd2D.png')
        
