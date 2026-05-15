#!/usr/bin/env python
"""
Compute diffusion coefficient from MSD data.
Optionally compute ionic conductivity via the Nernst-Einstein relation.

Usage:
  msd2diff.py [options] MSD_FILE [MSD_FILE...]

Options:
  -h, --help  Show this message and exit.
  --dim DIM   Dimension of diffusion. [default: 3]
  --offset OFFSET
              Offset of the given data. [default: 0]
  --main-SID MAINSID
              Species ID whose MSD is to be extracted. [default: 1]
  --subtract-SID SUBSID
              Species ID whose MSD is to be subtracted from that of main-SID.
              If this is less than 1, nothing is subtracted. [default: 0]
  --plot      Plot a fitted graph. [default: False]
  --out4fp    Output data file for fp.py in any-target mode.
  --out4fp-name OUTNAME
              File name for --out4fp. [default: data.pmd.D]
  --structure STRUCT
              Structure file (pmdini, POSCAR, etc.) from which cell volume and
              number of mobile ions (matching --main-SID) are read. [default: None]
  --natoms NATOMS
              Number of mobile ions. Overrides the value from --structure. [default: None]
  --volume VOLUME
              Cell volume in Angstrom^3. Overrides the value from --structure. [default: None]
  --temperature TEMP
              Temperature in K. Required for Nernst-Einstein conductivity. [default: None]
  --charge CHARGE
              Charge number (valence) of mobile ions. [default: 1]
"""
import os
import sys
from docopt import docopt
import numpy as np
from nappy.util import parse_option, gen_header

__author__ = "RYO KOBAYASHI"
__version__ = "260515"

def read_out_msd(fname='out.msd',offset=0,column=2):

    with open(fname,'r') as f:
        lines = f.readlines()
    ts = []
    msds = []
    n0 = 0
    t0 = 0.0
    msd0 = 0.0
    dt = -1.0
    for il,line in enumerate(lines):
        if line[0] == '#':
            opt = parse_option(line)
            if opt != None:
                if 'dt' in opt.keys():
                    dt = float(opt['dt'])
                continue
            else:
                continue
        data = line.split()
        if il < offset:
            n0 = int(data[0])
            t0 = float(data[1])
            msd0 = float(data[column])
            continue
        n = int(data[0])
        t = float(data[1])
        msd = float(data[column])
        ts.append(t - t0)
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
    # fac = 1.0e-5 /1.e-4
    a = a *fac /(2.0*dim)
    b = b *fac
    # print(res[0],xvar,np.mean(A[:,0]),len(ts))
    std = np.sqrt(res[0]/len(ts)/xvar) *fac /(2.0*dim)
    return a,b,std


def nernst_einstein(D, natoms, volume_ang3, temperature, charge=1):
    """
    Compute ionic conductivity from diffusion coefficient via Nernst-Einstein.

    Parameters
    ----------
    D : float
        Diffusion coefficient [cm^2/s]
    natoms : int
        Number of mobile ions in the cell
    volume_ang3 : float
        Cell volume [Angstrom^3]
    temperature : float
        Temperature [K]
    charge : int
        Charge number (valence) of mobile ions

    Returns
    -------
    sigma : float
        Ionic conductivity [S/cm]
    """
    e  = 1.60217663e-19  # C
    kB = 1.380649e-23    # J/K
    V_cm3 = volume_ang3 * 1.0e-24  # Ang^3 -> cm^3
    n = natoms / V_cm3              # number density [1/cm^3]
    sigma = n * (charge * e)**2 * D / (kB * temperature)
    return sigma


def main():

    #args = docopt(__doc__)
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    fnames = args['MSD_FILE']
    offset = int(args['--offset'])
    dim = int(args['--dim'])
    sidmain = int(args['--main-SID'])
    sidsub = int(args['--subtract-SID'])
    plot = args['--plot']
    out4fp = args['--out4fp']
    out4fpname = args['--out4fp-name']
    struct_fname = args['--structure']
    temperature = float(args['--temperature']) if args['--temperature'] != 'None' else None
    charge = int(args['--charge'])

    # Resolve natoms and volume from structure file or direct options
    natoms = None
    volume = None
    if struct_fname not in (None, 'None'):
        import nappy.io
        nsys = nappy.io.read(fname=struct_fname)
        if type(nsys) == list:
            nsys = nsys[0]
        volume = nsys.get_volume()  # Ang^3
        natoms = int((nsys.atoms.sid == sidmain).sum())
        print(f' Structure: {struct_fname}')
        print(f'   Cell volume = {volume:.4f} Ang^3')
        print(f'   Num of mobile ions (SID={sidmain}) = {natoms}')
    if args['--natoms'] != 'None':
        natoms = int(args['--natoms'])
    if args['--volume'] != 'None':
        volume = float(args['--volume'])

    do_conductivity = (temperature is not None and natoms is not None and volume is not None)

    Ds = []
    Bs = []
    MSDs = []
    Ts = []
    for fname in fnames:
        ts,msdmain = read_out_msd(fname,offset,column=sidmain+1)
        if sidsub > 0:
            tmp, msdsub = read_out_msd(fname,offset,column=sidsub+1)
            msdmain = msdmain -msdsub
        #...Assuming input MSD unit in A^2/fs and output in cm^2/s
        fac = 1.0e-16 /1.0e-15
        #...Least square
        D,b,std = msd2D(ts,msdmain,fac,dim=dim)
        print(f' MSD: {fname:s}')
        print(f'   Diffusion coefficient   = {D:0.4e}'+
              f' +/- {std:0.3e} [cm^2/s]')
        if do_conductivity:
            sigma = nernst_einstein(D, natoms, volume, temperature, charge)
            sigma_std = nernst_einstein(std, natoms, volume, temperature, charge)
            print(f'   Ionic conductivity (NE) = {sigma*1e3:.4e}'
                  f' +/- {sigma_std*1e3:.3e} [mS/cm]'
                  f'  (T={temperature:.1f} K, z={charge})')
        Ds.append(D)
        Bs.append(b)
        Ts.append(ts)
        MSDs.append(msdmain)

    if out4fp:
        from datetime import datetime
        with open(out4fpname,'w') as f:
            f.write(gen_header(sys.argv))
            f.write('# datatype: independent\n')
            f.write('#\n')
            ndat = len(Ds)
            f.write(f'   {ndat:d}     1.0\n')
            for i,D in enumerate(Ds):
                f.write(f'  {D:12.4e}  {D:12.4e}\n')

    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(context='talk',style='ticks')
        cmap = plt.get_cmap("tab10")
        #...Original time unit == fs
        unit = 'fs'
        tfac = 1.0
        if ts[-1] > 1.0e+5: #...if max t > 100ps, time unit in ps
            unit = 'ps'
            tfac = 1.0e-3
        plt.xlabel('Time ({0:s})'.format(unit))
        plt.ylabel(r'MSD ($\mathrm{\AA}^2$)')
        for i,fname in enumerate(fnames):
            #...cm^2/s --> A^2/fs
            D = Ds[i] /fac *(2*dim)
            #...cm^2 --> A^2
            b = Bs[i] /fac
            ts = Ts[i]
            MSD = MSDs[i]
            fvals = np.array([ t*D+b for t in ts ])
            c = cmap(i)
            plt.plot(ts*tfac,MSD/tfac,'-',color=c,label=f'MSD ({fname:s})')
            plt.plot(ts*tfac,fvals/tfac,'--',color=c,label=f'Fit ({fname:s})')
        #plt.legend(loc='best')
        plt.savefig("graph_msd2D.png", format='png',
                    dpi=300, bbox_inches='tight')
        print(' Wrote graph_msd2D.png')

    return None


if __name__ == "__main__":

    main()
