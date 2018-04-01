#!/usr/bin/env python
"""
Analyze energies and forces (stresses?) of specified samples.
The sample directory name should start with 'smpl_' and
the name follows before the second '_' is used for coloring.

Usage:
  analyze_samples.py [options] DIRS [DIRS...]

Options:
  -h, --help  Show this message and exit.
  --graph-format FORMAT
              Specify a graph format. [default: png]
  --energy-limit ELIM
              Extract sample names of which the energy has large than ELIM. [default: none]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
import glob

from nappy.napsys import NAPSystem

__author__ = "RYO KOBAYASHI"
__version__ = "180224"

_graph_name = 'graph.erg-vs-vol.'

def read_erg(fname='erg.ref'):
    with open(fname,'r') as f:
        erg = float(f.readline().split()[0])
    return erg
        
def read_frc(fname='frc.ref'):
    with open(fname,'r') as f:
        natm = int(f.readline().split()[0])
        frcs = []
        for i in range(natm):
            frcs.append([ float(d) for d in f.readline().split() ])
    return frcs

def read_strs(fname='strs.ref'):
    with open(fname,'r') as f:
        strs = [ float(d) for d in f.readline().split() ]
    return strs

def read_sample(dirname):
    """
    Read system, energy, forces and stress information of given DIRNAME.
    """
    #...The directory must have erg.ref, frc.ref, strs.ref, and pos files.
    files = ('pos','erg.ref','frc.ref','strs.ref')
    for f in files:
        if not os.path.exists(dirname+'/'+f):
            raise RuntimeError('The file '+f+' does not exist in '+dirname)
    #...Read pos first
    nsys = NAPSystem(fname=dirname+'/pos',ffmt='pmd')
    erg = read_erg(fname=dirname+'/erg.ref')
    frcs = read_frc(fname=dirname+'/frc.ref')
    strs = read_strs(fname=dirname+'/strs.ref')
    return nsys,erg,frcs,strs

def statistics(systems):
    e_ave = 0.0
    e_var = 0.0
    e_min = 1e+30
    e_max = -1e+30
    f_ave = 0.0
    f_var = 0.0
    f_min = 1e+30
    f_max = -1e+30
    nf = 0
    for s in systems:
        nsys = s['nsys']
        natm = len(nsys.atoms)
        erg = s['erg']/natm
        frcs = s['frcs']
        #...energy
        e_ave += erg
        e_var += erg*erg
        e_min = min(e_min,erg)
        e_max = max(e_max,erg)
        #...forces
        for i in range(natm):
            for j in range(3):
                nf += 1
                f_ave += frcs[i][j]
                f_var += frcs[i][j]*frcs[i][j]
                f_min = min(f_min,abs(frcs[i][j]))
                f_max = max(f_max,abs(frcs[i][j]))
    e_ave /= len(systems)
    e_var = e_var/len(systems) -e_ave**2
    f_ave /= nf
    f_var = f_var/nf -f_ave**2

    print('Energy per atom:')
    print('  Average:            {0:8.4f}'.format(e_ave))
    print('  Standard deviation: {0:8.4f}'.format(np.sqrt(e_var)))
    print('  Minimum:            {0:8.4f}'.format(e_min))
    print('  Maximum:            {0:8.4f}'.format(e_max))
    print('  Energy range:       {0:8.4f}'.format(e_max-e_min))

    print('Force component:')
    print('  Average:            {0:8.4f}'.format(f_ave))
    print('  Standard deviation: {0:8.4f}'.format(np.sqrt(f_var)))
    print('  Minimum:            {0:8.4f}'.format(f_min))
    print('  Maximum:            {0:8.4f}'.format(f_max))
    print('  Force range:        {0:8.4f}'.format(f_max-f_min))
    return

def angle2color(angle):
    """
    Convert angle in degree [0:360] (float) to color in RGB.
    """
    from matplotlib import colors
    sat = 1.0
    val = 0.8
    hue = angle/360
    return colors.hsv_to_rgb((hue,sat,val))
    

def uniq(arr):
    uniq_arr = []
    for a in arr:
        if a not in uniq_arr:
            uniq_arr.append(a)
    return uniq_arr
    
def draw_graph(systems,uniq_names,graph_format='png',
               graph_name='graph.png'):
    """
    Draw a graph of erv-vs-vol of given systems.
    """
    try:
        import matplotlib.pyplot as plt
    except:
        raise ImportError('Cannot import module matplotlib.pyplot')
    
    try:
        import seaborn as sns
        sns.set(context='poster',style='darkgrid')
    except:
        pass

    cmap = plt.get_cmap('tab10')

    markersize = 10
    if len(systems) > 2000:
        markersize = 5

    num_name = len(uniq_names)
    dangle = 360.0 /num_name
    for i,name in enumerate(uniq_names):
        ergs = []
        vols = []
        #angle = i*dangle
        #color = angle2color(angle)
        color = cmap(i)
        for s in systems:
            if s['name'] != name:
                continue
            nsys = s['nsys']
            natm = len(nsys.atoms)
            erg = s['erg'] /natm
            ergs.append(erg)
            vols.append(nsys.volume()/natm)
        plt.plot(vols,ergs,'o',color=color,mec='black',mew=0.5,
                 ms=markersize,label=name)
    plt.xlabel('Volume (Ang^3/atom)')
    plt.ylabel('Energy (eV/atom)')
    plt.legend(loc='best')
    plt.savefig(graph_name, format=graph_format,
                dpi=300, bbox_inches='tight')
    return

def get_high_energy_samples(systems,elim=1.0):
    emin = 0.0
    for s in systems:
        nsys = s['nsys']
        natm = len(nsys.atoms)
        erg = s['erg']/natm
        emin = min(emin,erg)
    print('Minimum energy = ',emin)
    dnames = []
    for s in systems:
        nsys = s['nsys']
        natm = len(nsys.atoms)
        erg = s['erg']/natm
        if np.abs(erg-emin) > elim:
            dnames.append(s['dname'])
    print('Num of samples over ELIM = ',len(dnames))
    return dnames

def arrange_dirs(dirs):

    #...If the dirname contains '/' at the end, remove it.
    for i in range(len(dirs)):
        if dirs[i][-1] == '/':
            dirs[i] = dirs[i][:-1]

    newdirs = []
    #...If the dirname does not contain 'smpl_', look for smpl_ dirs in it.
    for d in dirs:
        basename = os.path.basename(d)
        if not 'smpl_' in basename:
            ds_in_dirs = glob.glob(d+'/smpl_*')
            for dd in ds_in_dirs:
                newdirs.append(dd)
        else:
            newdirs.append(d)
    return newdirs
    

if __name__ == "__main__":

    args = docopt(__doc__)


    dirs = args['DIRS']
    graph_format = args['--graph-format']

    dirs = arrange_dirs(dirs)
    print('Number of dirs = {0:d}'.format(len(dirs)))
    
    systems = []
    uniq_names = []
    for d in dirs:
        name = d.split('/')[-1].split('_')[1]
        try:
            nsys,erg,frcs,strs = read_sample(d)
            if name not in uniq_names:
                uniq_names.append(name)
        except:
            continue
        s = {}
        s['nsys'] = nsys
        s['erg'] = erg
        s['frcs'] = frcs
        s['strs'] = strs
        s['name'] = name
        s['dname'] = d
        systems.append(s)

    statistics(systems)

    graph_name = 'graph.erg-vs-vol.'+graph_format
    draw_graph(systems,uniq_names,graph_format=graph_format,
               graph_name=graph_name)

    print('')
    print('- '+graph_name)

    elim = args['--energy-limit']
    if elim != 'none':
        print('')
        elim = float(elim)
        dnames = get_high_energy_samples(systems,elim)
        with open('out.high_energy_samples','w') as f:
            for d in dnames:
                f.write('{0:s}\n'.format(d))
        print('- out.high_energy_samples')
    
