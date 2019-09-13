#!/usr/bin/env python
"""
Fit parameters of classical interatomic potentials using metaheuristic approaches.
Because of computational efficiency, optimization of a lot of parameters using 
gradient-based approaches is out of focus in this script, and that should be 
performed in a Fortran program.

Usage:
  fp.py [options]

Options:
  -h, --help  Show this message and exit.
  --nproc NPROC
              Number of processes to be used. [default: 1]
  --pmddir-prefix PREFIX
              Prefix for pmd directory. [default: pmddir_]
  --pmd-script SCRIPT
              Name of script that performs pmd and post-processing. [default: pmd_rdf_adf_vol.sh]
"""
from __future__ import print_function

import os,sys
import shutil
from docopt import docopt
import numpy as np
from numpy import sin,cos,exp,sqrt
import subprocess
import time
import json

from nappy.fitpot.fp2prms import fp2BVSx, fp2BVS, fp2Morse, read_params_Coulomb

__author__ = "RYO KOBAYASHI"
__version__ = "rev190904"


def read_in_fitpot(fname='in.fitpot'):
    #...initialize
    infp = {}
    infp['rdf_match'] = True
    infp['adf_match'] = True
    infp['vol_match'] = True
    infp['lat_match'] = False
    infp['fval_upper_limit'] = 1.0e+10
    infp['print_level'] = 1
    infp['weights'] = {'rdf':1.0, 'adf':1.0, 'vol':1.0, 'lat':1.0}
    
    mode = None
    specorder = None
    interact = []
    rdf_pairs = []
    infp['interact'] = interact
    infp['rdf_pairs'] = rdf_pairs
    
    with open(fname,'r') as f:
        lines = f.readlines()

    for line in lines:
        if line[0] in ('!','#'):
            mode = None
            continue
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if data[0] == 'num_iteration':
            maxiter = int(data[1])
            infp['num_iteration'] = maxiter
            mode = None
        elif data[0] == 'print_level':
            print_level = int(data[1])
            infp['print_level'] = print_level
            mode = None
        elif data[0] == 'fitting_method':
            fit_method = data[1]
            infp['fitting_method'] = fit_method
            mode = None
        elif data[0] == 'sample_directory':
            sampledir = data[1]
            if '"' in sampledir:
                sampledir = sampledir.replace('"','')
            elif "'" in sampledir:
                sampledir = sampledir.replace("'",'')
            infp['sample_directory'] = sampledir
            mode = None
        elif data[0] == 'param_file':
            prmfile = data[1]
            infp['param_file'] = prmfile
            mode = None
        elif data[0] == 'potential':
            potential = data[1]
            infp['potential'] = potential
            mode = None
        elif data[0] == 'fval_upper_limit':
            fup_limit = float(data[1])
            infp['fval_upper_limit'] = fup_limit
            mode = None
        elif data[0] == 'specorder':
            specorder = data[1:]
            infp['specorder'] = specorder
            mode = None
        elif data[0] == 'interactions':
            mode = 'interactions'
            nint = int(data[1])
        elif data[0] == 'rdf_pairs':
            mode = 'rdf_pairs'
            nint = int(data[1])
        elif data[0] == 'sample_error':
            mode = 'sample_error'
        elif data[0] == 'rdf_match':
            rdf_match = True if data[1] in ('true', 'True', 'T', 'TRUE') else False
            infp['rdf_match'] = rdf_match
            if rdf_match and len(data) > 2:
                weight = float(data[2])
                infp['weights']['rdf'] = weight
            mode = None
        elif data[0] == 'adf_match':
            adf_match = True if data[1] in ('true', 'True', 'T', 'TRUE') else False
            infp['adf_match'] = adf_match
            if adf_match and len(data) > 2:
                weight = float(data[2])
                infp['weights']['adf'] = weight
            mode = None
        elif data[0] == 'vol_match':
            vol_match = True if data[1] in ('true', 'True', 'T', 'TRUE') else False
            infp['vol_match'] = vol_match
            if vol_match and len(data) > 2:
                weight = float(data[2])
                infp['weights']['vol'] = weight
            mode = None
        elif data[0] == 'lat_match':
            lat_match = True if data[1] in ('true', 'True', 'T', 'TRUE') else False
            infp['lat_match'] = lat_match
            if lat_match and len(data) > 2:
                weight = float(data[2])
                infp['weights']['lat'] = weight
            mode = None
        elif data[0] == 'de_num_individuals':
            nind = int(data[1])
            infp['de_num_individuals'] = nind
            mode = None
        elif data[0] == 'de_crossover_rate':
            cr = float(data[1])
            infp['de_crossover_rate'] = cr
            mode = None
        elif data[0] == 'de_fraction':
            frac = float(data[1])
            infp['de_fraction'] = frac
            mode = None
        elif data[0] == 'de_temperature':
            temp = float(data[1])
            infp['de_temperature'] = temp
            mode = None
        else:
            if mode == 'interactions' and len(data) in (2,3):
                interact.append(tuple(data))
                infp['interactions'] = interact
            elif mode == 'rdf_pairs' and len(data) == 2:
                rdf_pairs.append(tuple(data))
                infp['rdf_pairs'] = rdf_pairs
            else:
                mode = None
                pass
    
    return infp
    
def write_vars_fitpot(vs,vrs,fname='in.vars.fitpot',**kwargs):
    rc2 = kwargs['rc2']
    rc3 = kwargs['rc3']
    nv = len(vs)
    with open(fname,'w') as f:
        f.write(' {0:5d}  {1:7.3f} {2:7.3f}\n'.format(nv,rc2,rc3))
        for i in range(len(vs)):
            f.write(' {0:10.4f}  {1:10.4f}  {2:10.4f}\n'.format(vs[i],*vrs[i]))
    return None

def read_vars_fitpot(fname='in.vars.fitpot'):
    with open(fname,'r') as f:
        lines = f.readlines()
    iv = 0
    nv = -1
    rc2 = 5.0
    rc3 = 3.0
    vs = []
    vrs = []
    for line in lines:
        if line[0] in ('!','#'):
            continue
        data = line.split()
        if len(data) == 0:
            continue
        if nv < 0:
            nv = int(data[0])
            rc2 = float(data[1])
            rc3 = float(data[2])
            continue
        vs.append(float(data[0]))
        vrs.append([ float(data[1]), float(data[2])])
    vs = np.array(vs)
    vrs = np.array(vrs)
    return rc2,rc3,vs,vrs

def read_rdf(fname,specorder,pairs=[]):
    """
    Read reference/FF RDF from data.(ref/pmd).rdf.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    #...Make column indices for pairs
    icols = {}
    if len(pairs) == 0:
        raise ValueError("No pair !")
    for ip,p in enumerate(pairs):
        inc = 2
        for i in range(len(specorder)):
            si = specorder[i]
            for j in range(i,len(specorder)):
                sj = specorder[j]
                if set(p) == set([si,sj]):
                    icols[p] = inc
                inc += 1
    
    rs = []
    rdfs = {}
    for p in pairs:
        rdfs[p] = []
    
    for line in lines:
        if line[0] in ('!','#'):
            continue
        data = line.split()
        if len(data) == 0:
            continue
        rs.append(float(data[0]))
        for ip,p in enumerate(pairs):
            icol = icols[p]
            rdfs[p].append(float(data[icol]))
    return rs,rdfs

def read_adf(fname,specorder,triplets=[]):
    """
    Read reference/FF ADF from data.(ref/pmd).adf.X-Y-Z where X,Y,Z are element name.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    ths = []
    adfs = {}
    for t in triplets:
        adfs[t] = []
        
    for line in lines:
        if line[0] in ('!','#'):
            continue
        data = line.split()
        if len(data) == 0:
            continue
        ths.append(float(data[0]))
        for it,t in enumerate(triplets):
            adfs[t].append(float(data[it+1]))
    return ths, adfs

def read_vol(fname):
    """
    REead reference/FF volume from data.(ref/pmd).vol.
    """
    with open(fname,'r') as f:
        vol = float(f.readline())
    return vol

def read_lat(fname):
    """
    Read reference/FF lattice parameters from data.(ref/pmd).lat.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    a,b,c,alp,bet,gmm = [ float(x) for x in lines[0].split() ]
    return a,b,c,alp,bet,gmm

def get_data(basedir,prefix='ref',**kwargs):
    """
    Get rdf, adf, vol from a given basedir.
    The prefix should be either ref or pmd.
    """

    specorder = kwargs['specorder']
    rdf_pairs = kwargs['rdf_pairs']
    triplets = kwargs['triplets']

    rs = []
    rdfs = []
    ths = []
    adfs = []
    vol = 0.0
    a=b=c=alp=bet=gmm= 0.0
    if kwargs['rdf_match']:
        rs,rdfs = read_rdf(basedir+'/data.{0:s}.rdf'.format(prefix),specorder,rdf_pairs)
    if kwargs['adf_match']:
        ths,adfs = read_adf(basedir+'/data.{0:s}.adf'.format(prefix),specorder,triplets)
    if kwargs['vol_match']:
        vol = read_vol(basedir+'/data.{0:s}.vol'.format(prefix))
    if kwargs['lat_match']:
        a,b,c,alp,bet,gmm = read_lat(basedir+'/data.{0:s}.lat'.format(prefix))
    data = {'rs':rs,
            'rdfs':rdfs,
            'ths':ths,
            'adfs':adfs,
            'vol':vol,
            'lat':(a,b,c,alp,bet,gmm)}
    return data

def loss_func(pmddata,**kwargs):
    """
    Compute loss function value from reference and pmd data.
    """
    refdata = kwargs['refdata']
    wgts = kwargs['weights']
    rdf_pairs = kwargs['rdf_pairs']
    triplets = kwargs['triplets']

    #...RDF
    Lr = 0.0
    if kwargs['rdf_match'] and len(rdf_pairs) != 0:
        rs = refdata['rs']
        rdfs_ref = refdata['rdfs']
        rdfs_pmd = pmddata['rdfs']
        for p in rdf_pairs:
            # print(p)
            rdf_ref = rdfs_ref[p]
            rdf_pmd = rdfs_pmd[p]
            diff2sum = 0.0
            z = 0.0
            for i,r in enumerate(rs):
                ref = rdf_ref[i]
                pmd = rdf_pmd[i]
                diff = pmd -ref
                diff2sum += diff*diff
                z += ref*ref
                # print('i,r,ref,pmd,z,diff2sum=',i,r,ref,pmd,z,diff2sum)
            Lr += diff2sum/z #/len(rs)
        Lr /= len(rdf_pairs)

    #...ADF
    Lth = 0.0
    if kwargs['adf_match'] and len(triplets) != 0:
        ths = refdata['ths']
        adfs_ref = refdata['adfs']
        adfs_pmd = pmddata['adfs']
        for t in triplets:
            adf_ref = adfs_ref[t]
            adf_pmd = adfs_pmd[t]
            diff2sum = 0.0
            z = 0.0
            for i,th in enumerate(ths):
                ref = adf_ref[i]
                pmd = adf_pmd[i]
                diff = pmd -ref
                diff2sum += diff*diff
                z += ref*ref
            Lth += diff2sum /z #/len(ths)
        Lth /= len(triplets)

    #...volume
    Lvol = 0.0
    if kwargs['vol_match']:
        vol_ref = refdata['vol']
        vol_pmd = pmddata['vol']
        diff = vol_pmd -vol_ref
        Lvol = diff*diff /(vol_ref*vol_ref)

    #...lattice parameters
    Llat = 0.0
    if kwargs['lat_match']:
        a0,b0,c0,alp0,bet0,gmm0 = refdata['lat']
        # print('refdata=',refdata['lat'])
        a,b,c,alp,bet,gmm = pmddata['lat']
        # print('pmddata=',pmddata['lat'])
        #...Need to take into account the definition difference bet/ dump and vasp
        # a0l,b0l,c0l,alp0l,bet0l,gmm0l = lat_vasp2dump(a0,b0,c0,alp0,bet0,gmm0)
        diff = ((a0-a)/a0)**2 +((b0-b)/b0)**2 +((c0-c)/c0)**2 \
                +((alp0-alp)/alp0)**2 +((bet0-bet)/bet0)**2 +((gmm0-gmm)/gmm0)**2
        # diff /= 6
        Llat = diff

    lrw = Lr*wgts['rdf']
    lthw = Lth*wgts['adf']
    lvolw = Lvol*wgts['vol']
    llatw = Llat*wgts['lat']
    L = lrw +lthw +lvolw +llatw

    if kwargs['print_level'] > 0:
        print(' iid,Lr,Lth,Lvol,Llat,L= {0:8d}'.format(kwargs['iid'])
              +'{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f}'.format(lrw,lthw,lvolw,llatw,L))
    return L

def func_wrapper(variables, vranges, **kwargs):
    """
    Wrapper function for the above loss_func().
    This converts variables to be optimized to parameters for pmd,
    perform pmd, then get rdf, adf, and vol from the pmd result.
    Then give them to the above loss_func().
    """
    # infp = kwargs['infp']
    # pairs = kwargs['pairs']
    # triplets = kwargs['triplets']
    refdata = kwargs['refdata']
    wgts = kwargs['weights']
    specorder = kwargs['specorder']
    refdata = kwargs['refdata']
    pmdscript = kwargs['pmd-script']
    pmddir = kwargs['pmddir-prefix'] +'{0:03d}'.format(kwargs['index'])
    print_level = kwargs['print_level']
    # print('refdata=',refdata)
    # print('pairs=',pairs)
    # print('triplets=',triplets)

    #...Create in.params.XXX files in each pmddir
    varsfp = {}
    varsfp['rc2'] = kwargs['rc2']
    varsfp['rc3'] = kwargs['rc3']
    varsfp['variables'] = variables
    varsfp['vranges'] = vranges
    cwd = os.getcwd()
    if not os.path.exists(pmddir):
        os.mkdir(pmddir)
        shutil.copy(pmdscript,pmddir+'/')
    os.chdir(pmddir)
    if not 'vids' in kwargs.keys():
        print(kwargs.keys())

    if kwargs['potential'] == 'BVSx':
        fp2BVSx(varsfp, **kwargs)
    elif kwargs['potential'] == 'BVS':
        fp2BVS(varsfp, **kwargs)
    elif kwargs['Morse'] == 'Morse':
        fp2Morse(varsfp, **kwargs)

    #...Compute pmd
    L_up_lim = kwargs['fval_upper_limit']
    if print_level > 1:
        print('Running pmd and post-processing at '+pmddir, flush=True)
    try:
        cmd = "./{0:s} > log.iid_{1:d}".format(pmdscript,kwargs['iid'])
        # subprocess.run(cmd.split(),check=True)
        subprocess.run(cmd,shell=True,check=True)
        os.chdir(cwd)
        # print('Going to get_data from ',pmddir)
        pmddata = get_data(pmddir,prefix='pmd',**kwargs)
        L = min( loss_func(pmddata,**kwargs), L_up_lim )
    except:
        if print_level > 1:
            print('  Since pmd or post-process failed at {0:s}, '.format(pmddir)
                  +'the upper limit value is applied to its loss function.',
                  flush=True)
        os.chdir(cwd)
        L = L_up_lim
        
    return L

def get_triplets(interact):
    triplets = []
    for it in interact:
        if len(it) == 3:
            triplets.append(it)
    return triplets

def get_pairs(interact):
    pairs = []
    for it in interact:
        if len(it) == 2:
            pairs.append(it)
    return pairs

def latprms2hmat(a,b,c,alp,bet,gmm):
    """
    Convert lattice parameters to hmat.
    See https://arxiv.org/pdf/1506.01455.pdf
    """
    alpr = np.radians(alp)
    betr = np.radians(bet)
    gmmr = np.radians(gmm)
    # val = (cos(alpr) * cos(betr) - cos(gmmr))\
    #       / (sin(alpr) * sin(betr))
    # val = max(abs(val),1.0)
    # gmmstar = np.arccos(val)
    
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    a1[:] = [a, 0.0, 0.0]
    a2[:] = [b*cos(gmmr), b*sin(gmmr), 0.0]
    # a3[:] = [c*np.cos(betr),
    #          -c*np.sin(betr)*np.cos(gmmstar),
    #          c*np.sin(betr)*np.sin(gmmstar)]
    a3[:] = [c*cos(betr),
             c*(cos(alpr) -cos(betr)*cos(gmmr))/sin(gmmr),
             c*sqrt(sin(gmmr)**2 -cos(alpr)**2 -cos(betr)**2
                    +2.0*cos(alpr)*cos(betr)*cos(gmmr))/sin(gmmr)]
    hmat = np.zeros((3,3))
    hmat[:,0] = a1[:]
    hmat[:,1] = a2[:]
    hmat[:,2] = a3[:]
    return hmat

def lat_vasp2dump(a,b,c,alpha,beta,gamma):
    from nappy.napsys import to_lammps

    try:
        hmat = latprms2hmat(a,b,c,alpha,beta,gamma)
        # print('hmat=',hmat)
        xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,_ = to_lammps(hmat,[])
        # print('after to_lammps=',xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
    except:
        raise

    a1 = np.array([xhi-xlo,     0.0,     0.0])
    b1 = np.array([     xy, yhi-ylo,     0.0])
    c1 = np.array([     xz,      yz, zhi-zlo])
    al = np.linalg.norm(a1)
    bl = np.linalg.norm(b1)
    cl = np.linalg.norm(c1)
    alpl = np.arccos(np.dot(b1,c1)/bl/cl) /np.pi *180.0
    betl = np.arccos(np.dot(a1,c1)/al/cl) /np.pi *180.0
    gmml = np.arccos(np.dot(a1,b1)/al/bl) /np.pi *180.0

    return al,bl,cl,alpl,betl,gmml
    
def main(args):
    from fitpot.de import DE

    start = time.time()

    nproc = int(args['--nproc'])

    infp = read_in_fitpot('in.fitpot')
    pairs = get_pairs(infp['interactions'])
    rdf_pairs = infp['rdf_pairs']
    if len(rdf_pairs) == 0:  # if no rdf_pairs are specied, all the pairs are selected
        specorder = infp['specorder']
        rdf_pairs = []
        for i,si in enumerate(specorder):
            for j in range(i,len(specorder)):
                sj = specorder[j]
                rdf_pairs.append((si,sj))
    
    # print('pairs    =',pairs)
    # print('rdf_pairs=',rdf_pairs)
    triplets = get_triplets(infp['interactions'])
    rc2,rc3,vs,vrs = read_vars_fitpot(infp['param_file'])


    kwargs = infp
    # kwargs['infp'] = infp
    kwargs['rc2'] = rc2
    kwargs['rc3'] = rc3
    kwargs['pairs'] = pairs
    kwargs['rdf_pairs'] = rdf_pairs
    kwargs['triplets'] = triplets
    kwargs['pmddir-prefix'] = args['--pmddir-prefix']
    kwargs['pmd-script'] = args['--pmd-script']
    kwargs['start'] = start
    
    smpldir = infp['sample_directory']
    refdata = get_data(smpldir,prefix='ref',**kwargs)
    a0,b0,c0,alp0,bet0,gmm0 = refdata['lat']
    #...Need to take into account the definition difference bet/ dump and vasp
    a,b,c,alp,bet,gmm = lat_vasp2dump(a0,b0,c0,alp0,bet0,gmm0)
    # print('Reference lattice parameters:',a,b,c,alp,bet,gmm)
    refdata['lat'] = (a,b,c,alp,bet,gmm)
    kwargs['refdata'] = refdata

    fbvs, rads, vids, npqs = read_params_Coulomb('in.params.Coulomb')
    kwargs['fbvs'] = fbvs
    kwargs['rads'] = rads
    kwargs['vids'] = vids
    kwargs['npqs'] = npqs


    N = infp['de_num_individuals']
    F = infp['de_fraction']
    T = infp['de_temperature']
    CR = infp['de_crossover_rate']
    de = DE(N,F,CR,T, vs,vrs, func_wrapper, write_vars_fitpot, **kwargs)

    maxiter = kwargs['num_iteration']
    de.run(maxiter)

    print('elapsed time = {0:f} sec.'.format(time.time()-start))
    
    return None

if __name__ == "__main__":

    args = docopt(__doc__)

    main(args)