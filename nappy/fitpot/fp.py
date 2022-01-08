#!/usr/bin/env python
"""
Fit parameters of classical interatomic potentials using a metaheuristic algorithm.
Because of computational efficiency, optimization of a lot of parameters using
gradient-based approaches is out of focus in this program, and that should be
performed in a Fortran program.

Usage:
  fp.py [options]

Options:
  -h, --help  Show this message and exit.
  --nproc NPROC
              Number of processes to be used. If it's less than 1, use as many processes as possible. [default: 0]
  --subdir-prefix PREFIX
              Prefix for pmd directory. [default: subdir_]
  --subjob-script SCRIPT
              Name of script that performs MD and post-processing. [default: subjob.sh]
"""
from __future__ import print_function

import os
import sys
import shutil
from docopt import docopt
import numpy as np
from numpy import sin,cos,sqrt
import subprocess
import time
from datetime import datetime

from nappy.fitpot.fp2prms import fp2BVSx, fp2BVS, fp2Morse, read_params_Coulomb, fp2params
from nappy.fitpot.de import DE
from nappy.fitpot.cs import CS

__author__ = "RYO KOBAYASHI"
__version__ = "rev211111"

def read_in_fitpot(fname='in.fitpot'):
    #...initialize
    infp = {}
    infp['rdf_match'] = True
    infp['adf_match'] = True
    infp['vol_match'] = True
    infp['lat_match'] = False
    infp['fval_upper_limit'] = 100.0
    infp['missing_value'] = 1.0
    infp['print_level'] = 1
    infp['weights'] = {'rdf':1.0, 'adf':1.0, 'vol':1.0, 'lat':1.0}
    infp['update_vrange'] = -1
    infp['param_file'] = 'in.vars.fitpot'
    
    mode = None
    specorder = None
    infp['interactions'] = []
    infp['rdf_pairs'] = []
    infp['adf_triplets'] = []
    infp['match'] = []
    infp['param_files'] = []
    
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
        # elif data[0] == 'param_file':
        #     prmfile = data[1]
        #     infp['param_file'] = prmfile
        #     mode = None
        elif data[0] == 'param_files':
            infp[data[0]] = [ name for name in data[1:] ]
            mode = None
        elif data[0] == 'potential':
            potential = data[1]
            infp['potential'] = potential
            mode = None
        elif data[0] == 'fval_upper_limit':
            fup_limit = float(data[1])
            infp['fval_upper_limit'] = fup_limit
            mode = None
        elif data[0] == 'missing_value':
            misval = float(data[1])
            infp['missing_value'] = misval
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
        elif data[0] == 'adf_triplets':
            mode = 'adf_triplets'
            nint = int(data[1])
        elif data[0] == 'sample_error':
            mode = 'sample_error'
        elif data[0] == 'match':
            if len(data) < 2:
                raise RuntimeError('match entry requires at least one keyword.')
            for i in range(1,len(data)):
                infp['match'].append(data[i])
            mode = None
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
        elif data[0] == 'cs_num_individuals':
            nind = int(data[1])
            infp['cs_num_individuals'] = nind
            mode = None
        elif data[0] == 'cs_fraction':
            frac = float(data[1])
            infp['cs_fraction'] = frac
            mode = None
        elif data[0] == 'update_vrange':
            infp['update_vrange'] = int(data[1])
            mode = None
        else:
            if mode == 'interactions' and len(data) in (2,3):
                infp['interactions'].append(tuple(data))
            elif mode == 'rdf_pairs' and len(data) == 2:
                infp['rdf_pairs'].append(tuple(data))
            elif mode == 'adf_triplets' and len(data) == 3:
                infp['adf_triplets'].append(tuple(data))
            else:
                mode = None
                pass
    
    return infp

def write_info(infp,args):
    """
    Write out information on input parameters for fp.
    """

    print(' Inputs')
    print(' ----------')
    print('   num of processes (given by --nproc option)  ',int(args['--nproc']))
    try:
        if len(infp['param_files']) == 0:
            print('   potential       {0:s}'.format(infp['potential']))
        else:
            print('   param_files  ',end='')
            for fname in infp['param_files']:
                print(f'  {fname}', end='')
            print('')
    except:
        raise
    print('   specorder       ',infp['specorder'])
    fmethod = infp['fitting_method']
    print('   fitting_method  {0:s}'.format(fmethod))
    if fmethod in ('de','DE'):
        print('   num_individuals   ',infp['de_num_individuals'])
        print('   fraction          {0:7.4f}'.format(infp['de_fraction']))
        print('   temparature       {0:7.4f}'.format(infp['de_temperature']))
        print('   crossover_rate    {0:7.4f}'.format(infp['de_crossover_rate']))
    elif fmethod in ('cs','CS'):
        print('   num_individuals   ',infp['cs_num_individuals'])
        print('   fraction          {0:7.4f}'.format(infp['cs_fraction']))
    else:
        print('   There is no such fitting method...')
    print('   num_iteration   {0:d}'.format(infp['num_iteration']))
    print('   missing_value   {0:.1f}'.format(infp['missing_value']))
    print(' ----------')
    print()
    return None

def write_vars_fitpot(vs,vrs,fname='in.vars.fitpot',**kwargs):
    rc2 = kwargs['rc2']
    rc3 = kwargs['rc3']
    options = kwargs['options']
    hardlim = kwargs['hardlim']
    nv = len(vs)
    with open(fname,'w') as f:
        if 'hard-limit' in options.keys() and options['hard-limit']:
            f.write('!  hard-limit:  T\n')
            f.write('!\n')
        f.write(' {0:5d}  {1:7.3f} {2:7.3f}\n'.format(nv,rc2,rc3))
        for i in range(len(vs)):
            f.write(' {0:15.7f}  {1:15.7f}  {2:15.7f}'.format(vs[i],*vrs[i]))
            if 'hard-limit' in options.keys() and options['hard-limit']:
                f.write('  {0:10.4f}  {1:10.4f}\n'.format(*hardlim[i]))
            else:
                f.write('\n')
    return None

def parse_option(line):
    if line[0] not in ('#','!'):
        raise ValueError('Line is not a comment line.')
    data = line.split()
    key = None
    value = None
    if len(data) < 3:
        return key,value
    if 'hard-limit:' in data[1] or 'hard_limit:' in data[1]:
        key = 'hard-limit'
        value = False
        if data[2] not in ('No','no', 'NO', 'False', 'F', 'false'):
            value = True
    return key,value

def read_vars_fitpot(fname='in.vars.fitpot'):
    with open(fname,'r') as f:
        lines = f.readlines()
    iv = 0
    nv = -1
    rc2 = 5.0
    rc3 = 3.0
    vs = []
    vrs = []
    vrsh = []
    options = {}
    for line in lines:
        if line[0] in ('!','#'):
            k,v = parse_option(line)
            if k is not None:
                options[k] = v
                print(' option: ',k,v)
            continue
        data = line.split()
        if len(data) == 0:
            continue
        if nv < 0:
            nv = int(data[0])
            rc2 = float(data[1])
            rc3 = float(data[2])
            continue
        else:
            iv += 1
            if iv > nv:
                break
            if 'hard-limit' in options.keys() and options['hard-limit']:
                vs.append(float(data[0]))
                vrs.append([ float(data[1]), float(data[2])])
                vrsh.append([float(data[3]), float(data[4])])
                print(' iv,vrhmin,vrhmax= {0:3d} {1:11.3e} {2:11.3e}'.format(iv,
                                                                             float(data[3]),
                                                                             float(data[4])))
            else:
                vs.append(float(data[0]))
                vrs.append([ float(data[1]), float(data[2])])
                vrsh.append([-1e+30, 1e+30])
    vs = np.array(vs)
    vrs = np.array(vrs)
    vrsh = np.array(vrsh)
    print('')
    return rc2,rc3,vs,vrs,vrsh,options
    

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
    adf_triplets = kwargs['adf_triplets']

    rs = []
    rdfs = []
    ths = []
    adfs = []
    vol = 0.0
    a=b=c=alp=bet=gmm= 0.0
    if kwargs['rdf_match']:
        rs,rdfs = read_rdf(basedir+'/data.{0:s}.rdf'.format(prefix),specorder,rdf_pairs)
    if kwargs['adf_match']:
        ths,adfs = read_adf(basedir+'/data.{0:s}.adf'.format(prefix),specorder,adf_triplets)
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

def read_data(fname,):
    """
    General routine of reading data.
    
    Input file format
    -----------------
    ```
    #  Comment lines begins with '#' or '!'
    # 
    10    1.0
    0.1234  0.2345  0.3456  0.4567  0.5678  0.6789
    0.7890  0.8901  0.9012  0.0123
    ```
    - 1st line:  num of data (NDAT),  weight of the data (WDAT)
    - 2nd line-: data values (number of data should be equal to NDAT)
    """
    if not os.path.exists(fname):
        raise RuntimeError('File not exsits: ',fname)
    
    with open(fname,'r') as f:
        lines = f.readlines()

    ndat = 0
    wdat = 0.0
    data = None
    idat = 0
    done = False
    for line in lines:
        if line[0] in ('#','!'):
            continue
        ldat = line.split()
        if ndat < 1:
            ndat = int(ldat[0])
            wdat = float(ldat[1])
            data = np.zeros(ndat)
        else:
            if data is None:
                raise RuntimeError('data is None, which should not happen.')
            for i,d in enumerate(ldat):
                data[idat] = float(d)
                idat += 1
                if idat == ndat:
                    done = True
                    break
        if done:
            break
    return {'ndat':ndat, 'wdat':wdat, 'data':data}

def get_data2(basedir,prefix='ref',**kwargs):
    """
    New implementation of get_data, which loads data to be used to fit parameters.
    The prefix should be either ref or pmd.
    """

    matches = kwargs['match']
    # print('matches=',matches)

    data = {}
    for m in matches:
        fname = basedir+'/data.{0:s}.{1:s}'.format(prefix,m)
        # print('m,fname=',m,fname)
        try:
            data[m] = read_data(fname,)
        except:
            data[m] = None
            pass
    return data

def loss_func2(pmddata,eps=1.0e-8,**kwargs):
    """
    Compute loss function value using general get_data2 func.
    """
    refdata = kwargs['refdata']
    losses = {}
    L = 0.0
    misval = kwargs['missing_value']
    luplim = kwargs['fval_upper_limit']
    for name in refdata.keys():
        ref = refdata[name]
        wgt = ref['wdat']
        pmd = pmddata[name]
        if pmd == None:
            losses[name] = misval
            L += losses[name] *wgt
            continue
        num = ref['ndat']
        refd = ref['data']
        pmdd = pmd['data']
        z2 = 0.0
        sumdiff2 = 0.0
        for n in range(num):
            # print('n=',n)
            diff = pmdd[n] -refd[n]
            sumdiff2 += diff*diff
            z2 += refd[n]*refd[n]
        losses[name] = min(sumdiff2 /(z2+eps), luplim)
        L += losses[name] *wgt
        
    if kwargs['print_level'] > 0:
        print('   iid,losses= {0:8d}'.format(kwargs['iid']),end='')
        for k in losses.keys():
            loss = losses[k]
            print(' {0:10.4f}'.format(loss),end='')
        print(' {0:11.5f}'.format(L),flush=True)
    return L

def loss_func(pmddata,eps=1.0e-8,**kwargs):
    """
    Compute loss function value from reference and pmd data.
    """
    refdata = kwargs['refdata']
    wgts = kwargs['weights']
    rdf_pairs = kwargs['rdf_pairs']
    adf_triplets = kwargs['adf_triplets']

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
            Lr += diff2sum/(z+eps) #/len(rs)
        Lr /= len(rdf_pairs)

    #...ADF
    Lth = 0.0
    if kwargs['adf_match'] and len(adf_triplets) != 0:
        ths = refdata['ths']
        adfs_ref = refdata['adfs']
        adfs_pmd = pmddata['adfs']
        for t in adf_triplets:
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
            Lth += diff2sum /(z+eps) #/len(ths)
        Lth /= len(adf_triplets)

    #...volume
    Lvol = 0.0
    if kwargs['vol_match']:
        vol_ref = refdata['vol']
        vol_pmd = pmddata['vol']
        diff = vol_pmd -vol_ref
        Lvol = diff*diff /(vol_ref*vol_ref+eps)

    #...lattice parameters
    Llat = 0.0
    if kwargs['lat_match']:
        a0,b0,c0,alp0,bet0,gmm0 = refdata['lat']
        # print('refdata=',refdata['lat'])
        a,b,c,alp,bet,gmm = pmddata['lat']
        # print('pmddata=',pmddata['lat'])
        #...Need to take into account the definition difference bet/ dump and vasp
        # a0l,b0l,c0l,alp0l,bet0l,gmm0l = lat_vasp2dump(a0,b0,c0,alp0,bet0,gmm0)
        diff = ((a0-a)/a0)**2 \
            +((b0-b)/b0)**2 \
            +((c0-c)/c0)**2 \
            +((alp0-alp)/alp0)**2 \
            +((bet0-bet)/bet0)**2 \
            +((gmm0-gmm)/gmm0)**2
        # diff /= 6
        Llat = diff

    lrw = Lr*wgts['rdf']
    lthw = Lth*wgts['adf']
    lvolw = Lvol*wgts['vol']
    llatw = Llat*wgts['lat']
    L = lrw +lthw +lvolw +llatw

    if kwargs['print_level'] > 0:
        print('   iid,Lr,Lth,Lvol,Llat,L= {0:8d}'.format(kwargs['iid'])
              +'{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f}'.format(lrw,lthw,lvolw,llatw,L),
              flush=True)
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
    subjobscript = kwargs['subjob-script']
    subdir = kwargs['subdir-prefix'] +'{0:03d}'.format(kwargs['index'])
    print_level = kwargs['print_level']
    # print('refdata=',refdata)
    # print('pairs=',pairs)
    # print('triplets=',triplets)

    #...Create in.params.XXX files in each subdir
    varsfp = {}
    varsfp['rc2'] = kwargs['rc2']
    varsfp['rc3'] = kwargs['rc3']
    varsfp['variables'] = variables
    varsfp['vranges'] = vranges
    cwd = os.getcwd()
    if not os.path.exists(subdir):
        os.mkdir(subdir)
        shutil.copy(subjobscript,subdir+'/')
    os.chdir(subdir)

    if len(kwargs['param_files']) != 0:
        fp2params(varsfp['variables'], **kwargs)
    else:
        if 'vids' not in kwargs.keys():
            print(kwargs.keys())
        if kwargs['potential'] == 'BVSx':
            fp2BVSx(varsfp, **kwargs)
        elif kwargs['potential'] == 'BVS':
            fp2BVS(varsfp, **kwargs)
        elif kwargs['Morse'] == 'Morse':
            fp2Morse(varsfp, **kwargs)

    #...Compute pmd in the subdir_###
    L_up_lim = kwargs['fval_upper_limit']
    if print_level > 1:
        print('Running pmd and post-processing at '+subdir, flush=True)
    try:
        cmd = "./{0:s} > log.iid_{1:d}".format(subjobscript,kwargs['iid'])
        # subprocess.run(cmd.split(),check=True)
        # print('subdir,cmd=',subdir,cmd)
        subprocess.run(cmd,shell=True,check=True)
        os.chdir(cwd)
        # print('Going to get_data from ',subdir)
        if len(kwargs['match']) != 0:
            pmddata = get_data2(subdir,prefix='pmd',**kwargs)
            L = loss_func2(pmddata,**kwargs)
        else:
            pmddata = get_data(subdir,prefix='pmd',**kwargs)
            L = min( loss_func(pmddata,**kwargs), L_up_lim )
    except Exception as e:
        if print_level > 1:
            print('  Since pmd or post-process failed at {0:s}, '.format(subdir)
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
    # print(a,b,c,alp,bet,gmm)
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
    # print('alpr,cos(alpr)=',alpr,cos(alpr))
    # print('betr,cos(betr)=',betr,cos(betr))
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
    except Exception as e:
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

    start = time.time()

    nproc = int(args['--nproc'])

    infp = read_in_fitpot('in.fitpot')
    write_info(infp,args)

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
    adf_triplets = infp['adf_triplets']
    triplets = get_triplets(infp['interactions'])
    rc2,rc3,vs,vrs,vrsh,options = read_vars_fitpot(infp['param_file'])

    kwargs = infp
    kwargs['options'] = options
    kwargs['hardlim'] = vrsh
    # kwargs['infp'] = infp
    kwargs['rc2'] = rc2
    kwargs['rc3'] = rc3
    kwargs['pairs'] = pairs
    kwargs['rdf_pairs'] = rdf_pairs
    kwargs['triplets'] = triplets
    kwargs['adf_triplets'] = adf_triplets
    kwargs['subdir-prefix'] = args['--subdir-prefix']
    kwargs['subjob-script'] = args['--subjob-script']
    kwargs['start'] = start
    
    smpldir = infp['sample_directory']
    if len(infp['match']) != 0:
        refdata = get_data2(smpldir,prefix='ref',**kwargs)
    else:
        refdata = get_data(smpldir,prefix='ref',**kwargs)
        if kwargs['lat_match']:
            a0,b0,c0,alp0,bet0,gmm0 = refdata['lat']
            #...Need to take into account the definition difference bet/ dump and vasp
            a,b,c,alp,bet,gmm = lat_vasp2dump(a0,b0,c0,alp0,bet0,gmm0)
            # print('Reference lattice parameters:',a,b,c,alp,bet,gmm)
            refdata['lat'] = (a,b,c,alp,bet,gmm)
    kwargs['refdata'] = refdata

    if len(kwargs['param_files']) != 0: # New version of treating in.params.XXX files
        for fname in kwargs['param_files']:
            with open(fname,'r') as f:
                kwargs[fname] = f.read()
    else:
        fbvs, rads, vids, npqs, charges = read_params_Coulomb('in.params.Coulomb')
        kwargs['fbvs'] = fbvs
        kwargs['rads'] = rads
        kwargs['vids'] = vids
        kwargs['npqs'] = npqs
        kwargs['charges'] = charges

    print(' # iid,losses=      iid',end='')
    if len(kwargs['match']) > 0:
        for m in kwargs['match']:
            print('  {0:>9s}'.format(m),end='')
    else:
        if kwargs['rdf_match']: print('  {0:>9s}'.format('rdf'),end='')
        if kwargs['adf_match']: print('  {0:>9s}'.format('adf'),end='')
        if kwargs['vol_match']: print('  {0:>9s}'.format('vol'),end='')
        if kwargs['lat_match']: print('  {0:>9s}'.format('lat'),end='')
    print('      total')

    maxiter = kwargs['num_iteration']
    if kwargs['fitting_method'] in ('de','DE'):
        N = infp['de_num_individuals']
        F = infp['de_fraction']
        T = infp['de_temperature']
        CR = infp['de_crossover_rate']
        opt = DE(N,F,CR,T, vs,vrs,vrsh, func_wrapper, write_vars_fitpot,
                 nproc=nproc, **kwargs)
    elif kwargs['fitting_method'] in ('cs','CS','cuckoo','Cuckoo'):
        N = infp['cs_num_individuals']
        F = infp['cs_fraction']
        opt = CS(N,F, vs,vrs,vrsh, func_wrapper, write_vars_fitpot,
                 nproc=nproc, **kwargs)

    
    opt.run(maxiter)

    print('elapsed time = {0:f} sec.'.format(time.time()-start))
    
    return None

def headline():
    print('')
    print(' fp.py --- fit parameters to any target property ---')
    print('')
    cmd = ' '.join(s for s in sys.argv)
    print('   Executed as {0:s}'.format(cmd))
    hostname = subprocess.run(['hostname',], stdout=subprocess.PIPE).stdout.decode('utf-8')
    print('            on {0:s}'.format(hostname.strip()))
    print('            at {0:s}'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    print('   Please cite:')
    print('     1) R. Kobayashi, J. Open Source Software, 6(57), 2768 (2021)')
    print('     2) R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, APL Materials 8, 081111 (2020)')
    print()
    return None

if __name__ == "__main__":

    args = docopt(__doc__)

    headline()
    main(args)
