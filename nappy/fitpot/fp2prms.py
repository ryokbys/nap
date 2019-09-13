#!/usr/bin/env python
"""
Convert fitpot parameters `in.vars.fitpot` to pmd param file `in.params.XXX`.

Usage:
  fp2prms.py (Morse|BVSx|BVS) [options] IN_VARS_FITPOT

Options:
  -h, --help  Show this message and exit.
  --specorder SPECORDER
              Species order in comma-seperated format, e.g.) Li,P,O. [default: None]
  --pairs PAIRS
              Specify pairs for pair potential except for Coulomb by hyphen-connected and comma-separated,
              e.g.) Li-O,P-O. [default: None]
  --triplets TRIPLETS
              Specify triplets by hyphen-connected and comma-separated,
              e.g.) Li-O,P-O. 
              The 1st species (X1 in X1-X2-X3) is the center of two bonds. [default: None]
"""
from __future__ import print_function

import os
from docopt import docopt

__author__ = "Ryo KOBAYASHI"
__version__ = "rev190827"

def read_in_fitpot(infname='in.fitpot'):
    """
    Get specorder, pairs, triplets from in.fitpot.
    """
    if not os.path.exists(infname):
        raise FileNotFoundError(infname)
    
    with open(infname,'r') as f:
        lines = f.readlines()

    mode = None
    specorder = []
    interact = []
    for line in lines:
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if data[0] in ('#','!'):
            mode = None
            continue
        elif data[0] == 'specorder':
            specorder = [ x for x in data[1:] ]
            continue
        elif data[0] == 'interactions':
            num_interact = int(data[1])
            mode = 'interactions'
            continue
        else:
            if mode == 'interactions':
                if len(data) not in (2,3):
                    raise Exception('len(data) is not 2 nor 3.')
                interact.append(data)
                if len(interact) == num_interact:
                    mode = None
            else:
                mode = None

    return specorder, interact
    

def read_params_Coulomb(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    with open(infname,'r') as f:
        lines = f.readlines()

    mode = None
    fbvs = 0
    rads = {}
    vids = {}
    npqs = {}
    for line in lines:
        if line[0] in ('!','#'):
            continue
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if data[0] == 'charges':
            if not data[1] == 'fixed_bvs':
                raise ValueError('charges should be fixed_bvs in '+infname)
            mode = 'charges'
            continue
        elif data[0] == 'fbvs':
            fbvs = float(data[1])
            mode = None
            continue
        elif data[0] == 'interactions':
            mode = None
            continue
        elif mode == 'charges':
            if len(data) != 4:
                raise ValueError('format of {0:s} seems wrong.'.format(infname))
            csp = data[0]
            vid = float(data[1])
            rad = float(data[2])
            npq = int(data[3])
            vids[csp] = vid
            rads[csp] = rad
            npqs[csp] = npq

    return fbvs,rads,vids,npqs

def read_vars_fitpot(fname='in.vars.fitpot'):
    """
    Read in.vars.fitpot and return data.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    fpvars = []
    vranges = []
    il = -1
    nv = -1
    while True:
        il += 1
        line = lines[il]
        if line[0] in ('!','#'):  # skip comment line
            il += 1
            continue
        data = line.split()
        if nv < 0:
            nv = int(data[0])
            rc = float(data[1])
            rc3= float(data[2])
        else:
            fpvars.append(float(data[0]))
            vranges.append([ float(x) for x in data[1:3]])
            if len(fpvars) == nv:
                break
    varsfp = {}
    varsfp['rc2'] = rc
    varsfp['rc3'] = rc3
    varsfp['variables'] = fpvars
    varsfp['vranges'] = vranges
    return varsfp

def write_params_Morse(outfname,pairs,morse_prms):
    
    with open(outfname,'w') as f:
        f.write('# cspi, cspj,    D,      alpha,  rmin\n')
        inc = -1
        for pair in pairs:
            d0,alp,rmin = morse_prms[tuple(pair)]
            f.write('  {0:3s}   {1:3s}'.format(*pair))
            f.write('  {0:7.4f}  {1:7.4f}  {2:7.4f}\n'.format(d0,alp,rmin))
            
    return None

def write_params_Coulomb(outfname,specorder,pairs,fbvs,rads,vids=None,npqs=None):
    """
    Write in.params.Coulomb specific for BVS.
    """
    
    with open(outfname,'w') as f:
        f.write('terms   screened_cut\n')
        f.write('charges   fixed_bvs\n')
        if vids is not None and npqs is not None:
            for s in specorder:
                f.write('  {0:3s}  {1:4.1f}  {2:7.4f}  {3:2d}\n'.format(s,vids[s],
                                                                        rads[s],npqs[s]))
        else:
            for s in specorder:
                f.write('  {0:3s}  Vid   {1:7.4f}  npq\n'.format(s,rads[s]))
        f.write('\n')
        f.write('fbvs    {0:7.3f}\n'.format(fbvs))
        f.write('\n')
        f.write('interactions\n')
        for i in range(len(specorder)):
            si = specorder[i]
            for j in range(i,len(specorder)):
                sj = specorder[j]
                if not same_pair_exists([si,sj],pairs):
                    f.write('   {0:3s}  {1:3s}\n'.format(si,sj))
    return None

def write_params_angular(outfname,triplets,angular_prms):
    with open(outfname,'w') as f:
        f.write('# type,   cspi, cspj, cspk,  rc3,   alp,   bet,   gmm\n')
        for t in triplets:
            rc3,alp,bet,gmm = angular_prms[tuple(t)]
            f.write(' angular1   {0:3s}   {1:3s}   {2:3s} '.format(*t))
            f.write(' {0:6.2f}  {1:7.3f} {2:7.3f} {3:7.3f}\n'.format(rc3,alp,bet,gmm))
    
def sort_pairs(pairs,specorder):

    sorted_pairs = []
    for i in range(len(specorder)):
        si = specorder[i]
        for j in range(i,len(specorder)):
            sj = specorder[j]
            if same_pair_exists((si,sj),pairs):
                sorted_pairs.append((si,sj))
    return sorted_pairs

def same_pair(pair1,pair2):
    a1,a2 = pair1
    b1,b2 = pair2
    if (a1==b1 and a2==b2) or (a1==b2 and a2==b1):
        return True
    else:
        return False

def same_pair_exists(pair,pairs):
    for p in pairs:
        if same_pair(p,pair):
            return True
    return False
    
def fp2Morse(varsfp, **kwargs):

    rc2 = varsfp['rc2']
    rc3 = varsfp['rc3']
    vs = varsfp['variables']
    vrs = varsfp['vranges']
    nv = len(vs)

    pairs = kwargs['pairs']
    
    #...Check num of vars and pairs
    if len(pairs)*3 != nv:
        raise ValueError('Number of variables and pairs are inconsistent.')

    morse_prms = {}
    inc = -1
    for p in pairs:
        inc += 1
        d0 = vs[inc]
        inc += 1
        alp = vs[inc]
        inc += 1
        rmin = vs[inc]
        morse_prms[p] = (d0,alp,rmin)

    write_params_Morse('in.params.Morse',pairs,morse_prms)
        
    return

def fp2BVS(varsfp, **kwargs):

    rc2 = varsfp['rc2']
    rc3 = varsfp['rc3']
    vs = varsfp['variables']
    vrs = varsfp['vranges']
    nv = len(vs)

    specorder = kwargs['specorder']
    pairs = kwargs['pairs']

    #...Check num of vars
    if nv != len(pairs)*3 +len(specorder) +1:
        raise ValueError('Number of variables is wrong! nv,len(pairs),len(specorder)='
                         ,nv,len(pairs),len(specorder))

    inc = 0
    fbvs = vs[inc]
    rads = {}
    for s in specorder:
        inc += 1
        rads[s] = vs[inc]

    morse_prms = {}
    for p in pairs:
        inc += 1
        d0 = vs[inc]
        inc += 1
        alp = vs[inc]
        inc += 1
        rmin = vs[inc]
        morse_prms[tuple(p)] = (d0,alp,rmin)

    write_params_Morse('in.params.Morse',pairs,morse_prms)

    if 'vids' in kwargs.keys():
        vids0 = kwargs['vids']
        npqs0 = kwargs['npqs']
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads,
                             vids=vids0,npqs=npqs0)
    else:
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads)

    return None

def fp2BVSx(varsfp, **kwargs):

    rc2 = varsfp['rc2']
    rc3 = varsfp['rc3']
    vs = varsfp['variables']
    vrs = varsfp['vranges']
    nv = len(vs)

    specorder = kwargs['specorder']
    pairs = kwargs['pairs']
    triplets = kwargs['triplets']

    #...Check num of vars
    if nv != len(pairs)*3 +len(specorder) +1 +len(triplets)*3:
        raise ValueError('Number of variables is wrong.')

    inc = 0
    fbvs = vs[inc]
    rads = {}
    for s in specorder:
        inc += 1
        rads[s] = vs[inc]

    morse_prms = {}
    for p in pairs:
        inc += 1
        d0 = vs[inc]
        inc += 1
        alp = vs[inc]
        inc += 1
        rmin = vs[inc]
        morse_prms[tuple(p)] = (d0,alp,rmin)

    angular_prms = {}
    for t in triplets:
        inc += 1
        alp = vs[inc]
        inc += 1
        bet = vs[inc]
        inc += 1
        gmm = vs[inc]
        angular_prms[tuple(t)] = (rc3,alp,bet,gmm)

    write_params_Morse('in.params.Morse',pairs,morse_prms)
    write_params_angular('in.params.angular',triplets,angular_prms)
    if 'vids' in kwargs.keys():
        vids0 = kwargs['vids']
        npqs0 = kwargs['npqs']
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads,
                             vids=vids0,npqs=npqs0)
    else:
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads)

    return None


if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['IN_VARS_FITPOT']
    pairs = args['--pairs'].split(',')
    pairs = [ pair.split('-') for pair in pairs ]
    triplets = args['--triplets'].split(',')
    triplets = [ t.split('-') for t in triplets ]
    specorder = args['--specorder'].split(',')

    if specorder[0] == 'None':
        try:
            specorder, interact = read_in_fitpot('in.fitpot')
            pairs = []
            triplets = []
            for i in interact:
                if len(i) == 2:
                    pairs.append(i)
                elif len(i) == 3:
                    triplets.append(i)
            print(' specorder, pairs and triplets are loaded from in.fitpot')
        except:
            raise Exception('specorder and pair must be specified or loaded.')
    
    if len(pairs) == 0:
        raise ValueError('Pairs must be specified.')
    print(' Pairs to be extracted:')
    for pair in pairs:
        print('   {0:s}-{1:s}'.format(*pair))
    if len(triplets) != 0:
        print(' Triplets to be extracted:')
        for t in triplets:
            print('   {0:s}-{1:s}-{2:s}'.format(*t))
    
    pairs = sort_pairs(pairs,specorder)

    kwargs = {}
    kwargs['specorder'] = specorder
    kwargs['pairs'] = pairs
    if args['Morse']:
        varsfp = read_vars_fitpot(infname)
        fp2Morse(varsfp, **kwargs)
        print(' Wrote in.params.Morse')

    elif args['BVS']:
        """
        The term 'BVS' means that the in.var.fitpot file contains 
        fbvs, species radius and Morse parameters.
        Thus in this case, specorder should be specified.
        """
        #...If there is an old in.params.Coulomb file, get Vid and npq from it
        try:
            fbvs0,rads0,vids0,npqs0 = read_params_Coulomb('in.params.Coulomb')
            kwargs['vids'] = vids0
            kwargs['npqs'] = npqs0
        except:
            pass
        varsfp = read_vars_fitpot(infname)
        fp2BVS(varsfp, **kwargs)
        print(' Wrote in.params.{Morse,Coulomb}')
        
    elif args['BVSx']:
        """
        The term 'BVSx' means that the in.var.fitpot file contains 
        fbvs, species radius, Morse and angular parameters.
        Thus in this case, specorder, pairs and triplets should be specified.
        """
        kwargs['triplets'] = triplets
        #...If there is an old in.params.Coulomb file, get Vid and npq from it
        try:
            fbvs0,rads0,vids0,npqs0 = read_params_Coulomb('in.params.Coulomb')
            kwargs['vids'] = vids0
            kwargs['npqs'] = npqs0
        except:
            pass
        
        varsfp = read_vars_fitpot(infname)
        fp2BVSx(varsfp, **kwargs)
        print(' Wrote in.params.{Morse,Coulomb,angular}')
        
