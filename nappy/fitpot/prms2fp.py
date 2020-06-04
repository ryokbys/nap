#!/usr/bin/env python
"""
Convert pmd param file `in.params.XXX` to fitpot var file `in.vars.fitpot`.
By default `in.vars.fitpot` includes hard-limit of parameters.

Usage:
  prms2fp.py (Morse|BVSx|BVS) [options]

Options:
  -h, --help  Show this message and exit.
  --rc RC     Cutoff radius for pair interaction. [default: 6.0]
  --rc3 RC3   Cutoff radius for angular interaction. [default: 3.0]
  -o OUTFNAME
              Output file name. [default: in.vars.fitpot]
  --specorder SPECORDER
              Species order in comma-seperated format, e.g.) Li,P,O. [default: None]
"""
from __future__ import print_function

import os
from docopt import docopt

__author__ = "Ryo KOBAYASHI"
__version__ = "200425"

def write_vars_fitpot(outfname,fpvars,vranges,rc,rc3,hardlim=None):
    with open(outfname,'w') as f:
        #...By default, set hard-limit
        if hardlim is not None:
            f.write('! hard-limit: T\n')
            f.write('!\n')
        f.write('  {0:d}   {1:7.3f} {2:7.3f}\n'.format(len(fpvars),rc,rc3))
        for i in range(len(fpvars)):
            v = fpvars[i]
            vr = vranges[i]
            if hardlim is not None:
                hl = hardlim[i]
                f.write('  {0:9.4f}  {1:9.4f}  {2:9.4f} {3:8.3f} {4:8.3f}\n'.format(v,*vr,*hl))
            else:
                f.write('  {0:9.4f}  {1:9.4f}  {2:9.4f}\n'.format(v,*vr))
    print(' Wrote '+outfname)
    return None

def read_params_Coulomb(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    with open(infname,'r') as f:
        lines = f.readlines()

    mode = None
    fbvs = 0
    rads = {}
    for line in lines:
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if data[0] == 'charges':
            if not data[1] == 'fixed_bvs':
                raise ValueError('charges should be fixed_bvs in '+infname)
            mode = 'charges'
        elif data[0] == 'fbvs':
            fbvs = float(data[1])
        elif mode == 'charges':
            if len(data) != 4:
                raise ValueError('format of {0:s} seems wrong.'.format(infname))
            csp = data[0]
            rad = float(data[2])
            rads[csp] = rad

    return fbvs, rads

def read_params_Morse(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)
    
    with open(infname,'r') as f:
        lines = f.readlines()

    morse_prms = {}
    for l in lines:
        if l[0] in ('!','#'):
            continue
        data = l.split()
        if len(data) < 5:
            continue
        if not data[0].isdigit() and not data[1].isdigit():
            cspi = data[0]
            cspj = data[1]
            d0 = float(data[2])
            alp = float(data[3])
            rmin = float(data[4])
            morse_prms[(cspi,cspj)] = (d0,alp,rmin)
        else:
            continue
    return morse_prms

def read_params_angular(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)
    
    with open(infname,'r') as f:
        lines = f.readlines()

    angular_prms = {}
    for l in lines:
        if l[0] in ('!','#'):
            continue
        data = l.split()
        if len(data) < 8:
            continue
        if not data[0].isdigit() and not data[1].isdigit() \
           and not data[2].isdigit():
            atype = data[0]
            cspi = data[1]
            cspj = data[2]
            cspk = data[3]
            rc3 = float(data[4])
            alp = float(data[5])
            bet = float(data[6])
            gmm = float(data[7])
            angular_prms[(cspi,cspj,cspk)] = (alp,bet,gmm)
        else:
            continue
    return angular_prms
    
def prms_to_fpvars(specorder,prms):
    """
    Convert Morse/angular prms to fpvars taking specorder into account.

    Parameters
    ----------
    specorder : list
         Species order.
    prms : dict
         Dictionary of parameters. Each key stands for pair or triplet.

    Returns
    -------
    fpvars : list
         List of parameters.
    vranges : list
         List of parameter ranges.
    hardlim : list
         List of hard-limits of parameters.
    """

    keys = list(prms.keys())
    maxlen = -1
    for k in keys:
        maxlen = max(maxlen,len(k))

    fpvars = []
    vranges = []
    hardlim = []
    
    if maxlen == 2:  # 2-body
        for i in range(len(specorder)):
            si = specorder[i]
            for j in range(i,len(specorder)):
                sj = specorder[j]
                for k,vs in prms.items():
                    if same_pair((si,sj),k):
                        for v in vs:
                            fpvars.append(v)
                            vranges.append((max(v*.75,0.0),v*1.25))
                            hardlim.append((0.0, 10.0))
                        break
    
    elif maxlen == 3:  # 3-body
        for i in range(len(specorder)):
            si = specorder[i]
            for j in range(len(specorder)):
                sj = specorder[j]
                for k in range(j,len(specorder)):
                    sk = specorder[k]
                    for ks,vs in prms.items():
                        if same_triplet((si,sj,sk),ks):
                            for iv,v in enumerate(vs):
                                fpvars.append(v)
                                #...Code specific to angular1 potential
                                if iv == 2:
                                    vranges.append((max(-1.,v-0.5),min(1.0,v+0.5)))
                                    hardlim.append((-1., 1.))
                                elif iv == 1:
                                    vranges.append((1., 1.))
                                    hardlim.append((1., 1.))
                                else:
                                    vranges.append((max(v*0.75,0.0), v*1.25))
                                    hardlim.append((0.0, 10.0))
                            break
    return fpvars, vranges, hardlim
    

def Morse2fp(outfname,specorder,rc,rc3):

    morse_prms = read_params_Morse('in.params.Morse')

    fpvars, vranges = prms_to_fpvars(specorder,morse_prms)

    write_vars_fitpot(outfname,fpvars,vranges,rc,rc3)

    return None

def same_pair(pair1,pair2):
    a1,a2 = pair1
    b1,b2 = pair2
    if (a1==b1 and a2==b2) or (a1==b2 and a2==b1):
        return True
    else:
        return False

def same_triplet(triple1,triple2):
    a1,a2,a3 = triple1
    b1,b2,b3 = triple2
    if (a1==b1) and ((a2==b2 and a3==b3) or (a2==b3 and a3==b2)):
        return True
    else:
        return False

def BVS2fp(outfname,specorder,rc,rc3):

    morse_prms = read_params_Morse('in.params.Morse')
    fbvs, rads = read_params_Coulomb('in.params.Coulomb')

    fpvars = []
    vranges= []
    hardlim= []

    #...fbvs
    fpvars.append(fbvs)
    vranges.append((fbvs,fbvs))
    hardlim.append((fbvs,fbvs))
    #...rads
    for s in specorder:
        if s not in rads.keys():
            raise Exception('spcs not in rads.')
        rad = rads[s]
        fpvars.append(rad)
        vranges.append((rad*0.75,rad*1.25))
        hardlim.append((0.0, 3.0))
        
    #...Morse parameters
    v_morse, vr_morse, hl_morse = prms_to_fpvars(specorder,morse_prms)
    fpvars += v_morse
    vranges += vr_morse
    hardlim += hl_morse

    write_vars_fitpot(outfname,fpvars,vranges,rc,rc3,hardlim=hardlim)

    return None

def BVSx2fp(outfname,specorder,rc,rc3):

    morse_prms = read_params_Morse('in.params.Morse')
    fbvs, rads = read_params_Coulomb('in.params.Coulomb')
    angular_prms = read_params_angular('in.params.angular')

    fpvars = []
    vranges= []
    hardlim = []

    #...fbvs
    fpvars.append(fbvs)
    vranges.append((fbvs,fbvs))
    hardlim.append((fbvs,fbvs))
    #...rads
    for s in specorder:
        if s not in rads.keys():
            raise Exception('spcs not in rads.')
        rad = rads[s]
        fpvars.append(rad)
        vranges.append((rad*0.75,rad*1.25))
        hardlim.append((0.0,3.0))
    #...Morse parameters
    v_morse, vr_morse, hl_morse = prms_to_fpvars(specorder,morse_prms)
    fpvars += v_morse
    vranges += vr_morse
    hardlim += hl_morse
    #...angular parameters
    v_angular, vr_angular, hl_angl = prms_to_fpvars(specorder,angular_prms)
    fpvars += v_angular
    vranges += vr_angular
    hardlim += hl_angl
    write_vars_fitpot(outfname,fpvars,vranges,rc,rc3,hardlim=hardlim)

    return None



if __name__ == "__main__":

    args = docopt(__doc__)
    outfname = args['-o']
    rc = float(args['--rc'])
    rc3 = float(args['--rc3'])
    specorder = args['--specorder'].split(',')
    if specorder[0] == 'None':
        raise Exception('specorder must be specified.')

    if args['Morse']:
        Morse2fp(outfname,specorder,rc,rc3)

    elif args['BVS']:
        """
        The term 'BVS' means that the in.var.fitpot file contains 
        fbvs, species radius and Morse parameters.
        Thus in this case, specorder should be specified.
        """
        BVS2fp(outfname,specorder,rc,rc3)
        
    elif args['BVSx']:
        """
        The term 'BVSx' means that the in.var.fitpot file contains 
        fbvs, species radius, Morse and angular parameters.
        Thus in this case, specorder, pairs and triplets should be specified.
        """
        BVSx2fp(outfname,specorder,rc,rc3)
        
