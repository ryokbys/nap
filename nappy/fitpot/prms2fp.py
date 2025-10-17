#!/usr/bin/env python
"""
Convert pmd param file `in.params.XXX` to fitpot var file `in.vars.fitpot`.
By default `in.vars.fitpot` includes hard-limit of parameters.
If repulsive pairs are given, `in.vars.conditions` is also written.

Usage:
  prms2fp.py [options] POTENTIAL

Options:
  -h, --help  Show this message and exit.
  -v          Verbose output using icecream.
  --rc RC     Cutoff radius for pair interaction. [default: 6.0]
  --rc3 RC3   Cutoff radius for angular interaction. [default: 3.0]
  -o OUTFNAME
              Output file name. [default: in.vars.fitpot_YYMMDD]
  --specorder SPECORDER
              Species order in comma-seperated format, e.g.) Li,P,O. [default: None]
  --repul-pairs PAIRS
              Pairs that should be repulsive. Hyphen-connected, comma-separated, e.g.) Li-O,P-O  [default: None]
"""
import os
from docopt import docopt
import numpy as np
from datetime import datetime
from icecream import ic
ic.disable()

from uf3util import read_params_uf3, read_params_uf3l

__author__ = "Ryo KOBAYASHI"
__version__ = "241119"

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
                f.write('  {0:12.4f}  {1:12.4e}  {2:12.4e} {3:12.4e} {4:12.4e}\n'.format(v,*vr,*hl))
            else:
                f.write('  {0:12.4f}  {1:12.4e}  {2:12.4e}\n'.format(v,*vr))
    print(' --> '+outfname)
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

def uf32fp(outfname,specorder):
    """
    Create in.vars.fitpot file from parameter file for uf3 potential.
    Cut-off radii for 2- and 3-body are given in the parameter file.
    """
    uf3_prms = read_params_uf3('in.params.uf3')

    fpvars = []
    vranges= []

    for spi,erg in uf3_prms['1B'].items():
        fpvars.append(erg)
        vranges.append((-1e+10, 1e+10))

    rc2max = 0.0
    d2b = uf3_prms['2B']
    for pair in d2b.keys():
        ncoef = d2b[pair]['ncoef']
        coefs = d2b[pair]['coefs']
        nlead = d2b[pair]['nlead']
        ntrail= d2b[pair]['ntrail']
        rc2max = max(rc2max, d2b[pair]['rc2b'])
        for i in range(ncoef):
            fpvars.append(coefs[i])
            if i < nlead or i >= ncoef -ntrail:
                vranges.append((0.0, 0.0))
            else:
                vranges.append((-1e+10, 1e+10))

    rc3max = 0.0
    d3b = uf3_prms['3B']
    for trio in d3b.keys():
        print(trio)
        ncij = d3b[trio]['ncij']
        ncik = d3b[trio]['ncik']
        ncjk = d3b[trio]['ncjk']
        rc3max = max(rc3max, d3b[trio]['rcij'], d3b[trio]['rcik'])
        coefs = d3b[trio]['coefs']
        nlead = d3b[trio]['nlead']
        ntrail= d3b[trio]['ntrail']
        for i in range(ncij):
            for j in range(ncik):
                for k in range(ncjk):
                    fpvars.append(coefs[i,j,k])
                    if i<nlead or i>=ncij-ntrail \
                       or j<nlead or j>=ncik-ntrail \
                       or k<nlead or k>=ncjk-ntrail:
                        vranges.append((0.0, 0.0))
                    else:
                        vranges.append((-1e+10, 1e+10))
    write_vars_fitpot(outfname, fpvars, vranges, rc2max, rc3max)
    return None

def uf3l2fp(outfname, specorder, repul_pairs=[]):
    """
    Create in.vars.fitpot file from parameter file for uf3l potential.
    Cut-off radii for 2- and 3-body are given in the parameter file.
    """
    uf3l_prms = read_params_uf3l('in.params.uf3l')

    fpvars = []
    vranges = []

    for spi, erg in uf3l_prms['1B'].items():
        fpvars.append(erg)
        vranges.append((-1e+10, 1e+10))

    rc2max = 0.0
    d2b = uf3l_prms['2B']
    for pair in d2b.keys():
        ncoef = d2b[pair]['ncoef']
        coefs = d2b[pair]['coefs']
        ntrail = d2b[pair]['ntrail']
        rc2max = max(rc2max, d2b[pair]['rc2b'])
        vmin = -1e+10
        if any( set(pair) == set(rp) for rp in repul_pairs ):
            vmin = 0.0
        for i in range(ncoef-ntrail):
            fpvars.append(coefs[i])
            vranges.append((vmin, 1e+10))
        for i in range(ncoef-ntrail, ncoef):
            fpvars.append(coefs[i])
            vranges.append((0.0, 0.0))

    rc3max = 0.0
    d3b = uf3l_prms['3B']
    for trio in d3b.keys():
        print(trio)
        ncoef = d3b[trio]['ncoef']
        rc3max = max(rc3max, d3b[trio]['rc'])
        coefs = d3b[trio]['coefs']
        # gmj, gmk are parameters in exp.
        # Too small values cause violation of energy conservation,
        # and too large values make them negligibly small.
        # Here, we set the range at (0.5, 2.0).
        fpvars.append(d3b[trio]['gmj'])
        vranges.append((0.5, 2.0)) # gmj, gmk should be greater than 0.0
        fpvars.append(d3b[trio]['gmk'])
        vranges.append((0.5, 2.0))
        for i in range(ncoef):
            fpvars.append(coefs[i])
            vranges.append((0.0, 1e+10))
    write_vars_fitpot(outfname, fpvars, vranges, rc2max, rc3max)
    if repul_pairs != []:
        write_vars_conditions(uf3l_prms, repul_pairs)
    return None


def write_vars_conditions(uf3l_prms, repul_pairs):
    """
    Write in.vars.conditions for vars2cond if repul_pairs != [].
    """
    msgline = []
    msgline.append('# var-ID-LHS,  operator,  var-ID-RHS\n')
    d2b = uf3l_prms['2B']
    inc = 0
    nconds = 0
    for spi, erg in uf3l_prms['1B'].items():
        inc += 1
    for pair in d2b.keys():
        ncoef = d2b[pair]['ncoef']
        ntrail = d2b[pair]['ntrail']
        for i in range(ncoef-ntrail):
            inc += 1
            if any( set(pair) == set(rp) for rp in repul_pairs ):
                if i != 0:
                    msgline.append(f'   {inc:d}  <  {inc-1:d}\n')
                    nconds += 1
        for i in range(ncoef-ntrail, ncoef):
            inc += 1
    wgt = 1.0
    msgline.insert(1, f'  {nconds:d}  {wgt:.3f}\n')
    print(' --> in.vars.conditions')
    with open('in.vars.conditions', 'w') as f:
        for l in msgline:
            f.write(l)
    return None


def main():
    args = docopt(__doc__)
    if args['-v']:
        ic.enable()
    potname = args['POTENTIAL']
    outfname = args['-o']
    outfname = outfname.replace('YYMMDD', datetime.now().strftime('%y%m%d'))
    rc = float(args['--rc'])
    rc3 = float(args['--rc3'])
    specorder = args['--specorder'].split(',')
    if specorder[0] == 'None':
        raise Exception('specorder must be specified.')
    repul_pairs0 = args['--repul-pairs'].split(',')
    repul_pairs = []
    if repul_pairs0[0] != 'None':
        for rp in repul_pairs0:
            spi, spj = rp.split('-')
            repul_pairs.append((spi, spj))
    print(' repulsive pairs = ', repul_pairs)

    ic(args)
    if potname == 'Morse':
        Morse2fp(outfname, specorder, rc, rc3)

    elif potname in ('BVS', 'BVSx'):
        """
        The term 'BVS' means that the in.var.fitpot file contains
        fbvs, species radius and Morse parameters,
        Thus in this case, specorder should be specified.
        In case of 'BVSx' contains angular parameters, so triplets should be specified.
        """
        BVS2fp(outfname, specorder, rc, rc3)

    elif potname in ('UF3', 'uf3'):
        """
        Ultra-fast force-field.
        """
        print(' rc, rc3 are given from in.params.uf3, even if --rc or --rc3 is given.')
        uf32fp(outfname, specorder)

    elif potname in ('UF3L', 'uf3l'):
        """
        UF3L (light)
        """
        print(' rc, rc3 are given from in.params.uf3l, even if --rc or --rc3 is given.')
        uf3l2fp(outfname, specorder, repul_pairs=repul_pairs)

    return None

if __name__ == "__main__":

    main()
