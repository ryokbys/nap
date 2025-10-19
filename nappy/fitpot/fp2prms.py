#!/usr/bin/env python
"""
Convert fitpot parameters `in.vars.fitpot` to pmd param file `in.params.XXX`.

Usage:
  fp2prms.py [options] <pot-type> <in-file>

Options:
  -h, --help  Show this message and exit.
  --specorder SPECORDER
              Species order in comma-separated format, e.g.) Li,P,O. [default: None]
  --pairs PAIRS
              Specify pairs for pair potential except for Coulomb by hyphen-connected and comma-separated,
              e.g.) Li-O,P-O. [default: None]
  --triplets TRIPLETS
              Specify triplets by hyphen-connected and comma-separated,
              e.g.) Li-O,P-O.
              The 1st species (X1 in X1-X2-X3) is the center of two bonds. [default: None]
"""
import os

from docopt import docopt
#from prms2fp import read_params_uf3
from uf3util import read_params_uf3, write_params_uf3, read_params_uf3l, write_params_uf3l
from datetime import datetime

__author__ = "Ryo KOBAYASHI"
__version__ = "241120"

def read_in_fitpot2(infname='in.fitpot'):
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
    param_files = []
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
        elif data[0] == 'param_files':
            param_files = [ name for name in data[1:] ]
            mode = None
        else:
            if mode == 'interactions':
                if len(data) not in (2,3):
                    raise Exception('len(data) is not 2 nor 3.')
                interact.append(data)
                if len(interact) == num_interact:
                    mode = None
            else:
                mode = None

    return specorder, interact, param_files


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
    charges = 'None'
    for line in lines:
        data = line.split()
        if len(data) == 0:
            mode = None
            continue
        if line[0] in ('!','#'):
            continue
        if data[0] == 'charges':
            if not data[1] in ('fixed_bvs','fixed'):
                raise ValueError('charges should be fixed_bvs in '+infname)
            mode = 'charges:'+data[1]
            charges = data[1]
            continue
        elif data[0] == 'fbvs':
            fbvs = float(data[1])
            mode = None
            continue
        elif data[0] == 'interactions':
            mode = None
            continue
        elif data[0] == 'rad_screened_cut':
            csp = data[1]
            rad = float(data[2])
            rads[csp] = rad
            mode = None
            continue
        elif mode is not None and 'charges' in mode:
            if '_bvs' in mode:
                if len(data) != 4:
                    raise ValueError('format of {0:s} seems wrong.'.format(infname))
                csp = data[0]
                vid = float(data[1])
                rad = float(data[2])
                npq = int(data[3])
                vids[csp] = vid
                rads[csp] = rad
                npqs[csp] = npq
            elif 'fixed' in mode:
                if len(data) != 2:
                    raise ValueError('format of {0:s} seems wrong.'.format(infname))
                csp = data[0]
                vid = float(data[1])
                vids[csp] = vid
        else:
            pass
    return fbvs,rads,vids,npqs,charges

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

def write_params_Coulomb(outfname,specorder,pairs,fbvs,rads,
                         vids=None,npqs=None,charges=None):
    """
    Write in.params.Coulomb specific for BVS.
    """

    with open(outfname,'w') as f:
        f.write('terms   screened_cut\n')
        f.write('charges   {0:s}\n'.format(charges))
        if charges == 'fixed_bvs':
            if vids is None or npqs is None:
                raise ValueError('vids and npqs should not be None in the case of charges==fixed_bvs.')
            for s in specorder:
                f.write('  {0:3s}  {1:4.1f}  {2:7.4f}  {3:2d}\n'.format(s,vids[s],
                                                                        rads[s],npqs[s]))
        elif charges == 'fixed':
            if vids is None:
                raise ValueError('vids should not be None in the case of charges==fixed.')
            for s in specorder:
                f.write('  {0:3s}  {1:7.4f}\n'.format(s,vids[s]))
            f.write('\n')
            for s in specorder:
                f.write('rad_screened_cut  {0:3s}  {1:7.4f}\n'.format(s,rads[s]))
        else:
            raise ValueError('No charges specified.')
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
    return None

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
    nv = len(vs)

    specorder = kwargs['specorder']
    pairs = kwargs['pairs']

    #...Check num of vars
    if nv != len(pairs)*3 +len(specorder) +1:
        raise ValueError('Number of variables is wrong! nv,len(pairs),len(specorder)=',
                         nv,len(pairs),len(specorder))

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
        charges = kwargs['charges']
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads,
                             vids=vids0,npqs=npqs0,charges=charges)
    else:
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads)

    return None

def fp2BVSx(varsfp, **kwargs):

    rc2 = varsfp['rc2']
    rc3 = varsfp['rc3']
    vs = varsfp['variables']
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
        charges = kwargs['charges']
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads,
                             vids=vids0,npqs=npqs0,charges=charges)
    else:
        print('no vids?')
        print(kwargs.keys())
        write_params_Coulomb('in.params.Coulomb',specorder,pairs,fbvs,rads)

    return None

def fp2uf3(outfname, vs, uf3_prms, **kwargs):

    iv = -1
    for spi in uf3_prms['1B'].keys():
        iv += 1
        uf3_prms['1B'][spi] = vs[iv]

    for pair in uf3_prms['2B'].keys():
        dic = uf3_prms['2B'][pair]
        for i in range(len(dic['coefs'])):
            iv += 1
            dic['coefs'][i] = vs[iv]
        uf3_prms['2B'][pair] = dic

    for trio in uf3_prms['3B'].keys():
        dic = uf3_prms['3B'][trio]
        ncij = dic['ncij']
        ncik = dic['ncik']
        ncjk = dic['ncjk']
        for i in range(ncij):
            for j in range(ncik):
                for k in range(ncjk):
                    iv += 1
                    dic['coefs'][i,j,k] = vs[iv]
        uf3_prms['3B'][trio] = dic

    write_params_uf3(uf3_prms,
                     outfname=outfname,
                     overwrite=True)
    print(f' --> {outfname}')
    return None

def fp2uf3l(outfname, vs, uf3l_prms, **kwargs):

    iv = -1
    for spi in uf3l_prms['1B'].keys():
        iv += 1
        uf3l_prms['1B'][spi] = vs[iv]

    for pair in uf3l_prms['2B'].keys():
        dic = uf3l_prms['2B'][pair]
        ncoef = dic['ncoef']
        nlead = dic['nlead']
        ntrail= dic['ntrail']
        for i in range(ncoef):
            iv += 1
            dic['coefs'][i] = vs[iv]
        uf3l_prms['2B'][pair] = dic

    for trio in uf3l_prms['3B'].keys():
        dic = uf3l_prms['3B'][trio]
        ncoef = dic['ncoef']
        iv += 1
        dic['rcij'] = vs[iv]
        iv += 1
        dic['rcik'] = vs[iv]
        iv += 1
        dic['gmj'] = vs[iv]
        iv += 1
        dic['gmk'] = vs[iv]
        for i in range(ncoef):
            iv += 1
            dic['coefs'][i] = vs[iv]
        uf3l_prms['3B'][trio] = dic

    write_params_uf3l(uf3l_prms,
                      outfname=outfname,
                      overwrite=True)
    print(f' --> {outfname}')
    return None

def fp2params(vs,**kwargs):
    """
    Conversion from fp-vars to files specified in param_files in in.fitpot.
    The param_files should contain key-phrases such as '{p[0]}' that are converted from fp-vars,
    and the indices in the key-phrases must correspond to those in fp-vars..
    """
    import re
    try:
        param_fnames = kwargs['param_files']
    except:
        raise
    if type(param_fnames) != list or len(param_fnames) < 1:
        raise ValueError('MD_param_files may not be specified correctly...')

    for fname in param_fnames:
        try:
            fcontents = kwargs[fname]
            #...If the format is '{0:.1f}'-style, replace them to '{p[0]:.2f}'-style
            res = re.search(r'\{[0-9]+:',fcontents)
            if res != None:
                fcontents = re.sub(r'\{([0-9]+):',r'{p[\1]:',fcontents)
            new_contents = fcontents.format(p=vs)
        except:
            print('ERROR: Failed to replace the parameters in param_files !!!')
            print(fcontents)
            raise
        with open(fname,'w') as f:
            f.write(new_contents)
    return None

def main():
    import os,sys
    from nappy.fitpot.fp import read_in_fitpot
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])))
    # print(args)
    infname = args['<in-file>']
    pot_type = args['<pot-type>']
    pairs = args['--pairs'].split(',')
    pairs = [ pair.split('-') for pair in pairs ]
    triplets = args['--triplets'].split(',')
    triplets = [ t.split('-') for t in triplets ]
    specorder = args['--specorder'].split(',')

    if specorder[0] == 'None':
        try:
            #specorder, interact, param_files = read_in_fitpot('in.fitpot')
            infp = read_in_fitpot('in.fitpot')
            interact = infp['interactions']
            specorder = infp['specorder']
            pairs = []
            triplets = []
            for i in interact:
                if len(i) == 2:
                    pairs.append(i)
                elif len(i) == 3:
                    triplets.append(i)
            print(' Loaded some information from in.fitpot')
        except:
            raise Exception('Something wrong with reading in.fitpot.')
    else:
        try:
            #specorder, interact, param_files = read_in_fitpot('in.fitpot')
            infp = read_in_fitpot('in.fitpot')
            # if specorder[0] == 'None':
            #     specorder = infp['specorder']
            print(' Loaded some info from in.fitpot')
        except:
            raise Exception('Something wrong with reading in.fitpot.')

    # pairs = sort_pairs(pairs,specorder)

    kwargs = infp
    varsfp = read_vars_fitpot(infname)
    if pot_type == 'map':
        """
        The keyword map specifies that the in.params.XXX files contain entries such as '{p[0]}'
        that are to be replaced.
        """
        param_files = infp['param_files']
        for pfile in param_files:
            with open(pfile,'r') as f:
                kwargs[pfile] = f.read()
            os.system(f'cp -f {pfile} {pfile}.bak')
        fp2params(varsfp['variables'], **kwargs)
        #print(' Wrote the following files:')
        for pfile in param_files:
            print(f' --> {pfile}')
        print('')

    elif pot_type in ('UF3', 'uf3'):
        prmfname = 'in.params.uf3'
        if not os.path.exists(prmfname):
            raise Exception(f'Cannot find {prmfname}.')
        uf3_prms = read_params_uf3(prmfname)
        outfname = prmfname +'_'+datetime.now().strftime('%y%m%d')
        fp2uf3(outfname, varsfp['variables'], uf3_prms, **kwargs)
    elif pot_type in ('UF3L', 'uf3l'):
        prmfname = 'in.params.uf3l'
        if not os.path.exists(prmfname):
            raise Exception(f'Cannot find {prmfname}.')
        uf3l_prms = read_params_uf3l(prmfname)
        outfname = prmfname +'_'+datetime.now().strftime('%y%m%d')
        fp2uf3l(outfname, varsfp['variables'], uf3l_prms, **kwargs)
    else:
        if len(pairs) == 0:
            raise ValueError('Pairs must be specified.')
        print(' Pairs to be extracted:')
        for pair in pairs:
            print('   {0:s}-{1:s}'.format(*pair))
        if len(triplets) != 0:
            print(' Triplets to be extracted:')
            for t in triplets:
                print('   {0:s}-{1:s}-{2:s}'.format(*t))

                kwargs['pairs'] = pairs
        if pot_type == 'Morse':
            fp2Morse(varsfp, **kwargs)
            print(' --> in.params.Morse')

        elif pot_type == 'BVS':
            """
            The term 'BVS' means that the in.var.fitpot file contains
            fbvs, species radius and Morse parameters.
            Thus in this case, specorder should be specified.
            """
            #...If there is an old in.params.Coulomb file, get Vid and npq from it
            try:
                fbvs0,rads0,vids0,npqs0,charges = read_params_Coulomb('in.params.Coulomb')
                kwargs['vids'] = vids0
                kwargs['npqs'] = npqs0
                kwargs['charges'] = charges
            except Exception as e:
                pass
            fp2BVS(varsfp, **kwargs)
            print(' --> in.params.{Morse,Coulomb}')

        elif pot_type == 'BVSx':
            """
            The term 'BVSx' means that the in.var.fitpot file contains
            fbvs, species radius, Morse and angular parameters.
            Thus in this case, specorder, pairs and triplets should be specified.
            """
            kwargs['triplets'] = triplets
            #...If there is an old in.params.Coulomb file, get Vid and npq from it
            try:
                fbvs0,rads0,vids0,npqs0,charges = read_params_Coulomb('in.params.Coulomb')
                kwargs['vids'] = vids0
                kwargs['npqs'] = npqs0
                kwargs['charges'] = charges
            except Exception as e:
                pass

            fp2BVSx(varsfp, **kwargs)
            print(' --> in.params.{Morse,Coulomb,angular}')
        else:
            print(f' No such pot_type: {pot_type}')
    return None

if __name__ == "__main__":
    main()
