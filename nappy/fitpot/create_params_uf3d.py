#!/usr/bin/env python
"""
Create in.params.uf3d and in.vars.fitpot file.

Usage:
  {0:s} [options] CONFIG_UF3D_FILE

Options:
  -h, --help       Show this message and exit.
  -v, --verbose    Verbose output. [default: False]
  --author AUTHOR  Author of the param file. [default: None]
  --no-correct     Apply prior-knowledge correction. [default: False]
  --specorder SPECORDER
                   Species order in comma-seperated format, e.g.) Li,P,O. [default: None]
  --init3b INIT3B  Initial value of 3-body potential coefficients. [default: 1.0]
"""
import os, sys
from docopt import docopt
import numpy as np
from nappy.fitpot.uf3util import write_params_uf3d
from nappy.fitpot.prms2fp import write_vars_fitpot, write_vars_conditions
from nappy.fitpot.zbl import dddzbl
import yaml
from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "260208"


def load_yaml(file_path):
    import yaml
    with open(file_path, 'r', encoding='utf-8') as file:
        data = yaml.safe_load(file,)
    return data


def create_params(config,
                  outfname='in.params.uf3d',
                  author=None,
                  init3b=0.0):
    pairs = config['pairs']
    pair_cutoffs = [ p['cutoff'] for p in pairs ]
    trios = config['trios']
    trio_cutoffs = [ rc for t in trios for rc in t['cutoffs'] ]

    prms = {}
    prms['1B'] = []
    prms['2B'] = []
    prms['3B'] = []

    ## 1B
    species = config['species']
    for spi in species:
        d1b = {'species': spi,
               'erg': 0.0}
        prms['1B'].append(d1b)

    ## Params common to 2B and 3B
    leading_trim = 0
    trailing_trim = 3
    spacing = 'nk' # non-uniform knots

    ## 2B
    rmin = 0.1
    rmax = max(pair_cutoffs)
    #pair_res = config['pair_resolutions']
    #resolution = 20  # coef数はこれ＋３， knot数はこれ＋７
    for ip, p in enumerate(pairs):
        d2b = {'pair': p['pair']}
        d2b['repulsive'] = p.get('repulsive', False)
        d2b['nlead'] = leading_trim
        d2b['ntrail'] = trailing_trim
        d2b['spacing'] = spacing
        d2b['rc2b'] = p['cutoff']
        res = p['resolution']
        d2b['ncoef'] = res +3  # カットオフで値が０となるように３を加える
        d2b['nknot'] = d2b['ncoef'] +4  # 最短距離にpaddingするため４を加える
        d2b['coefs'] = [ 0.0 for i in range(d2b['ncoef'])]
        dk = (d2b['rc2b']-rmin)/(d2b['nknot']-8+1)
        d2b['knots'] = [ rmin for i in range(d2b['nknot']) ]
        d2b['knots'][-4:] = [ d2b['rc2b'] for i in range(4) ]
        d2b['knots'][3:-4] = [ rmin +dk*i for i in range(res)]
        prms['2B'].append(d2b)
    print(f' Max pair_cutoff = {rmax:.2f}')

    ## 3B
    rmax3 = max(trio_cutoffs)
    #trio_res = config['trio_resolutions']
    cosmin, cosmax = (-1.0, 1.0)
    for it, t in enumerate(trios):
        d3b = {'trio': t['trio']}
        d3b['repulsive'] = t.get('repulsive', False)
        d3b['nlead'] = leading_trim
        d3b['ntrail'] = trailing_trim
        d3b['spacing'] = spacing
        rcij, rcik = t['cutoffs'][0], t['cutoffs'][1]
        d3b['rcij'] = rcij
        d3b['rcik'] = rcik
        resij = t['resolutions'][0]
        resik = t['resolutions'][1]
        rescs = t['resolutions'][2]
        #...r_ij
        d3b['ncfij'] = resij +3
        d3b['nknij'] = d3b['ncfij'] +4
        dr = (rcij - rmin)/(d3b['nknij']-8+1) 
        d3b['knij'] = [ rmin for i in range(d3b['nknij']) ]
        d3b['knij'][-4:] = [ rcij for i in range(4) ]
        d3b['knij'][3:-4] = [ rmin +dr*i for i in range(resij) ]
        d3b['cfij'] = np.zeros(d3b['ncfij'])
        d3b['cfij'][:-3] = 1.0
        #...r_ik
        d3b['ncfik'] = resik +3
        d3b['nknik'] = d3b['ncfik'] +4
        dr = (rcik - rmin)/(d3b['nknik']-8+1)
        d3b['knik'] = [ rmin for i in range(d3b['nknik']) ]
        d3b['knik'][-4:] = [ rcik for i in range(4) ]
        d3b['knik'][3:-4] = [ rmin +dr*i for i in range(resik) ]
        d3b['cfik'] = np.zeros(d3b['ncfik'])
        d3b['cfik'][:-3] = 1.0
        #...cos
        d3b['ncfcs'] = rescs +3
        d3b['nkncs'] = d3b['ncfcs'] +4
        dk = (cosmax - cosmin)/(d3b['nkncs']-6-1) # 両端３つづつはノード点が同じ値
        d3b['kncs'] = [ cosmin for i in range(d3b['nkncs']) ]
        d3b['kncs'][-4:] = [ cosmax for i in range(4) ]
        d3b['kncs'][3:-4] = [ cosmin +dk*i for i in range(rescs) ]
        d3b['cfcs'] = np.array([ init3b for i in range(d3b['ncfcs'])])

        prms['3B'].append(d3b)
    print(f' Max trio_cutoff = {rmax3:.2f}')

    ans = 'y'
    if os.path.exists(outfname):
        ans = input(f' Overwrite {outfname}? [y/N] ').strip().lower()
    if ans != 'y':
        print(' Overwrite canceled.')
    else:
        print(f" --> {outfname}")
        write_params_uf3d(prms, outfname=outfname,
                          author=author, overwrite=True)
    return prms


def prms_to_fp(prms,
               outfname='in.vars.fitpot',
               specorder=[],
               correct=True):
    from nappy.elements import get_number_from_symbol
    
    fpvars = []
    vranges = []

    for d1b in prms['1B']:
        spi = d1b['species']
        erg = d1b['erg']
        fpvars.append(erg)
        vranges.append((-1e+10, 1e+10))

    rc2max = 0.0
    for d2b in prms['2B']:
        pair = d2b['pair']
        si, sj = pair
        #print(pair)
        knots = d2b['knots']
        ncoef = d2b['ncoef']
        coefs = d2b['coefs']
        ntrail = d2b['ntrail']
        rc2max = max(rc2max, d2b['rc2b'])
        if correct:
            zi = get_number_from_symbol(si)
            zj = get_number_from_symbol(sj)
            coefs = correct_wZBL(knots[-1], knots, coefs,
                                 zi=zi, zj=zj)
        vmin = -1e+10
        if d2b.get('repulsive',False) and correct:
            vmin = 0.0
        for i in range(ncoef-ntrail):
            fpvars.append(coefs[i])
            vranges.append((vmin, 1e+10))
        for i in range(ncoef-ntrail, ncoef):
            fpvars.append(coefs[i])
            vranges.append((0.0, 0.0))

    rc3max = 0.0
    for d3b in prms['3B']:
        trio = d3b['trio']
        ntrail = d3b['ntrail']
        #print(trio)
        ncfij = d3b['ncfij']
        ncfik = d3b['ncfik']
        ncfcs = d3b['ncfcs']
        #...rcij, rcik: cutoff for each pair
        rcij, rcik = d3b['rcij'], d3b['rcik']
        rc3max = max(rc3max, rcij, rcik)
        cfij, cfik, cfcs = d3b['cfij'], d3b['cfik'], d3b['cfcs']
        for i in range(ncfij):
            fpvars.append(cfij[i])
            if i < ncfij -ntrail:
                vranges.append((0.0, 1e+10))
            else:
                vranges.append((0.0, 0.0))
        for i in range(ncfik):
            fpvars.append(cfik[i])
            if i < ncfik -ntrail:
                vranges.append((0.0, 1e+10))
            else:
                vranges.append((0.0, 0.0))
        for i in range(ncfcs):
            fpvars.append(cfcs[i])
            vranges.append((0.0, 1e+10))
    write_vars_fitpot(outfname, fpvars, vranges, rc2max, rc3max)
    return None


def create_conditions(prms,
                      outfname='in.vars.conditions'):
    msgline = []
    msgline.append('# var-ID-LHS,  operator,  var-ID-RHS\n')

    inc = 0
    nconds = 0
    for d1 in prms['1B']:
        inc += 1
    for d2 in prms['2B']:
        pair = d2['pair']
        ncoef = d2['ncoef']
        ntrail = d2['ntrail']
        for i in range(ncoef-ntrail):
            inc += 1
            if d2.get('repulsive', False):
                if i != 0:
                    msgline.append(f'   {inc:d}  <  {inc-1:d}\n')
                    nconds += 1
        for i in range(ncoef-ntrail, ncoef):
            inc += 1
    wgt = 1.0
    msgline.insert(1, f'  {nconds:d}  {wgt:.3f}\n')
    print(f' --> {outfname}')
    with open(outfname, 'w') as f:
        for l in msgline:
            f.write(l)
    return None
    

def knot_index(r, knots):
    n = 0
    for i, t in enumerate(knots):
        if r > t:
            n = i
        else:
            return n
    return n


def correct_wZBL(point, knots, coefs, zi=1.0, zj=1.0):
    """
    ZBL potentialに漸近するようにshort-range補正を行う．

    ＊See Goodnote 2025-04-22
    C[n-3] = (2 dt)^3/3! *phi_ZBL'''(r) +3*C[n-2] -3*C[n-1] +C[n]
    phi_ZBL'''(r)のrは，short-rangeを行う各ノード点t[n]とでもする？
    """
    nr = knot_index(point, knots)
    n = nr-3
    dknot = knots[len(knots)//2] -knots[len(knots)//2-1]
    nodes = []
    while True:
        n -= 1
        if n < 0:
            break
        # ZBL contribution
        d3zbl = dddzbl(knots[n+3],zi,zj)
        target = 3.0*(coefs[n+1]-coefs[n+2]) +coefs[n+3] -(2.0*dknot)**3/(3*2)*d3zbl
        coefs[n] = max(coefs[n], target)
        nodes.append(n)
        #print('n,knot =',n,knots[n+4])
    return coefs


def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    verbose = args['--verbose']
    if verbose:
        ic.enable()
        ic("Verbose mode ON.")

    confname = args['CONFIG_UF3D_FILE']
    ic(f"Config file: {confname}")
    config = load_yaml(confname)

    correct = not args['--no-correct']

    specorder = args['--specorder'].split(',')
    if specorder[0] == 'None':
        raise Exception('specorder must be specified.')

    outfprms = 'in.params.uf3d'
    outfpvars = 'in.vars.fitpot'
    outfcond = 'in.vars.conditions'
    author = args['--author']
    init3b = float(args['--init3b'])
    
    prms_uf3d = create_params(config,
                              outfname=outfprms,
                              author=author,
                              init3b=init3b)

    prms_to_fp(prms_uf3d,
               outfname=outfpvars,
               specorder=specorder,
               correct=correct)
    if correct:
        create_conditions(prms_uf3d,
                          outfname=outfcond)


if __name__ == "__main__":

    main()
