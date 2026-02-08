#!/usr/bin/env python
"""
Create in.params.uf3l file.

Usage:
  {0:s} [options] CONFIG_UF3L_FILE

Options:
  -h, --help            Show this message and exit.
  -v, --verbose         Verbose output. [default: False]
  --author AUTHOR  Author of the param file. [default: None]
  --no-correct     Apply prior-knowledge correction. [default: False]
  --specorder SPECORDER
                   Species order in comma-seperated format, e.g.) Li,P,O. [default: None]
  --init3b INIT3B  Initial value of 3-body potential coefficients. [default: 1.0]
"""
import os, sys
from docopt import docopt
import numpy as np
from nappy.fitpot.uf3util import write_params_uf3l
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
                  outfname='in.params.uf3l',
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
               'erg': 0.0 }
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
        d3b['repulsiv'] = t.get('repulsive', False)
        d3b['nlead'] = leading_trim
        d3b['ntrail'] = trailing_trim
        d3b['spacing'] = spacing
        d3b['rcij'] = t['cutoffs'][0]
        d3b['rcik'] = t['cutoffs'][1]
        d3b['gmj'] = 1.0
        d3b['gmk'] = 1.0
        res = t['resolution']
        d3b['ncoef'] = res +3
        d3b['nknot'] = d3b['ncoef'] +4
        dk = (cosmax - cosmin)/(d3b['nknot']-6-1) # 両端３つづつはノード点が同じ値
        d3b['knots'] = [ cosmin for i in range(d3b['nknot']) ]
        d3b['knots'][-4:] = [ cosmax for i in range(4) ]
        d3b['knots'][3:-4] = [ cosmin +dk*i for i in range(res) ]
        d3b['coefs'] = np.array([ init3b for i in range(d3b['ncoef'])])
        prms['3B'].append(d3b)
    print(f' Max trio_cutoff = {rmax3:.2f}')

    ans = 'y'
    if os.path.exists(outfname):
        ans = input(f' Overwrite {outfname}? [y/N] ').strip().lower()
    if ans != 'y':
        print(' Overwrite canceled.')
    else:
        print(f' --> {outfname}')
        write_params_uf3l(prms, outfname=outfname,
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
        #...rcij, rcik: cutoff for each pair
        rcij, rcik = d3b['rcij'], d3b['rcik']
        rc3max = max(rc3max, rcij, rcik)
        ncoef = d3b['ncoef']
        coefs = d3b['coefs']
        # rcij, rcik are cutoff parameters for each pair
        fpvars.append(d3b['rcij'])
        vranges.append((rc3max, rc3max))
        fpvars.append(d3b['rcik'])
        vranges.append((rc3max, rc3max))
        # gmj, gmk should be greater than 0.0
        # gmj, gmk are parameters in exp.
        # Too small values cause violation of energy conservation,
        # and too large values make them negligibly small.
        # But for now, we fix it to 1.0.
        fpvars.append(d3b['gmj'])
        vranges.append((1.0, 1.0)) # gmj, gmk should be greater than 0.0
        fpvars.append(d3b['gmk'])
        vranges.append((1.0, 1.0))
        for i in range(ncoef):
            fpvars.append(coefs[i])
            vranges.append((0.0, 1e+10))
        
    write_vars_fitpot(outfname, fpvars, vranges, rc2max, rc3max)
    return None


def create_conditions(prms,
                      outfname='in.vars.conditions'):
    """
    Write in.vars.conditions.
    """
    msgline = []
    msgline.append('# var-ID-LHS,  operator,  var-ID-RHS\n')
    
    inc = 0
    nconds = 0
    for d1b in prms['1B']:
        inc += 1
    for d2b in prms['2B']:
        ncoef = d2b['ncoef']
        ntrail = d2b['ntrail']
        for i in range(ncoef-ntrail):
            inc += 1
            if d2b.get('repulsive', False):
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
        ic("Verbose mode on.")

    confname = args['CONFIG_UF3L_FILE']
    ic(f"Config file: {confname}")
    config = load_yaml(confname)

    correct = not args['--no-correct']

    specorder = args['--specorder'].split(',')
    if specorder[0] == 'None':
        raise Exception('specorder must be specified.')

    outfprms = 'in.params.uf3l'
    outfpvars = 'in.vars.fitpot'
    outfcond = 'in.vars.conditions'
    author = args['--author']
    init3b = float(args['--init3b'])
    
    prms_uf3l = create_params(config,
                              outfname=outfprms,
                              author=author,
                              init3b=init3b)
    prms_to_fp(prms_uf3l,
               outfname=outfpvars,
               specorder=specorder,
               correct=correct)
    if correct:
        create_conditions(prms_uf3l,
                          outfname=outfcond)
    return None


if __name__ == "__main__":

    main()
