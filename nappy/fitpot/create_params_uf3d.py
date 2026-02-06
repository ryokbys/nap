#!/usr/bin/env python
"""
Create in.params.uf3d file.

Usage:
  {0:s} [options] CONFIG_UF3D_FILE

Options:
  -h, --help            Show this message and exit.
  -v, --verbose         Verbose output. [default: False]
  -o, --output OUTPUT   Output file name. [default: in.params.uf3d]
  --author AUTHOR       Author of the param file. [default: None]
"""
import os, sys
from docopt import docopt
import numpy as np
from nappy.fitpot.uf3util import write_params_uf3d
import yaml
from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "260206"


def load_yaml(file_path):
    import yaml
    with open(file_path, 'r', encoding='utf-8') as file:
        data = yaml.safe_load(file,)
    return data


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

    outfname = args['--output']
    ic(f"Output file: {outfname}")

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
               'epot': 0.0}
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
        #...r_ik
        d3b['ncfik'] = resik +3
        d3b['nknik'] = d3b['ncfik'] +4
        dr = (rcik - rmin)/(d3b['nknik']-8+1)
        d3b['knik'] = [ rmin for i in range(d3b['nknik']) ]
        d3b['knik'][-4:] = [ rcik for i in range(4) ]
        d3b['knik'][3:-4] = [ rmin +dr*i for i in range(resik) ]
        d3b['cfik'] = np.zeros(d3b['ncfik'])
        #...cos
        d3b['ncfcs'] = rescs +3
        d3b['nkncs'] = d3b['ncfcs'] +4
        dk = (cosmax - cosmin)/(d3b['nkncs']-6-1) # 両端３つづつはノード点が同じ値
        d3b['kncs'] = [ cosmin for i in range(d3b['nkncs']) ]
        d3b['kncs'][-4:] = [ cosmax for i in range(4) ]
        d3b['kncs'][3:-4] = [ cosmin +dk*i for i in range(rescs) ]
        d3b['cfcs'] = np.zeros(d3b['ncfcs'])

        prms['3B'].append(d3b)
    print(f' Max trio_cutoff = {rmax3:.2f}')

    ans = 'y'
    if os.path.exists(outfname):
        ans = input(f' Overwrite {outfname}? [y/N] ').strip().lower()
    if ans != 'y':
        print(' Overwrite canceled.')
    else:
        author = args['--author']
        if author == 'None':
            author = __author__
        write_params_uf3d(prms, outfname=outfname,
                          author=author, overwrite=True)
        print(f' --> {outfname}')

if __name__ == "__main__":

    main()
