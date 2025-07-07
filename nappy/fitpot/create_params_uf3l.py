#!/usr/bin/env python
"""
Create in.params.uf3l file.

Usage:
  {0:s} [options] CONFIG_UF3L_FILE

Options:
  -h, --help            Show this message and exit.
  -v, --verbose         Verbose output. [default: False]
  -o, --output OUTPUT   Output file name. [default: in.params.uf3l]
  --author AUTHOR       Author of the param file. [default: None]
"""
import os, sys
from docopt import docopt
import numpy as np
from nappy.fitpot.uf3util import write_params_uf3l
import yaml
from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "250705"


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
        ic("Verbose mode on.")

    confname = args['CONFIG_UF3L_FILE']
    ic(f"Config file: {confname}")
    config = load_yaml(confname)

    outfname = args['--output']
    ic(f"Output file: {outfname}")

    pairs = config['pairs']
    pair_cutoffs = config['pair_cutoffs']
    trios = config['trios']
    trio_cutoffs = config['trio_cutoffs']

    prms = {}
    prms['1B'] = {}
    prms['2B'] = {}
    prms['3B'] = {}

    ## 1B
    species = config['species']
    for spi in species:
        prms['1B'][spi] = 0.0

    ## Params common to 2B and 3B
    leading_trim = 0
    trailing_trim = 3
    spacing = 'nk' # non-uniform knots

    ## 2B
    rmin = 0.1
    rmax = max(pair_cutoffs)
    pair_res = config['pair_resolutions']
    #resolution = 20  # coef数はこれ＋３， knot数はこれ＋７
    for ip, p in enumerate(pairs):
        d2b = {}
        d2b['nlead'] = leading_trim
        d2b['ntrail'] = trailing_trim
        d2b['spacing'] = spacing
        d2b['rc2b'] = pair_cutoffs[ip]
        res = pair_res[ip]
        d2b['ncoef'] = res +3  # カットオフで値が０となるように３を加える
        d2b['nknot'] = d2b['ncoef'] +4  # 最短距離にpaddingするため４を加える
        d2b['coefs'] = [ 0.0 for i in range(d2b['ncoef'])]
        dk = (d2b['rc2b']-rmin)/(d2b['nknot']-8+1)
        d2b['knots'] = [ rmin for i in range(d2b['nknot']) ]
        d2b['knots'][-4:] = [ d2b['rc2b'] for i in range(4) ]
        d2b['knots'][3:-4] = [ rmin +dk*i for i in range(res)]
        prms['2B'][tuple(p)] = d2b
    print(f' Max pair_cutoff = {rmax:.2f}')

    ## 3B
    rmax3 = max(trio_cutoffs)
    trio_res = config['trio_resolutions']
    cosmin, cosmax = (-1.0, 1.0)
    for it, t in enumerate(trios):
        d3b = {}
        d3b['nlead'] = leading_trim
        d3b['ntrail'] = trailing_trim
        d3b['spacing'] = spacing
        d3b['rc'] = trio_cutoffs[it]
        d3b['gmj'] = 1.0
        d3b['gmk'] = 1.0
        res = trio_res[it]
        d3b['ncoef'] = res +3
        d3b['nknot'] = d3b['ncoef'] +4
        dk = (cosmax - cosmin)/(d3b['nknot']-6-1) # 両端３つづつはノード点が同じ値
        d3b['knots'] = [ cosmin for i in range(d3b['nknot']) ]
        d3b['knots'][-4:] = [ cosmax for i in range(4) ]
        d3b['knots'][3:-4] = [ cosmin +dk*i for i in range(res) ]
        d3b['coefs'] = np.zeros(d3b['ncoef'])
        prms['3B'][tuple(t)] = d3b
    print(f' Max trio_cutoff = {rmax3:.2f}')

    if os.path.exists(outfname):
        ans = input(f' Overwrite {outfname}? [y/N] ').strip().lower()
        if ans != 'y':
            print(' Overwrite canceled.')
        else:
            author = args['--author']
            if author == 'None':
                author = __author__
            write_params_uf3l(prms, outfname=outfname,
                              author=author, overwrite=True)
            print(f' --> {outfname}')

if __name__ == "__main__":

    main()
