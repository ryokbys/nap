#!/usr/bin/env python
"""
Correct 2- and 3-body terms of UF3L potential.

Usage:
  {0:s} [options] FILENAME

Options:
  -h, --help      Show this message and exit.
  -v, --verbose   Verbose output. [default: False]
  -c, --config CONFIG
                  Config yaml file name. [default: config_uf3l.yaml]
  --adf-path PATH
                  Path to the out.adf. [default: out.adf]
"""
import os, sys
from docopt import docopt
import numpy as np
from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "250705"

azbl = (0.1818, 0.5099, 0.2802, 0.02817)
bzbl = (3.2, 0.9423, 0.4029, 0.2016)

def load_yaml(file_path):
    import yaml
    with open(file_path, 'r', encoding='utf-8') as file:
        data = yaml.safe_load(file,)
    return data


def knot_index(r, knots):
    n = 0
    for i, t in enumerate(knots):
        if r > t:
            n = i
        else:
            return n
    return n


def correct_coefs_wZBL(point, knots, coefs, zi=1.0, zj=1.0):
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


def dzbl(r, zi, zj):
    k = 14.4  # 1/(4*pi*eps0) in eV*Ang/e^2 unit
    rs = 0.4683766 / ( zi**0.23 + zj**0.23 )
    return (k*zi*zj) * (-1.0/r**2 *phi(r/rs) +1.0/r/rs *dphi(r/rs))


def dddzbl(r, zi, zj):
    k = 14.4  # 1/(4*pi*eps0) in eV*Ang/e^2 unit
    rs = 0.4683766 / ( zi**0.23 + zj**0.23 )
    return (k*zi*zj) *( -6.0/r**4 *phi(r/rs) +2.0/r**3/rs*dphi(r/rs) \
                       +2.0/r**3/rs*dphi(r/rs) -1.0/r**2/rs**2*ddphi(r/rs) \
                       +2.0/r**3/rs*dphi(r/rs) -1.0/r**2/rs**2*ddphi(r/rs) \
                       -1.0/r**2/rs**2*ddphi(r/rs) +1.0/r/rs**3*dddphi(r/rs) )


def phi(x):
    #return 0.1818*np.exp(-3.2*x) +0.5099*np.exp(-0.9423*x) \
    #    +0.2802*np.exp(-0.4029*x) +0.02817*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += azbl[i]*np.exp(-bzbl[i]*x)
    return s


def dphi(x):
    #return -0.58176*np.exp(-3.2*x) -0.48047877*np.exp(-0.9423*x) \
    #    -0.11289258*np.exp(-0.4029*x) -0.005679072*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += -azbl[i] * bzbl[i] * np.exp(-bzbl[i]*x)
    return s


def ddphi(x):
    #return 1.861632*np.exp(-3.2*x) +0.452755144971*np.exp(-0.9423*x) \
    #    +0.045484420482*np.exp(-0.4029*x) +0.0011449009152*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += azbl[i] * bzbl[i]**2 * np.exp(-bzbl[i]*x)
    return s


def dddphi(x):
    #return -5.9572224*np.exp(-3.2*x) -0.426631173106*np.exp(-0.9423*x) \
    #    -0.0183256730122*np.exp(-0.4029*x) -0.000230812024504*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += -azbl[i] * bzbl[i]**3 * np.exp(-bzbl[i]*x)
    return s


def get_comb_index(comb, comb_list):
    """
    Get index of the given COMB in the COMB_LIST.
    If there is not, return -1.
    """
    if len(comb) != len(comb_list[0]):
        raise TypeError('comb and comb_list is not consistent')
    lenc = len(comb)
    for idx,ci in enumerate(comb_list):
        identical = True
        for i in range(lenc):
            if comb[i] != ci[i]:
                identical = False
        if identical:
            return idx
    return -1


def correct2b(uf3l_prms, config):
    from nappy.elements import get_number_from_symbol
    pairs = config['pairs']
    pair_rins = config['pair_rins']

    uf32b = uf3l_prms['2B']

    for pair in uf32b.keys():
        si, sj = pair
        ipair = get_comb_index(pair, pairs)
        zi = get_number_from_symbol(si)
        zj = get_number_from_symbol(sj)
        coefs = uf32b[pair]['coefs']
        knots = uf32b[pair]['knots']
        new_coefs = correct_coefs_wZBL(pair_rins[ipair],
                                       knots, coefs,
                                       zi=zi, zj=zj)
        uf3l_prms['2B'][pair]['coefs'] = new_coefs

    return uf3l_prms


def correct3b(uf3l_prms, config, adf_file_path,
              sgm=15, ):
    """
    Correct 3B so that the V(-cos) becomes an inverse of ADF.
    """
    from nappy.adf import read_adf
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(context='talk', style='ticks')
    uf33b = uf3l_prms['3B']
    trios = config['trios']
    trio_vmins = config['trio_vmins']
    trio_vmaxs = config['trio_vmaxs']

    #...Create V(theta) from ADF data in out.adf
    assert os.path.exists(adf_file_path), f' {adf_file_path} does not exists!'
    ts, adfs = read_adf(adf_file_path)
    for trio in uf33b.keys():
        assert trio in adfs.keys(), f'There is no such a trio, {trio}, in {adf_file_path}.'
        adf = adfs[trio]
        coefs = uf33b[trio]['coefs']
        knots = uf33b[trio]['knots']
        bounds = np.array([ (0.0, 1e+10) for c in coefs ])
        itrio = get_comb_index(trio, trios)
        v3min = trio_vmins[itrio]
        v3max = trio_vmaxs[itrio]
        cs, vtgt = gen_vtgt_from_adf(ts, adf, sgm=sgm,
                                     vmax=v3max, vmin=v3min)
        # popt, ier = leastsq(chi, coefs, args=(cs, vtgt, knots),
        #                    full_output=False)
        result = minimize(objective, x0=coefs, args=(cs, vtgt, knots),
                          bounds=bounds)
        uf3l_prms['3B'][trio]['coefs'] = result.x

        fig, ax = plt.subplots(figsize=(5,4))
        ax.plot(cs, vtgt, 'b-', label=f'{trio}')
        ax.set_ylabel('V(-cos)')
        ax.set_xlabel('-cos')
        ax.legend(frameon=False)
        fname = f'graph_vtgt_{itrio}.png'
        plt.savefig(fname, bbox_inches='tight', format='png')
        print(f' --> {fname}')

    return uf3l_prms


def bspl_at(r, coefs, knots,):
    from nappy.fitpot.uf3util import b_spl
    er = 0.0
    nr, b = b_spl(r, knots)
    for l in range(-3, 1):
        n = nr + l
        if n < 0 or n >= len(coefs):
            continue
        c = coefs[n]
        er += c * b[l+3]
    return er


def objective(coefs, xs, ys, knots):
    obj = 0.0
    for i, xi in enumerate(xs):
        yi = ys[i]
        obj += (yi - bspl_at(xi, coefs, knots))**2
    return obj


def gen_vtgt_from_rdf(rs, rdf, sgm=0,):
    """
    UNDER CONSTSTRUCTION!

    Genenrate target V(r) from RDF, p(r).
    """
    from nappy.gaussian_smear import gsmear
    import copy
    eps = 1e-15
    grdf = copy.copy(rdf)
    if sgm > 0:
        grdf = gsmear(rs, rdf, sgm)
    irmax = grdf.argmax()
    for i in range(irmax, 0, -1):
        pass
    #cs = np.array([ -np.cos(t / 180 * np.pi) for t in ts ])
    return rs, vtgt


def gen_vtgt_from_adf(ts, adf, sgm=15, vmax=10.0, vmin=0.0):
    """
    Genenrate target V(-cos) from ADF, p(theta).
    """
    from nappy.gaussian_smear import gsmear
    eps = 1e-15
    gadf = gsmear(ts, adf, sgm)
    gadf /= gadf.max()
    #vtgt = np.array([ -np.log(a + eps) for a in gadf ])
    #adf /= adf.max()
    vtgt = -(gadf - 1.0) * (vmax - vmin) + vmin
    cs = np.array([ -np.cos(t / 180 * np.pi) for t in ts ])
    return cs, vtgt


def main():
    from nappy.fitpot.uf3util import read_params_uf3l, write_params_uf3l
    import copy
    from datetime import datetime

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    verbose = args['--verbose']
    if verbose:
        ic.enable()
        ic(' Verbose mode on.')

    infname = args['FILENAME']
    ic(' read_params_uf3l...')
    uf3l_prms = read_params_uf3l(infname)
    uf3l_prms_tmp = copy.deepcopy(uf3l_prms)

    confname = args['--config']
    ic(f' load_yaml({confname})...')
    config = load_yaml(confname)

    adf_path = args['--adf-path']

    ic(' correct2b...')
    uf3l_prms_tmp = correct2b(uf3l_prms_tmp, config)
    ic(' correct3b...')
    uf3l_prms_tmp = correct3b(uf3l_prms_tmp, config, adf_path,)
    today = datetime.now().strftime('%y%m%d')
    outfname = f'{infname}_corr_{today}'
    write_params_uf3l(uf3l_prms_tmp, outfname=outfname, overwrite=True)
    print(f" --> {outfname}")
    return None


if __name__ == "__main__":

    main()
