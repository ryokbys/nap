#!/usr/bin/env python
"""
Plot 2B and 3B functions of UF3L potential.

Usage:
  {0:s} [options] UF3L_PARAM_FILE

Options:
  -h, --help          Show this message and exit.
  -v, --verbose       Verbose output. [default: False]
  --postfix POSTFIX   Postfix name added to each output png file name. [default: None]
  --theta             Whether or not plot 3B with theta instead of -cos(theta). [default: False]
"""
import os,sys
from docopt import docopt
import numpy as np
from icecream import ic
from nappy.fitpot.uf3util import b_spl, get_bspl_curve

__author__ = "RYO KOBAYASHI"
__version__ = "260208"


def plot_b_spl(ax, rs, knots, coefs, plot=True, fill=False, **kwargs):
    es = np.zeros(len(rs))
    for i,r in enumerate(rs):
        nr, b = b_spl(r, knots)
        tmp = 0.0
        for l in range(-3,1):  # -3,-2,-1,0
            n = nr +l
            if n < 0 or n > len(knots)-4: continue
            c2t = coefs[n]
            #if n < 0:
            #    c2t = coefs[0]
            #else:
            #    c2t = coefs[n]
            tmp += c2t *b[l+3]
        es[i] = tmp
    if plot:
        ax.plot(rs, es, **kwargs)
    if fill:
        ax.fill_between(rs, es, **kwargs)
    return None


def plot_uf3l_2b(prms,
                 postfix=''):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import copy

    sns.set_theme(context='talk', style='ticks')
    # Analyze 2b
    miny = -1.0
    maxy =  5.0
    minx = 10.0
    maxx = 0.0
    for p in prms['2B']:
        pair = p['pair']
        coefs = p['coefs']
        knots = p['knots']
        rs = np.linspace(0.11, max(knots),100)
        minx = min(rs.min(), minx)
        maxx = max(rs.max(), maxx)
        es = get_bspl_curve(rs, knots, coefs)
        min1st = -1
        max1st = -1
        for i,ei in enumerate(es):
            if i == 0: continue
            if ei > es[i-1] and min1st < 0:
                min1st = i-1
            if ei < es[i-1] and min1st >= 0 and max1st < 0:
                max1st = i-1
        miny = min(miny, es[min1st])
        maxy = max(maxy, es[max1st])
        ic(pair, min1st, miny)
        ic(pair, max1st, maxy)

    # Plot
    xlim = (minx, maxx)
    ylim = (miny*(1.0+0.1), maxy*(1.0+0.1))
    ic(xlim)
    ic(ylim)
    for p in prms['2B']:
        pair = p['pair']
        si, sj = pair
        fname = f'graph_2b_{si}-{sj}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        coefs = p['coefs']
        knots = p['knots']
        rs = np.linspace(0.11, max(knots),100)
        ic(pair, coefs)
        fig, ax = plt.subplots(figsize=(5,4))
        for i in range(len(coefs)):
            coefs_tmp = copy.copy(coefs)
            coefs_tmp[:] = 0.0
            coefs_tmp[i] = coefs[i]
            plot_b_spl(ax, rs, knots, coefs_tmp, plot=False, fill=True,
                       alpha=0.2, color='gray')
        plot_b_spl(ax, rs, knots, coefs,
                   alpha=1.0, label=f'{si}-{sj}',linewidth=3.0,
                   linestyle='-', color='black')
        ax.set_xlabel('Distance (Å)')
        ax.set_ylabel('Energy (eV)')
        ax.set_ylim(*ylim)
        ax.set_xlim(*xlim)
        ax.legend(frameon=False)
        #plt.show()
        plt.savefig(fname, dpi=150,
                    bbox_inches='tight',
                    format='png')
        print(f' --> {fname}')
    return None


def plot_uf3l_3b(prms,
                 postfix = '',
                 theta = False):
    """
    Plot 3-body terms of UF3L potential.
    If theta == True, plot as funcs of theta instead of (-cos(theta)).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import copy

    sns.set_theme(context='talk', style='ticks')

    tiny = 1.0e-8

    # Analyze 3B
    miny = -1.0
    maxy = 1.0
    minx = 10.0
    maxx = 0.0
    for t in prms['3B']:
        trio = t['trio']
        coefs = np.array(t['coefs'])
        knots = t['knots']
        cs = np.linspace(-1.0+tiny, 1.0-tiny, 100)
        ts = np.arccos(cs)
        es = get_bspl_curve(cs, knots, coefs)
        miny = min(miny, es.min())
        maxy = max(maxy, es.max())
        ic(trio,miny,maxy)
        minx = min(minx, cs.min())
        maxx = max(maxx, cs.max())

    # Plot
    xlim = (minx, maxx)
    ylim = (miny*(1.0-0.1), maxy*(1.0+0.1))
    for t in prms['3B']:
        trio = t['trio']
        si, sj, sk = trio
        fname = f'graph_3b_{si}-{sj}-{sk}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        coefs = np.array(t['coefs'])
        ic(trio, coefs)
        knots = t['knots']
        cs = np.linspace(-1.0+tiny, 1.0-tiny, 100)
        ts = np.arccos(cs)
        fig, ax = plt.subplots(figsize=(5,4))
        for i in range(len(coefs)):
            coefs_tmp = copy.copy(coefs)
            coefs_tmp[:] = 0.0
            coefs_tmp[i] = coefs[i]
            plot_b_spl(ax, cs, knots, coefs_tmp,
                       alpha=0.2, color='gray', fill=True)
        plot_b_spl(ax, cs, knots, coefs,
                   alpha=1.0, label=f'{si}-{sj}-{sk}',linewidth=3.0,
                   linestyle='-', color='black')
        ax.set_xlabel(r'$\cos \theta_{ijk}$')
        ax.set_ylabel('Energy (eV)')
        ax.set_ylim(*ylim)
        tick_indices = np.linspace(0, len(cs)-1, 5, dtype=int)
        tick_positions = cs[tick_indices]
        tick_labels = [ -1.0, -0.5, 0, 0.5, 1.0]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.legend(frameon=False)
        if theta:
            ax_top = ax.secondary_xaxis('top',)
            ax_top.set_xlabel(r'$\theta_{ijk}$ (rad)')
            tick_indices = np.linspace(0, len(cs)-1, 3, dtype=int)
            tick_positions = cs[tick_indices]
            tick_labels = [  '0', r'$\pi/2$', r'$\pi$',]
            ax_top.set_xticks(tick_positions)
            ax_top.set_xticklabels(tick_labels)
        #plt.show()
        plt.savefig(fname, dpi=150,
                    bbox_inches='tight',
                    format='png')
        print(f' --> {fname}')
    return


def main():
    from nappy.fitpot.uf3util import read_params_uf3l, write_params_uf3l
    ic.disable()
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    verbose = args['--verbose']
    if verbose:
        ic.enable()
        ic("Verbose mode on.")

    infname = args['UF3L_PARAM_FILE']
    print(" UF3L param file: ", infname)
    postfix = args['--postfix']
    if postfix == 'None': postfix = ''
    ic(postfix)

    # Load param file
    uf3l_prms = read_params_uf3l(infname)

    # Plot 2-bodies
    plot_uf3l_2b(uf3l_prms)

    theta = args['--theta']
    # Plot 3-bodies
    plot_uf3l_3b(uf3l_prms, theta=theta)

if __name__ == "__main__":

    main()
#
