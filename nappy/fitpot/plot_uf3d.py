#!/usr/bin/env python
"""
Plot 2B and 3B functions of UF3D potential.

Usage:
  {0:s} [options] UF3D_PARAM_FILE

Options:
  -h, --help          Show this message and exit.
  -v, --verbose       Verbose output. [default: False]
  --postfix POSTFIX   Postfix name added to each output png file name. [default: None]
  --theta             Whether or not plot 3B with theta instead of -cos(theta). [default: False]
"""
import os,sys
from docopt import docopt
import numpy as np
from nappy.fitpot.uf3util import b_spl, get_bspl_curve
from icecream import ic
ic.disable()

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


def plot_2b(fname, label, knots, coefs, rs, xlim=None, ylim=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import copy
    sns.set_theme(context='talk', style='ticks')
    
    fig, ax = plt.subplots(figsize=(5,4))
    for i in range(len(coefs)):
        coefs_tmp = copy.copy(coefs)
        coefs_tmp[:] = 0.0
        coefs_tmp[i] = coefs[i]
        plot_b_spl(ax, rs, knots, coefs_tmp, plot=False, fill=True,
                   alpha=0.2, color='gray')
    plot_b_spl(ax, rs, knots, coefs,
               alpha=1.0, linewidth=3.0, label=label,
               linestyle='-', color='black')
    #ax.set_title(title)
    ax.set_xlabel('Distance (Å)')
    ax.set_ylabel('Energy (eV)')
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(frameon=False, handlelength=0, handletextpad=0)
    #plt.show()
    plt.savefig(fname, dpi=150,
                bbox_inches='tight',
                format='png')
    print(f' --> {fname}')
    return None
    


def plot_uf3d_2b(prms,
                 vline_at=None,
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
        si, sj = p['pair']
        fname = f'graph_2B_{si}-{sj}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        coefs = p['coefs']
        knots = p['knots']
        rs = np.linspace(0.11, max(knots),100)
        ic(pair, coefs)
        label = f'{si}-{sj}'
        plot_2b(fname, label, knots, coefs,rs,xlim,ylim)
    return None


def plot_uf3d_3b(prms,
                 postfix = '',
                 theta = False):
    """
    Plot 3-body terms of UF3D potential.
    If theta == True, plot as funcs of theta instead of (-cos(theta)).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import copy

    sns.set_theme(context='talk', style='ticks')

    tiny = 1.0e-8

    # Analyze 3B
    minyij = -0.1
    maxyij = 1.0
    minxij = 10.0
    maxxij = 0.0
    minycs = -0.1
    maxycs = 1.0
    minxcs = 10.0
    maxxcs = 0.0
    for t in prms['3B']:
        cfij = np.array(t['cfij'])
        knij = np.array(t['knij'])
        cfik = np.array(t['cfik'])
        knik = np.array(t['knik'])
        cfcs = np.array(t['cfcs'])
        kncs = np.array(t['kncs'])
        
        rsij = np.linspace(0.11, max(knij),100)
        esij = get_bspl_curve(rsij,knij,cfij)
        rsik = np.linspace(0.11, max(knik),100)
        esik = get_bspl_curve(rsik,knik,cfik)
        minxij = min(minxij, knij.min(), knik.min())
        maxxij = max(maxxij, knij.max(), knik.max())
        minyij = min(minyij, esij.min(), esik.min())
        maxyij = max(maxyij, esij.max(), esik.max())
        cs = np.linspace(-1.0+tiny, 1.0-tiny, 100)
        ts = np.arccos(cs)
        escs = get_bspl_curve(cs, kncs, cfcs)
        minxcs = min(minxcs, cs.min())
        maxxcs = max(maxxcs, cs.max())
        minycs = min(minycs, escs.min())
        maxycs = max(maxycs, escs.max())
        #ic(trio,miny,maxy)
        
    # Plot
    xlimij = (minxij, maxxij)
    ylimij = (minyij -0.1*(maxyij-minyij),
              maxyij +0.1*(maxyij-minyij))
    xlimcs = (minxcs, maxxcs)
    ylimcs = (minycs -0.1*(maxycs-minycs),
              maxycs +0.1*(maxycs-minycs))
    for t in prms['3B']:
        si, sj, sk = t['trio']
        #...bond ij part
        fname = f'graph_3B_{si}-{sj}_in_{sj}-{si}-{sk}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        cfij = t['cfij']
        knij = t['knij']
        rs = np.linspace(0.11, max(knij), 100)
        #print(xlimij, ylimij, rs)
        label = f'{si}-{sj}\nin {sj}-{si}-{sk}'
        plot_2b(fname, label, knij, cfij, rs, xlimij, ylimij)
        print(f' --> {fname}')
        #...bond ik part
        fname = f'graph_3B_{si}-{sk}_in_{sj}-{si}-{sk}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        cfik = t['cfik']
        knik = t['knik']
        rs = np.linspace(0.11, max(knik), 100)
        label = f'{si}-{sk}\nin {sj}-{si}-{sk}'
        plot_2b(fname, label, knik, cfik, rs, xlimij, ylimij)
        print(f' --> {fname}')
        #...angular part
        fname = f'graph_3B_{sj}-{si}-{sk}'
        if len(postfix) > 0:
            fname += f'_{postfix}'
        fname += '.png'
        label = f'{sj}-{si}-{sk}'
        coefs = np.array(t['cfcs'])
        ic((si,sj,sk), coefs)
        knots = t['kncs']
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
                   alpha=1.0, linewidth=3.0, label=label,
                   linestyle='-', color='black')
        #ax.set_title(title)
        ax.set_xlabel(r'$\cos \theta_{ijk}$')
        ax.set_ylabel('Energy (eV)')
        ax.set_ylim(*ylimcs)
        tick_indices = np.linspace(0, len(cs)-1, 5, dtype=int)
        tick_positions = cs[tick_indices]
        tick_labels = [ 1.0, 0.5, 0, -0.5, -1.0]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.legend(frameon=False, handlelength=0, handletextpad=0)
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
    from nappy.fitpot.uf3util import read_params_uf3d, write_params_uf3d
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    verbose = args['--verbose']
    if verbose:
        ic.enable()
        ic("Verbose mode on.")

    infname = args['UF3D_PARAM_FILE']
    print(" UF3D param file: ", infname)
    postfix = args['--postfix']
    if postfix == 'None': postfix = ''
    ic(postfix)

    # Load param file
    uf3d_prms = read_params_uf3d(infname)

    # Plot 2-bodies
    plot_uf3d_2b(uf3d_prms)

    theta = args['--theta']
    # Plot 3-bodies
    plot_uf3d_3b(uf3d_prms, theta=theta)

if __name__ == "__main__":

    main()
#
