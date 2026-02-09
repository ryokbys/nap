#!/usr/bin/env python
"""
Make parity plots of energy, force, and stress using fitpot outputs.

Usage:
  {0:s} [options] 

Options:
  -h, --help     Show this message and exit.
  -v, --verbose  Verbose output. [default: False]
  --E-range ERANGE
                 Specify the range of energy partiy plot,
                 by comma-separated two values, e.g., -3.0,3.0. [default: None]
  --F-range FRANGE
                 Specify the range of force partiy plot,
                 by comma-separated two values, e.g., -3.0,3.0. [default: None]
  --S-range SRANGE
                 Specify the range of stress partiy plot,
                 by comma-separated two values, e.g., -3.0,3.0. [default: None]
  --refname REFNAME
                 Specify the name of the reference data. [default: DFT]
  --max-plot-num NUM
                 Max num of plotting points. If the data is over it,
                 use KDE plot instead of plotting all the data. [default: 5000]
"""
import os,sys
from docopt import docopt
import numpy as np

from icecream import ic
ic.disable()

__author__ = "RYO KOBAYASHI"
__version__ = "251020"


def read_data(fname):
    import numpy as np
    with open(fname,'r') as f:
        lines = f.readlines()
    data = []
    for i,line in enumerate(lines):
        if line[0] == '#':
            continue
        tmp = line.split()
        data.append([ float(tmp[0]), 
                      float(tmp[1]),
                      float(tmp[3]) ])
    return np.array(data)


def read_trn_tst_data(target):
    if target not in ('erg', 'frc', 'strs'):
        raise ValueError('target must be either erg, frc, or strs.')
    trnfile = f'out.{target}.trn.fin'
    tstfile = f'out.{target}.tst.fin'
    trndata = read_data(trnfile)
    tstdata = read_data(tstfile)
    return trndata, tstdata


def calc_stats(reference, predicted):
    if reference.shape != predicted.shape:
        raise ValueError('Inconsistent data!')
    rmse = np.sqrt(np.mean((reference - predicted)**2))
    mean = np.mean(reference)
    res_sum = np.sum((reference - predicted)**2)
    sq_sum =  np.sum((reference - mean)**2)
    r2 = 1.0 - (res_sum / sq_sum)
    return rmse, r2


def plot(target, trn_data, tst_data, limit_ratio=1.0, xylim=None,
         loc='best', bbox_to_anchor=None, outfname="graph_parity.png",
         refname = "DFT", max_plot_num=5000):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(context='talk', style='ticks')
    
    from scipy.stats import gaussian_kde
    if target == 'erg':
        label = 'energy (eV/atom)'
    elif target == 'frc':
        label = r'force (eV/$\mathrm{\AA}$)'
    elif target == 'strs':
        label = 'stress (GPa)'
    
    fig, ax = plt.subplots(figsize=(5,5), )
    if xylim:
        limits = np.array(xylim)
    else:
        minval = min(trn_data[:,0].min(), trn_data[:,1].min(), 
                     tst_data[:,0].min(), tst_data[:,1].min())
        maxval = max(trn_data[:,0].max(), trn_data[:,1].max(), 
                     tst_data[:,0].max(), tst_data[:,1].max())
        padding = ( maxval - minval ) * 0.05
        limits = np.array([ minval - padding, maxval + padding ]) * limit_ratio
    ax.plot(limits, limits, linestyle="--", color="k", linewidth=2, 
            label="", zorder=1)
    tst_rmse, tst_r2 = calc_stats(tst_data[:,0], tst_data[:,1])
    
    if len(tst_data) < max_plot_num:
        # データ数が少ない場合のみ，trainingとtestの両方をフルデータで表示．
        ax.plot(trn_data[:,0], trn_data[:,1], 
                'ro', mec='k', ms=6, label='training', zorder=2)
        ax.plot(tst_data[:,0], tst_data[:,1], 
                'bs', mec='k', ms=6, label='test', zorder=3)
        margin = 0.05
        ax.text(1-margin, 0+margin, f'RMSE = {tst_rmse:0.3f}\nR^2 = {tst_r2:0.3f}', 
                transform=ax.transAxes, ha='right', va='bottom')
        ax.legend(frameon=False, handletextpad=0, 
                  loc=loc, bbox_to_anchor=bbox_to_anchor)

    else:
        # データ数が多い場合は，testデータだけをKDEを使って表示．
        sample_size = max_plot_num  # サンプルする点の数
        indices = np.random.choice(len(tst_data), sample_size, replace=False)  # ランダムサンプリング
        x_sampled = tst_data[indices,0]
        y_sampled = tst_data[indices,1]

        # ② KDE（カーネル密度推定）を計算
        xy = np.vstack([x_sampled, y_sampled])
        kde = gaussian_kde(xy)(xy)  # 各点の密度を計算
        
        # ③ KDE の値に基づいて散布図をプロット（色付け）
        ax.scatter(x_sampled, y_sampled, c=kde, 
                   cmap='jet', s=10, edgecolors='none', zorder=2)
        margin = 0.05
        ax.text(1-margin, 0+margin, 
                f'RMSE = {tst_rmse:0.3f}\nR^2 = {tst_r2:0.3f}', 
                transform=ax.transAxes, ha='right', va='bottom')
    
    ax.set_xlabel(f'{refname} {label}')
    ax.set_ylabel(f'FF {label}')
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal")
    plt.savefig(outfname, format='png', dpi=300, bbox_inches='tight')
    print(f' --> {outfname}')
    return None


def main():

    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    verbose = args['--verbose']
    if verbose:
        ic.enable()
        ic("Verbose mode ON.")

    maxnum = int(args['--max-plot-num'])

    erg_trn, erg_tst = read_trn_tst_data('erg')
    frc_trn, frc_tst = read_trn_tst_data('frc')
    strs_trn, strs_tst = read_trn_tst_data('strs')

    erange = args['--E-range']
    if erange == 'None':
        erange = None
    else:
        erange = [ float(x) for x in erange.split(',')]
        assert len(erange) == 2, "--E-range should be comma-separated two values."
    frange = args['--F-range']
    if frange == 'None':
        frange = None
    else:
        frange = [ float(x) for x in frange.split(',')]
        assert len(frange) == 2, "--F-range should be comma-separated two values."
    srange = args['--S-range']
    if srange == 'None':
        srange = None
    else:
        srange = [ float(x) for x in srange.split(',')]
        assert len(Srange) == 2, "--S-range should be comma-separated two values."

    refname = args['--refname']
    
    plot('erg', erg_trn, erg_tst, loc='best', xylim=erange,
         outfname='graph_parity_E.png', refname=refname, max_plot_num=maxnum)
    plot('frc', frc_trn, frc_tst, xylim=frange,
         outfname='graph_parity_F.png', refname=refname, max_plot_num=maxnum)
    plot('strs', strs_trn, strs_tst, loc='center right', xylim=srange,
         outfname='graph_parity_S.png', refname=refname, max_plot_num=maxnum)
    return None

        
if __name__ == "__main__":

    main()
#
