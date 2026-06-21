#!/usr/bin/env python
"""
Plot MSD (mean square displacement) vs. time from out.msd file(s)
produced by nappy/msd.py.

Behaviour:
  - Multiple files : compute mean ± SEM across files; draw shaded error band.
  - Single file with error columns (produced with --err) : use those directly.
  - Single file without error columns : plot mean only, no band.

Usage:
  {0:s} [options] MSD_FILE [MSD_FILE...]

Options:
  -h, --help              Show this help message and exit.
  --specorder=SPECORDER   Species in order, comma-separated (e.g. Li,Ge,P,S).
                          Sets subplot order. Inferred from first file if
                          omitted. [default: None]
  -o PREFIX               Output image filename prefix. [default: graph_msd]
  --dpi=DPI               Image resolution in DPI. [default: 150]
  --fit                   Overlay linear fit line and annotate each subplot
                          with the diffusion coefficient D (cm²/s).
                          [default: False]
  --fit-start=FSTART      Fraction of time steps to skip at the start before
                          fitting (to exclude the ballistic regime).
                          [default: 0.2]
"""

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from docopt import docopt
from scipy import stats

__author__ = "nappy-msd skill"
__version__ = "260620b"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_msd_file(fname):
    """Parse an out.msd file (non-xyz, non-com mode) from nappy/msd.py.

    Returns
    -------
    times   : ndarray (nt,)          time in fs
    msds    : dict {spc: ndarray}    total MSD per species in Å²
    errs    : dict {spc: ndarray or None}   SEM per species, or None
    species : list[str]              species in file column order
    """
    species = []
    col_types = {}   # 0-based col index -> ('msd'|'err', spc_str)

    with open(fname) as f:
        raw = f.readlines()

    # The column-header line starts with '#' and contains 'data_ID'
    for line in raw:
        if line.startswith('#') and 'data_ID' in line:
            for m in re.finditer(r'(\d+):msd_(\w+)', line):
                idx = int(m.group(1)) - 1   # 0-based
                spc = m.group(2)
                if spc not in species:
                    species.append(spc)
                col_types[idx] = ('msd', spc)
            for m in re.finditer(r'(\d+):err_(\w+)', line):
                idx = int(m.group(1)) - 1
                spc = m.group(2)
                col_types[idx] = ('err', spc)
            break

    if not species:
        raise ValueError(
            f'No MSD columns found in {fname}. '
            'Make sure msd.py was run without --xyz.'
        )

    # Read numeric data rows
    rows = []
    for line in raw:
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        rows.append([float(v) for v in s.split()])

    if not rows:
        raise ValueError(f'No data rows in {fname}')

    data = np.array(rows)
    times = data[:, 1]   # column 1 = time(fs)

    msds = {}
    errs = {spc: None for spc in species}
    for idx, (typ, spc) in col_types.items():
        if typ == 'msd':
            msds[spc] = data[:, idx]
        else:
            errs[spc] = data[:, idx]

    return times, msds, errs, species


# ---------------------------------------------------------------------------
# Statistics across multiple files
# ---------------------------------------------------------------------------

def multi_file_stats(parsed, all_species):
    """Align datasets to shortest length; return mean ± SEM per species.

    Parameters
    ----------
    parsed      : list of (times, msds, errs, species)
    all_species : ordered list of all species across files

    Returns
    -------
    times     : ndarray (nt,)
    mean_msds : dict {spc: ndarray}
    sem_msds  : dict {spc: ndarray}
    """
    nt = min(len(p[0]) for p in parsed)
    times = parsed[0][0][:nt]

    mean_msds = {}
    sem_msds = {}
    for spc in all_species:
        arrays = [p[1][spc][:nt] for p in parsed if spc in p[1]]
        stacked = np.array(arrays)          # (n_files, nt)
        mean_msds[spc] = stacked.mean(axis=0)
        sem_msds[spc] = stats.sem(stacked, axis=0)

    return times, mean_msds, sem_msds


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plot(times_fs, mean_msds, sem_msds, species, prefix, dpi,
              fit=False, fit_start=0.2):
    """Save MSD vs time (ps) with shaded error band to {prefix}.png.

    When fit=True, overlays a linear fit line on each subplot and annotates
    the title with the diffusion coefficient D computed from the slope:
      D [cm²/s] = slope [Å²/fs] × 1e-1 / (2 × 3)
    The first fit_start fraction of time steps is excluded from the fit to
    avoid the ballistic regime.
    """
    sns.set_theme(context='talk', style='ticks')
    palette = sns.color_palette('tab10', n_colors=max(len(species), 1))

    t = times_fs / 1000.0   # fs -> ps
    n = len(species)

    if n == 1:
        fig, axes = plt.subplots(1, 1, figsize=(5, 5), dpi=dpi)
        axes = [axes]
    else:
        ncols = min(n, 3)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(5 * ncols, 5 * nrows), dpi=dpi)
        axes = np.array(axes).flatten().tolist()

    # Conversion factor: Å²/fs → cm²/s
    _fac = 1.0e-16 / 1.0e-15   # = 0.1
    _dim = 3

    for i, spc in enumerate(species):
        ax = axes[i]
        color = palette[i]
        mean = mean_msds[spc]
        sem  = sem_msds[spc]

        ax.plot(t, mean, color=color, lw=1.8)
        if np.any(sem > 0):
            ax.fill_between(t, mean - sem, mean + sem,
                            color=color, alpha=0.25)

        title = spc
        if fit and len(times_fs) > 2:
            n_skip = max(1, int(len(times_fs) * fit_start))
            t_fit  = times_fs[n_skip:]          # fs, for regression
            msd_fit = mean[n_skip:]
            slope, intercept = np.polyfit(t_fit, msd_fit, 1)   # Å²/fs
            D = slope * _fac / (2.0 * _dim)    # cm²/s
            fit_line = slope * times_fs + intercept             # Å²
            ax.plot(t, fit_line, '--', color=color, lw=1.2, alpha=0.7)
            title = f'{spc}   D = {D:.3e} cm²/s'
            print(f' {spc}: D = {D:.4e} cm²/s  (fit from {fit_start*100:.0f}% of data)')

        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('MSD (Å²)')
        ax.set_title(title)
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(bottom=0)
        ax.set_box_aspect(1)

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    outname = f'{prefix}.png'
    fig.savefig(outname, dpi=dpi)
    plt.close(fig)
    print(f' --> {outname}')



# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    fnames = args['MSD_FILE']
    prefix = args['-o']
    dpi    = int(args['--dpi'])
    fit    = args['--fit']
    fit_start = float(args['--fit-start'])
    specorder_str = args['--specorder']
    specorder = specorder_str.split(',') if specorder_str != 'None' else []

    if not fnames:
        print('Error: no MSD file specified.', file=sys.stderr)
        sys.exit(1)

    # Parse all files
    parsed = []
    for fname in fnames:
        t, msds, errs, species = parse_msd_file(fname)
        parsed.append((t, msds, errs, species))
        print(f' Read {fname}: {len(t)} steps, species={species}')

    # Collect ordered species union
    all_species: list = []
    for _, _, _, sp in parsed:
        for s in sp:
            if s not in all_species:
                all_species.append(s)

    # Apply user-supplied specorder (reorder; append extras)
    if specorder:
        ordered = [s for s in specorder if s in all_species]
        for s in all_species:
            if s not in ordered:
                ordered.append(s)
    else:
        ordered = all_species

    # Decide mode
    if len(parsed) == 1:
        t, msds, errs, _ = parsed[0]
        mean_msds = msds
        has_err = all(errs.get(s) is not None for s in all_species)
        if has_err:
            sem_msds = {s: errs[s] for s in all_species}
            print(' Using error columns from file (--err mode).')
        else:
            sem_msds = {s: np.zeros(len(t)) for s in all_species}
            print(' Warning: no error columns; plotting without error band.')
    else:
        t, mean_msds, sem_msds = multi_file_stats(parsed, all_species)
        print(f' Computed mean ± SEM across {len(parsed)} files.')

    make_plot(t, mean_msds, sem_msds, ordered, prefix, dpi,
              fit=fit, fit_start=fit_start)


if __name__ == '__main__':
    main()
