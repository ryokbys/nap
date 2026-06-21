#!/usr/bin/env python
"""
Plot RDF data from an out.rdf file produced by nappy/rdf.py.
Generates graph_rdf_total.png (total g(r)) and, for multi-species systems,
graph_rdfs.png (pairwise g(r) subplot grid).

Usage:
  plot_rdf.py [options] RDFFILE

Options:
  -h, --help              Show this help message and exit.
  --specorder=SPECORDER   Species in order, comma-separated, e.g. Li,Ge,P,S.
                          Used to arrange the pairwise subplot grid. If omitted,
                          order is inferred from the file. [default: None]
  -o PREFIX               Prefix for output image filenames. [default: graph_rdf]
  --dpi=DPI               Image resolution in dots per inch. [default: 150]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from docopt import docopt

__author__ = "nappy-rdf skill"
__version__ = "260620"


def parse_rdf(fname):
    """
    Parse an out.rdf file written by nappy/rdf.py write_rdf_normal().

    File format
    -----------
    Header (line 0):
        # 1:rd[i],  2:all-all,  3:Li-Li,  4:Li-Ge, ...
    Data lines:
        r  g_total  g_pair1  g_pair2  ...

    Returns
    -------
    rd   : 1-D numpy array of r values (Angstrom)
    pairs: list of (sp1, sp2) string tuples; first entry is ('all', 'all')
    rdfs : dict mapping each pair tuple to a 1-D numpy array of g(r)
    """
    with open(fname) as f:
        lines = f.readlines()

    # Header tokens: ['#', '1:rd[i],', '2:all-all,', '3:Li-Li,', ...]
    header_tokens = lines[0].split()
    pairs = []
    for tok in header_tokens[2:]:
        label = tok.strip(',').split(':')[-1]   # 'all-all', 'Li-Li', ...
        pairs.append(tuple(label.split('-')))

    data = np.loadtxt(fname, comments='#')
    rd = data[:, 0]
    rdfs = {pair: data[:, i + 1] for i, pair in enumerate(pairs)}

    return rd, pairs, rdfs


def plot_total(rd, rdfs, prefix, dpi):
    """Plot total g(r) and save to {prefix}_total.png."""
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.plot(rd, rdfs[('all', 'all')], color='steelblue', lw=1.5)
    ax.set_xlabel('Distance (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title('Total RDF')
    ax.set_xlim(rd[0], rd[-1])
    ax.set_box_aspect(1)
    plt.tight_layout()
    outname = f'{prefix}_total.png'
    fig.savefig(outname, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f' --> {outname}')


def plot_pairwise(rd, rdfs, pairs, specorder, prefix, dpi):
    """
    Plot pairwise g_{ij}(r) in an upper-triangular subplot grid and save
    to {prefix}s.png.
    """
    pairwise = [p for p in pairs if p != ('all', 'all')]
    if not pairwise:
        return

    if specorder is None:
        # Infer unique species from pair labels, preserving first-seen order
        seen = []
        for p in pairwise:
            for s in p:
                if s not in seen:
                    seen.append(s)
        specorder = seen

    nsp = len(specorder)
    fig, axes = plt.subplots(
        nsp, nsp,
        figsize=(5 * nsp, 5 * nsp),
        dpi=dpi,
        sharex=True,
    )
    if nsp == 1:
        axes = np.array([[axes]])

    for i, si in enumerate(specorder):
        for j, sj in enumerate(specorder):
            ax = axes[i, j]
            if j < i:
                ax.axis('off')
                continue
            grdata = rdfs.get((si, sj))
            if grdata is None:
                grdata = rdfs.get((sj, si))
            if grdata is None:
                ax.axis('off')
                continue
            ax.plot(rd, grdata, color='steelblue', lw=1.5)
            ax.set_box_aspect(1)
            ax.text(0.05, 0.85, f'{si}-{sj}',
                    transform=ax.transAxes, ha='left', fontsize=10)

    fig.supxlabel('Distance (Å)')
    fig.supylabel('g(r)')
    plt.tight_layout()
    outname = f'{prefix}s.png'
    fig.savefig(outname, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f' --> {outname}')


def main():
    args = docopt(__doc__, version=__version__)
    fname  = args['RDFFILE']
    prefix = args['-o']
    dpi    = int(args['--dpi'])
    sp_str = args['--specorder']

    specorder = None if sp_str == 'None' else sp_str.split(',')

    sns.set_theme(context='talk', style='ticks')

    rd, pairs, rdfs = parse_rdf(fname)
    plot_total(rd, rdfs, prefix, dpi)
    plot_pairwise(rd, rdfs, pairs, specorder, prefix, dpi)

    return None


if __name__ == '__main__':
    main()
