#!/usr/bin/env python
"""
Extract an equilibrated cell from an extxyz MD trajectory.

Usage:
  {0:s} [options] FILENAME

Options:
  -h, --help            Show this message and exit.
  -o, --output OUTPUT   Output file name. [default: pmdini_equil]
"""
import copy
import os
import sys

from docopt import docopt
import numpy as np

import nappy.io

__author__ = "RYO KOBAYASHI"
__version__ = "260707"


def load_trajectory(fname):
    nsyss = nappy.io.read(fname=fname, format='extxyz')
    if isinstance(nsyss, list):
        return nsyss
    return [nsyss]


def make_equilibrated_system(nsyss):
    if len(nsyss) == 0:
        raise ValueError('No structure was read from the input trajectory.')

    start = len(nsyss) // 2
    hmat_avg = np.mean([nsys.get_hmat() for nsys in nsyss[start:]], axis=0)

    nsys_out = copy.deepcopy(nsyss[-1])
    nsys_out.set_hmat(hmat_avg)
    return nsys_out, start


def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    infname = args['FILENAME']
    outfname = args['--output']

    nsyss = load_trajectory(infname)
    nsys_out, start = make_equilibrated_system(nsyss)
    nappy.io.write(nsys_out, fname=outfname, format='pmd')

    print(f'Input file: {infname}')
    print(f'Number of frames: {len(nsyss)}')
    print(f'Averaged hmat frames: {start + 1} - {len(nsyss)}')
    print(f'Output file: {outfname}')
    return None


if __name__ == "__main__":
    main()
