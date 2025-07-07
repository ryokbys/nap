#!/usr/bin/env python
"""
Convert pmd files to one extxyz file.

Usage:
  {0:s} [options] FILES [FILES...]

Options:
  -h, --help            Show this message and exit.
  -o, --output OUTPUT   Output file name. [default: from_pmds.extxyz]
  --in-format FORMAT    Input file format. [default: pmd]
"""
import os, sys
from docopt import docopt
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "250611"

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    files = args['FILES']
    outfname = args['--output']
    fmt = args['--in-format']

    nsyss = []
    for f in files:
        print('.',end='',flush=True)
        nsyss.append(nappy.io.read(f, format=fmt))

    with open(outfname, 'w') as f:
        for nsys in nsyss:
            nappy.io.write_extxyz(f, nsys)
    return None

if __name__ == "__main__":

    main()
#
#
#
