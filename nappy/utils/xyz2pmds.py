#!/usr/bin/env python
"""
Convert extxyz file to series of pmd files.

Usage:
  {0:s} [options] XYZFILE

Options:
  -h, --help       Show this message and exit.
  --postfix PFIX   Postfix for output pmd file like pmd_POSTFIX_####. [default None]
  --specorder SPECORDER
                   Specorder for all the systems read from the extxyz file.
                   The list should be given as a comma-separated list. [default: None]
"""
import os, sys
from docopt import docopt
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "250614"

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)

    xyzfile = args['XYZFILE']
    postfix = args['--postfix']
    if postfix == None:
        postfix = ''
    else:
        if postfix[-1] != '_':
            postfix += '_'

    specorder = args['--specorder'].split(',')
    assert type(specorder) in (list, tuple), ' Error: specorder must be specified as comma-separated list.'

    nsyss = nappy.io.read(xyzfile, format='extxyz', specorder=specorder)

    for i, nsys in enumerate(nsyss):
        fname = f'pmd_{postfix}{i:05d}'
        nappy.io.write(nsys, fname=fname, format='pmd')

    return None

if __name__ == "__main__":

    main()
#
#
#
