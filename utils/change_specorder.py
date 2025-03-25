#!/usr/bin/env python
"""
Change specorder of the given pmd-format file.

Usage:
  {0:s} [options] FILENAME

Options:
  -h, --help             Show this message and exit.
  -v, --verbose          Verbose output. [default: False]
  -o, --output OUTFNAME  Output file name. [default: None]
  --specorder SPECORDER  New specorder in comma-separated format (e.g., Si,O). [default: None]
"""
import os,sys
from docopt import docopt
import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "250325"

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    fname = args['FILENAME']
    ofname = args['--output']
    if ofname == 'None':
        ofname = fname
    print(f" Input file: {fname}")
    print(f" Output file: {ofname}")
    if args['--verbose']:
        print(" Verbose mode on.")

    try:
        nsys = nappy.io.read(fname, format='pmd')
    except:
        raise
    print(f" Old specorder: {nsys.specorder}")
    specorder = args['--specorder']
    if specorder == 'None':
        raise ValueError('--specorder option must be set.')
    specorder = [ s for s in specorder.split(',') ]
    print(f" New specorder: {specorder}")
    nsys.set_specorder(*specorder)
    auxs = []
    if 'fx' in nsys.atoms.columns.values:
        auxs = ['fx','fy','fz']
    nappy.io.write(nsys, ofname, format='pmd', auxs=auxs)
    return None

if __name__ == "__main__":

    main()
#
#
#
