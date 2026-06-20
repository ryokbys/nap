#!/usr/bin/env python
"""
What is this python script about...

Usage:
  {0:s} [options] FILENAME

Options:
  -h, --help            Show this message and exit.
  -v, --verbose         Verbose output. [default: False]
  -o, --output OUTPUT   Output file name. [default: outfile]
"""

import os,sys
from docopt import docopt
# import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "YYMMDD"

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])),
                  version=__version__)
    
    print(f"Input file: {args['FILENAME']}")
    if args['--verbose']:
        print("Verbose mode on.")

    print(f"Output file: {args['--output']}")


if __name__ == "__main__":

    main()
#
#
#
