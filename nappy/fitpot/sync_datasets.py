#!/usr/bin/env python
"""
Syncronize given two datasets.

Usage:
  sync_datasets [options] DSET1 DSET2

Options:
  -h, --help  Show this message and exit.
  --dry       Dry run.
"""
import os
from docopt import docopt
import glob

__author__ = "RYO KOBAYASHI"
__version__ = "170608"

if __name__ == "__main__":

    args = docopt(__doc__)
    dset1 = args['DSET1']
    dset2 = args['DSET2']
    dry = args['--dry']

    list1 = [ os.path.basename(d) for d in glob.glob(dset1+'/smpl_*')] 
    list2 = [ os.path.basename(d) for d in glob.glob(dset2+'/smpl_*')] 

    from1to2 = []
    for d1 in list1:
        if d1 not in list2:
            from1to2.append(d1)
    from2to1 = []
    for d2 in list2:
        if d2 not in list1:
            from2to1.append(d2)

    print('Copy files from {0:s} to {1:s}: {2:d}'.format(dset1,dset2,len(from1to2)))
    for d12 in from1to2:
        print('  '+d12)
        if not dry:
            os.system('cp -r {0:s}/{1:s} {2:s}/'.format(dset1,d12,dset2))
    print('Copy files from {0:s} to {1:s}: {2:d}'.format(dset2,dset1,len(from2to1)))
    for d21 in from2to1:
        print('  '+d21)
        if not dry:
            os.system('cp -r {0:s}/{1:s} {2:s}/'.format(dset2,d21,dset1))

    print('done')
