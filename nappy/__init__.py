"""
`nappy` provides python scripts related to programs in the *nap* package.
"""
from __future__ import print_function

__all__ = ['napsys','atom','common','io','rdf','adf','msd','gaussian_smear',
           'util','vasp','mkcell']

from . import *

_nappy_dir = '.nappy'

def get_nappy_dir():
    import os
    homedir = os.environ['HOME']
    return homedir + '/' +_nappy_dir

