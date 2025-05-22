"""
`nappy` provides python scripts related to programs in the *nap* package.
"""
__all__ = ['pmd','napsys','atom','common','io','rdf','adf','msd','gaussian_smear',
           'util','vasp','mkcell','manipulate','units','database']

__author__ = 'RYO KOBAYASHI'
__version__ = '250522'

from . import napsys
from . import io

_nappy_dir = '.nappy'

def get_nappy_dir():
    import os
    homedir = os.environ['HOME']
    return homedir + '/' +_nappy_dir

