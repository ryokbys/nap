"""
`nappy` provides python scripts related to programs in the *nap* package.
"""

import os

import vasp
import scheduler
import clutil
import espresso

_nappy_dir = '.nappy'

def get_nappy_dir():
    homedir = os.environ['HOME']
    return homedir + '/' +_nappy_dir

