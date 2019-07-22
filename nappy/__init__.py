"""
`nappy` provides python scripts related to programs in the *nap* package.
"""
from __future__ import print_function

import os
from . import *
# from . import vasp
# from . import scheduler
# from . import clutil
# from . import espresso
# from . import mkcell

_nappy_dir = '.nappy'

def get_nappy_dir():
    homedir = os.environ['HOME']
    return homedir + '/' +_nappy_dir

