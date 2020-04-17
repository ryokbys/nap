import fnmatch
import functools
import os
import re

def find(path, pattern):
    """
    Emulation of Linux command `find`.
    See:  http://d.hatena.ne.jp/fgshun/20080901/1220272713
    """
    if isinstance(pattern, str):
        match = functools.partial(fnmatch.fnmatch, pat=pattern)
    elif isinstance(pattern, re._pattern_type):
        match = pattern.match
    elif callable(pattern):
        match = pattern
    else:
        raise TypeError

    for root, dirs, files in os.walk(path):
        for file_ in files:
            if match(file_):
                yield os.path.join(root, file_)

