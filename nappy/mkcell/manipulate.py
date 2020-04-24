#!/usr/bin/env python
"""
Manipulation functions.

Usage:
  manipulate.py [options]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

from docopt import docopt

from nappy.napsys import NAPSystem

__author__ = "RYO KOBAYASHI"
__version__ = "200424"

def replace(nsys0,spc1,spc2,num=1):
    """
    Replace an atom of SPC1 in the system NSYS with SPC2.

    Parameters
    ----------
    nsys0 : NAPSystem object
           System in which atoms are to be replaced. (not inplace)
    spc1 : integer or string
           An atom ID or atom species to be replaced.
           If it is integer, the atom specified by SPC1 is replaced.
           If it is string, NUM atoms of species SPC1 are to be replaced with atoms of SPC2.
    spc2 : integer or string
           An atom ID to be replaced with spcs1, which is available only if spcs1 is integer.
           Atom species that replace atoms of spc1.
    num : integer, optional
           Number of atoms to be replaced.
    
    Returns
    -------
    nsys_out : NAPSystem object
           Result of the system in which atoms are replaced.
    """

    if type(spc1) not in (int,str):
        raise TypeError('SPC1 must be either int or str.')
    if type(spc2) not in (int,str):
        raise TypeError('SPC2 must be either int or str.')
    if type(spc2) is int and type(spc1) is not int:
        raise ValueError('SPC2 must be string if SPCS1 is string.')
    if not isinstance(nsys0,NAPSystem):
        raise ValueError('NSYS0 is not an instance of NAPSystem.')
    
    import copy
    nsys = copy.deepcopy(nsys0)

    if type(spc1) is int:
        if spc1 > nsys.num_atoms():
            raise ValueError('An integer SPC1 is greater than num of atoms in NSYS0')
        sid1 = nsys.atoms.sid[spc1]
        if type(spc2) is int:
            sid2 = nsys.atoms.sid[spc2]
            if sid1 == sid2:
                # do nothing
                return nsys
            else:
                nsys.atoms.at[spc1,'sid'] = sid2
                nsys.atoms.at[spc2,'sid'] = sid1
                return nsys
        else:
            if spc2 not in nsys.specorder:
                nsys.specorder.append(spc2)
            sid2 = nsys.specorder.index(spc2) +1
            if sid1 == sid2:
                # do nothing
                return nsys
            else:
                nsys.atoms.at[spc1,'sid'] = sid2
                return nsys
    else:  # type(spc1) is str
        import random
        if spc1 not in nsys.specorder:
            raise ValueError('The species {0:s} is not in the system.'.format(spc1))
        sid1 = nsys.specorder.index(spc1) +1
        if spc2 not in nsys.specorder:
            nsys.specorder.append(spc2)
        sid2 = nsys.specorder.index(spc2) +1
        indices = nsys.atoms[nsys.atoms.sid == sid1].index.tolist()
        for n in range(num):
            if len(indices) == 0:
                break
            idx = random.choice(indices)
            nsys.atoms.at[idx,'sid'] = sid2
            indices.remove(idx)
        return nsys
    
    raise RuntimeError('Something is wrong.')


if __name__ == "__main__":

    args = docopt(__doc__)

    print('manipulate.py does nothing...')
