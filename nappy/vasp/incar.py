# -*- coding: utf-8 -*-
"""
Utility functions related to INCAR file.
"""
from __future__ import print_function

_keywords = (
    'SYSTEM',
    'ISTART', 'ICHARG','PREC','ALGO','ISPIN',
    'NWRITE','LWAVE','LCHARGE',
    'NBAND','MAMGMON',
    'ENCUT','EDIFF','NELM','NELMIN','LREAL',
    'IMIX','AMIX','BMIX','AMIX_MAG','BMIX_MAG','AMIN',
    'IBRION','NSW','ISIF','ISYM','EDIFFG',
    'ISMEAR','SIGMA','LORBIT',
    'NPAR','NCORE',
)

_comment_char = ('#',)

def parse_INCAR(fname='INCAR'):
    """
    Parse INCAR file and return dictionary of INCAR parameters.
    """

    incar_dic = {}
    with open(fname,'r') as f:
        lines = f.readlines()
    for line in lines:
        entry = line.split()
        #...Check if the line contains any keyword
        contains_keyword = False
        for kw in _keywords:
            if kw in line:
                contains_keyword = True
                break
        if not contains_keyword:
            continue
        #...check if it is a comment line
        if entry[0][0] is '#':
            continue
        print('entry =',entry)
        #...If else, start parsing the line
        for kw in _keywords:
            print(' entry[0],kw =',entry[0],kw)
            if entry[0] == kw:
                value = parse_entry(entry[1])
                if value is None:
                    value = parse_entry(entry[2])
                incar_dic[kw] = value
    return incar_dic


def parse_entry(entry):
    """
    Parse entry value and return int, float, str, or bool.
    Return None if the value is '='.
    """
    if entry[0] is '=':
        return None
    if entry[0].isdigit():  # int or float
        if '.' in entry or 'e' in entry:  # float
            return float(entry)
        else:
            return int(entry)
    elif entry in ('True','False'):  # bool
        return bool(entry)
    else:  # str
        return entry
    raise ValueError('entry seems wrong: ',entry)


