#!/bin/env python
from __future__ import print_function

import re

def read_POTCAR(fname='POTCAR'):
    species=[]
    valence=[]
    encut=[]
    isp=0

    with open(fname,'r') as f:
        lines= f.readlines()
        for iline in range(len(lines)):
            line= lines[iline]
            if (re.match(r'^\s+US ',line) or
                re.match(r'^\s+PAW_(PBE|GGA) ',line) or re.match(r'^\s+PAW ',line)) and \
               'radial sets' not in line:
                isp=isp +1
                data= line.split()
                species.append(data[1])  # species name
                valence.append(float(lines[iline+1].rstrip()))
            if 'ENMAX' in line:
                data= line.split()
                encut.append(float(data[2].rstrip(';')))

    potcar={}
    potcar['num_species']= isp
    potcar['species']= species
    potcar['valence']= valence
    potcar['encut']= encut
    return potcar


if __name__ == '__main__':

    potcar=read_POTCAR()
    print(potcar)

