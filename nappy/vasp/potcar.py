#!/bin/env python

import re

def read_POTCAR(fname='POTCAR'):
    file= open(fname,'r')
    species=[]
    valence=[]
    encut=[]
    isp=0
    enmax='ENMAX'
    lines= file.readlines()
    for iline in range(len(lines)):
        line= lines[iline]
        if re.match(r'^\s+US ',line) or re.match('^\s+PAW_PBE ',line):
            isp=isp +1
            data= line.split()
            species.append(data[1]) # species name
            valence.append(float(lines[iline+1].rstrip()))
        if enmax in line:
            data= line.split()
            encut.append(float(data[2].rstrip(';')))
    file.close()

    potcar={}
    potcar['num_species']= isp
    potcar['species']= species
    potcar['valence']= valence
    potcar['encut']= encut
    return potcar

if __name__ == '__main__':

    potcar=read_POTCAR()
    print potcar

