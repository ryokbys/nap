#!/usr/bin/env python
"""
Reduce the number of descriptors in `in.params.desc` by looking at 
`out.fsmask.##`.

Usage:
  reduce_desc.py [options] DESC_FILE MASK_FILE

Options:
  -h, --help  Show this message and exit.
"""
from docopt import docopt
import nappy.nn.desc as desc

__author__ = "RYO KOBAYASHI"
__version__ = "180926"

def read_fsmask(fname='out.fsmask'):
    with open(fname,'r') as f:
        lines = f.readlines()
    masks = []
    for il,line in enumerate(lines):
        if il == 0:
            continue
        data = line.split()
        m = int(data[1])
        masks.append(m)
    return masks

def get_reduced_descs(nsp,nsf,descs,masks):
    descs_new = []
    for i,d in enumerate(descs):
        if masks[i] == 1:
            continue
        descs_new.append(d)
    return descs_new


if __name__ == "__main__":

    args = docopt(__doc__)
    f_desc = args['DESC_FILE']
    f_mask = args['MASK_FILE']
    of_desc = f_desc+'.new'

    nsp,nsf,descs,r_inner = desc.read_desc(f_desc)
    masks = read_fsmask(f_mask)
    descs_new = get_reduced_descs(nsp,nsf,descs,masks)
        
    print(' Number of chosen descriptors  = {0:d}'.format(len(descs_new)))
    desc.write_desc(nsp,len(descs_new),descs_new,r_inner,fname=of_desc)
    
    
