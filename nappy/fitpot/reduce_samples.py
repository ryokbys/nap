"""
Reduce degenerate samples by looking at erg.ref.
Put redundant samples to reduced/ directory.

Usage:
  reduce.py DIRS [DIRS...]

"""
from __future__ import print_function

import os
from docopt import docopt

reduce_dir = 'reduced'


if __name__ == "__main__":
    
    args = docopt(__doc__)
    dirs = args['DIRS']
    names = []
    print('num of dirs = ',len(dirs))
    os.system('mkdir -p ./'+reduce_dir)
    nreduced = 0
    for d in dirs:
        with open(d+'/erg.ref','r') as f:
            erg= f.readline().split()[0]
        dbase = d[:-5]
        name= dbase+'_'+erg
        #print(dbase,erg,name)
        if name in names:
            redundant = True
        else:
            names.append(name)
            redundant = False
        if redundant:
            cmd = 'mv '+d+' '+reduce_dir+'/'
            print(cmd)
            os.system(cmd)
            nreduced += 1
    print('num of reduced samples = ',nreduced)
