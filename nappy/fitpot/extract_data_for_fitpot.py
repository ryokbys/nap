#!/usr/bin/env python

"""
Find data directories, which include *pos*, *erg.ref*, and *frc.ref*,
under the directory specified in an argument,
and copy those directories to *data_set* directory with the names
of ascending order, such as 00001, 00002, ...

usage:
    extract_data_for_fitpot.py [options] DIR

options:
    -h, --help      show this help message and exit
    --offset=<nof>  offset number of directory name [default: 0]
    
"""
import os,shutil
import commands
from docopt import docopt

__author__= "Ryo KOBAYASHI"
__version__= "0.1"

if __name__ == "__main__":

    args= docopt(__doc__)


    dir= args['DIR']
    offset= int(args['--offset'])
    
    if not os.path.isdir(dir):
        print(' {0} is not a directory !!!'.format(dir))
        exit

    _fpos="pos"
    _ferg="erg.ref"
    _ffrc="frc.ref"
    numDir= offset
    
    sdirs= commands.getoutput('find {0} -name "{1}" | xargs -n 1 dirname'.format(dir,_fpos)).split("\n")
    print(' {0:d} dirs will be extracted...'.format(len(sdirs)))

    for sdir in sdirs:
        #...only if the directory contains erg.ref and frc.ref files
        if not ( os.path.exists(sdir+'/'+_ferg)
                 and os.path.exists(sdir+'/'+_ffrc) ):
            continue
        numDir += 1
        ndname='{0:05d}'.format(numDir)
        print('.',end='')
        os.mkdir(ndname)
        shutil.copy('{0}/{1}'.format(sdir,_fpos),ndname)
        shutil.copy('{0}/{1}'.format(sdir,_ferg),ndname)
        shutil.copy('{0}/{1}'.format(sdir,_ffrc),ndname)
        os.mkdir(ndname+'/smd')
