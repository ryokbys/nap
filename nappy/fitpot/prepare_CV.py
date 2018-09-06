#!/usr/bin/env python
"""
Prepare cross validation for fitpot.

Usage:
  prepare_CV.py [options]

Options:
  -h, --help  Show this message and exit.
  -n N        Number for N-fold CV. [default: 10]
  -b BATCH    Batch file name to be copied to each CV directory. [default: batch.fitpot.sh]
  -p PREFIX   Prefix to be added to the batch title. [default: ]
"""
from __future__ import print_function

__author__ = "RYO KOBAYASHI"
__version__ = "180614"

def maindir_from_input(fname='in.fitpot'):
    """
    Read in.fitpot to get main_directory.
    """
    with open(fname,'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'main_directory' in line:
            maindir = line.split()[1]
    return maindir

def paramfname_from_input(fname='in.fitpot'):
    with open(fname,'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'param_file' in line:
            param_file = line.split()[1]
    return param_file

def make_new_input(fname='in.fitpot',dname='CV_0'):
    """
    Make new in.fitpot in the given directory by changing as follows:
      - add "../" before the current main_directory
      - add "sample_list  dir_list.txt" entry
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    if not os.path.exists(dname):
        raise IOError('There is not a directory: {0:s}.'.format(dname))

    with open(dname+'/'+fname,'w') as f:
        for il,line in enumerate(lines):
            data = line.split()
            if len(data) == 0:
                f.write(line)
            elif data[0] == 'main_directory':
                tmp = data[1].strip('"')
                tmp = '"../'+tmp+'"'
                f.write('{0:s}   {1:s}\n'.format(data[0],tmp))
            elif 'sample_list' in line:
                f.write('sample_list   dir_list.txt\n')
            else:
                f.write(line)
    return

def make_new_batch(fname='batch.fitpot.sh',dname='CV_0',prefix=''):
    """
    Make new batch file in the given directory by replacing the batch titile.
    """
    with open(fname,'r') as f:
        lines = f.readlines()

    if not os.path.exists(dname):
        raise IOError('There is not a directory: {0:s}.'.format(dname))

    with open(dname+'/'+fname,'w') as f:
        for il,line in enumerate(lines):
            data = line.split()
            if len(data) == 0:
                f.write(line)
            elif data[0] == '#PJM' and data[1] == '-N':
                f.write('{0:s} {1:s} {2:s}\n'.format(data[0],data[1],
                                                     prefix+dname))
            else:
                f.write(line)
    return
    

if __name__ == "__main__":
    
    import os
    import glob
    import random
    from docopt import docopt
    
    args = docopt(__doc__)
    ncv = int(args['-n'])
    batchname = args['-b']
    prefix = args['-p']

    if not os.path.exists('./in.fitpot'):
        raise IOError('There is not an in.fitpot file.')
    if not os.path.exists(batchname):
        raise IOError('There is not batch file: {0:s}.'.format(batchname))
    print('Prepare for {0:d}-fold CV.'.format(ncv))

    maindir = maindir_from_input('in.fitpot').replace('"','')
    paramfname = paramfname_from_input('in.fitpot')
    
    if not os.path.exists(maindir):
        raise ValueError('There is no directory, ',maindir)

    smpldirs = glob.glob(maindir+'/smpl_*')
    random.shuffle(smpldirs)
    for i in range(len(smpldirs)):
        smpldirs[i] = smpldirs[i].replace(maindir+'/','')
    icv = []
    for i in range(len(smpldirs)):
        j = i%ncv
        icv.append(j)
    for i in range(ncv):
        dname = 'CV_{0:d}'.format(i)
        fname = 'dir_list.txt'
        os.system('mkdir -p {0:s}'.format(dname))
        with open(dname+'/'+fname,'w') as f:
            for j,s in enumerate(smpldirs):
                ic = 1
                if icv[j] == i:
                    ic = 2
                f.write('{0:s}  {1:d}\n'.format(s,ic))
        print('Make {0:s}/{1:s}'.format(dname,fname))

        make_new_input('in.fitpot',dname)
        make_new_batch(batchname,dname,prefix)
        os.system('cp {0:s} {1:s}/'.format(paramfname,dname))
        
    msg = """
Check if in.fitpot file in each CV_# directory is correct:
- main_directory (to add "../")
- sample_list   dir_list.txt

And confirm that the following files are copied to each CV_# directory:
- in.fitpot
- in.vars.fitpot (or in.params.....)
- batch.fitpot.sh (if needed)
"""
    print(msg)



    
