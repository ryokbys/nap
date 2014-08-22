#!/usr/local/bin/python

import glob,random,os

#.....number of samples in learning_set
num_learning_set= 300
dir_learning_set= '../learning_set'
dir_test_set= '../test_set'

dirs= glob.glob('[0-9]????')
dirs.sort()

ratio= float(num_learning_set)/len(dirs)
print ' num_learning_set  =', num_learning_set
print ' num_dirs          =', len(dirs)
print ' learning_set_ratio=', ratio

num= 0
for dir in dirs:
    rnd= random.random()
    if rnd < ratio:
        dest= dir_learning_set
        num += 1
    else:
        dest= dir_test_set
    print ' {0} ---> {1}'.format(dir,dest)
    os.system('cp -r {0} {1}/'.format(dir,dest))

print ' num of samples in learning set =',num
