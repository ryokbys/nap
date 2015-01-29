#!/usr/local/bin/python

import glob,random,os
import optparse

usage= '%prog [options] '

parser= optparse.OptionParser(usage=usage)
parser.add_option("-n",dest="nlearn",type="int",default=1,
                  help="number of samples as a learning set.")
parser.add_option("--dir-learn",dest="dir_learn",type="string",
                  default="learning_set",
                  help="name of the directory for a learning set.")
parser.add_option("--dir-test",dest="dir_test",type="string",
                  default="test_set",
                  help="name of the directory for a test set.")
(options,args)= parser.parse_args()

#.....number of samples in learning_set
num_learning_set= options.nlearn
dir_learning_set= options.dir_learn
dir_test_set= options.dir_test

if not os.path.exists(dir_learning_set):
    os.system('mkdir '+dir_learning_set)
if not os.path.exists(dir_test_set):
    os.system('mkdir '+dir_test_set)

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
