#!/usr/bin/python

import glob,random,os


dirs= glob.glob('[0-9]????')
dirs.sort()

nas = [ 0 for i in range(len(dirs))]

print 'len(nas)=',len(nas)

newdirs= [ '{0:05d}'.format(i+1) for i in range(len(dirs))]
newtmpdir='newtmpdir/'
os.system('mkdir -p '+newtmpdir)

length= len(nas)
for i in range(length):
    r= random.random()
    j= int(r*len(dirs))
    print 'mv '+dirs[j]+' '+newtmpdir+newdirs[i]
    os.system('mv '+dirs[j]+' '+newtmpdir+newdirs[i])
    dirs.pop(j)
    
