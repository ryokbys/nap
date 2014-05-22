#!/usr/bin/python
#-----------------------------------------------------------------------
# Rename 'fin###' files to 'ini###'
#-----------------------------------------------------------------------
# Usage:
#   ./rename_remote.py numNodes workdir

import sys,glob,os,pwd

if len(sys.argv) != 3:
    print "usage: ./rename_remote.py numNodes workdir"
    sys.exit()

numNodes= int(sys.argv[1])
workdir = sys.argv[2]

#-----username
username= pwd.getpwuid(os.getuid())[0]
#-----hosts
hostfile= open("hosts.list",'r')
lines= map(lambda x:x[:-1], hostfile.readlines())

n=0
for host in lines:
    print "ssh %s@%s \"cd %s; mv fin%03d ini%03d\"" % (username,host,workdir,n,n)
    os.system("ssh %s@%s \"cd %s; mv fin%03d ini%03d\"" % (username,host,workdir,n,n))
    n=n+1
    if n >= numNodes:
        sys.exit()
