#!/usr/bin/python
#-----------------------------------------------------------------------
# Gather calculation results from remote hosts
#-----------------------------------------------------------------------
# Usage:
#   ./gather_remote_data.py numNodes /fullpath/to/workdir

import sys,os,pwd

if len(sys.argv) != 3:
    print "usage: ./gather_remote_data.py numNodes /fullpath/to/workdir"
    sys.exit()

numNodes= int(sys.argv[1])
workdir= sys.argv[2]
username= pwd.getpwuid(os.getuid())[0]

#-----list of hosts
if not os.path.exists("./hosts.list"):
    print "error: hosts.list does not exist!!"
    sys.exit()

hostfile= open("hosts.list",'r')
lines= map(lambda x:x[:-1], hostfile.readlines())
# print lines

n= 0
for host in lines:
    print "scp %s@%s:%s/fin%03d ./" % (username,host,workdir,n)
    os.system("scp %s@%s:%s/fin%03d ./" % (username,host,workdir,n))
    print "scp \"%s@%s:%s/pmd%03d-*\" ./" % (username,host,workdir,n)
    os.system("scp \"%s@%s:%s/pmd%03d-*\" ./" % (username,host,workdir,n))
    #-----if node #0, receive out and out_erg
    if n==0:
        print "scp \"%s@%s:%s/out*\" ./" % (username,host,workdir)
        os.system("scp \"%s@%s:%s/out*\" ./" % (username,host,workdir))
    n=n+1
    if n>= numNodes:
        sys.exit()
