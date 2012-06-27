#!/usr/bin/python
#-----------------------------------------------------------------------
# Make work directory at /tmp/ in each to-be-used host
#-----------------------------------------------------------------------
# Usage:
#   ./mk_workdir.py numNodes workdir

import sys,os,pwd
from datetime import date

if len(sys.argv) != 3:
    print "usage: ./mk_workdir.py numNodes workdir"
    sys.exit()

numNodes= int(sys.argv[1])
workdir = sys.argv[2]

#-----decide work dir name
username= pwd.getpwuid(os.getuid())[0]
today= date.today()
#workdir= "/tmp/"+username+"/"+today.strftime("%y%m%d")
# workdir= "~/"+today.strftime("%y%m%d")

#-----list of hosts
hostfile= open("hosts.list",'r')
lines= map(lambda x:x[:-1], hostfile.readlines())
# print lines

n= 0
for host in lines:
    print "ssh %s@%s \"mkdir -p %s\"" % (username,host,workdir)
    os.system("ssh %s@%s \"mkdir -p %s\"" % (username,host,workdir))
    print "scp ./[1234]* pmd pmd.in hosts.list ini%03d %s@%s:%s/" % (n,username,host,workdir)
    os.system("scp ./[1234]* pmd pmd.in hosts.list ini%03d %s@%s:%s/" % (n,username,host,workdir))
    # os.system("scp ./* %s@%s:%s/" % (username,host,workdir))
    n=n+1
    if n >= numNodes:
        sys.exit()
