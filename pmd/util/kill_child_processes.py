#!/usr/bin/python
#-----------------------------------------------------------------------
# Kill MPI zombie processes in remote hosts
#-----------------------------------------------------------------------
# Usage:
#   ./kill_childe_processes.py numNodes workdir

import sys,os,pwd
from datetime import date

if len(sys.argv) != 3:
    print "usage: ./kill_child_processes.py numNodes workdir"
    sys.exit()

numNodes= int(sys.argv[1])
workdir = sys.argv[2]

#-----decide work dir name
username= pwd.getpwuid(os.getuid())[0]

#-----list of hosts
hostfile= open("hosts.list",'r')
lines= map(lambda x:x[:-1], hostfile.readlines())
# print lines

n= 0
for host in lines:
    print "ssh %s@%s \"kill \\`ps ux | grep %s | grep -v \'grep\' | awk \'{print \\$2}\'\\`\"" % (username,host,workdir)
    os.system("ssh %s@%s \"kill \\`ps ux | grep %s | grep -v \'grep\' | awk \'{print \\$2}\'\\`\"" % (username,host,workdir))
    #ps ux | awk '/pmd/&& !/awk/ {print $2}'
    n=n+1
    if n >= numNodes:
        sys.exit()
