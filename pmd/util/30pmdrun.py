#!/usr/bin/python
#-----------------------------------------------------------------------
# Execute PMD using MPIRUN command
#   with reading parallel nodes configuration from 'in.control'
#-----------------------------------------------------------------------
# Usage:
#   ./30pmdrun.py

import os

# file= open("pmd.in",'r')
# line= file.readline()
# file.close()
# 
# nodes= line.split()
# nx= int(nodes[0])
# ny= int(nodes[1])
# nz= int(nodes[2])
nx= 1
ny= 1
nz= 1
nxyz= nx*ny*nz

print "---> Execute PMD using MPIRUN command with nodes..."
print "  nx,ny,nz,nyxz= %d %d %d %d" % (nx,ny,nz,nxyz)

#-----execute
print " mpirun -np %d -machinefile hosts.list ./pmd" % (nxyz)
os.system("mpirun -np %d -machinefile hosts.list ./pmd" % (nxyz))
