#!/usr/bin/python
#-----------------------------------------------------------------------
# Execute PMD2POT using MPIRUN command
#   with reading parallel nodes configuration from 'pmd.in'
#-----------------------------------------------------------------------
# Usage:
#   ./41pmd2pot.py

import os

file= open("pmd.in",'r')
line= file.readline()
file.close()

nodes= line.split()
nx= int(nodes[0])
ny= int(nodes[1])
nz= int(nodes[2])
nxyz= nx*ny*nz

print "---> Execute PMD2POT using MPIRUN command with nodes..."
print "  nx,ny,nz,nyxz= %d %d %d %d" % (nx,ny,nz,nxyz)

#-----execute
os.system("mpirun -np %d -machinefile hosts.list ./pmd2pot" % (nxyz))
