#!/bin/env python
# -*- coding: utf-8 -*-

import os,sys,glob,subprocess

if len(sys.argv) != 2:
    print '[Error] Number of arguments was wrong.'
    print ' Usage: ./parallel_run_pmd.py in.params.SW_Si'
    sys.exit()

fparam= sys.argv[1]

dirs= glob.glob('0????')
dirs.sort()
#print dirs

#.....read offset and num_nodes from node-info file
nodeinfo='node-info'
f=open(nodeinfo,'r')
for line in f.readlines():
    if 'offset' in line:
        offset= int(line.split()[1])
    if 'num_nodes' in line:
        nnodes= int(line.split()[1])
f.close()
#.....assign to-be-computed directories to each node
dir_per_node= []
ndir_per_node= len(dirs)/nnodes +1
idir= 0
done= False
for inode in range(nnodes):
    arr= []
    for i in range(ndir_per_node):
        if idir >= len(dirs):
            done= True
            break
        arr.append(dirs[idir])
        idir += 1
    if len(arr) != 0:
        dir_per_node.append(arr)
    if done:
        break

fconf= open('multi-prog.conf','w')
off= offset
for inode in range(nnodes):
    dlist= dir_per_node[inode]
    s= ""
    for d in dlist:
        s += " "+d
    fconf.write('{0:4d}'.format(off))
    fconf.write(' {0} {1}'.format('./run_pmd.sh',fparam))
    fconf.write(' {0}\n'.format(s))
    off += 1
fconf.close()

os.system("srun -n{0} --multi-prog multi-prog.conf".format(nnodes))
