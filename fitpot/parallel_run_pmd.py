#!/usr/local/bin/python

import os,sys,glob,subprocess

if len(sys.argv) != 2:
    print '[Error] Number of arguments was wrong.'
    print ' Usage: ./parallel_run_pmd.py in.params.SW_Si'
    sys.exit()

fparam= sys.argv[1]

dirs= glob.glob('0????')
dirs.sort()
#print dirs

nodefname= 'nodelist.txt'
if not os.path.exists(nodefname):
    print ' [Error] {0} does not exist !!!'.format(nodefname)
    sys.exit()
nodefile=open(nodefname,'r')
nodes=[]
for line in nodefile.readlines():
    nodes.append(line.split()[0])
nodefile.close()

#...assign to-be-computed directories to each node
dir_per_node= []
ndir_per_node= len(dirs)/len(nodes)
nrem= len(dirs)%len(nodes)
ndirs=[]
for i in range(len(nodes)):
    ndirs.append( ndir_per_node )
    if i < nrem:
        ndirs[i] += 1
idir= 0
done= False
for n in range(len(nodes)):
    node= nodes[n]
    arr= []
    for i in range(ndirs[n]):
        if idir >= len(dirs):
            done= True
            break
        arr.append(dirs[idir])
        idir += 1
    if len(arr) != 0:
        dir_per_node.append(arr)
    if done:
        break

procs= []
for inode in range(len(nodes)):
    node= nodes[inode]
    dir_list= dir_per_node[inode]
    str= ""
    for dir in dir_list:
        str += " "+dir
    #...create node file for pmd run
    fname='/tmp/nodefile_{0}'.format(node)
    f= open(fname,'w')
    f.write(node+'\n')
    f.close()
    #...run run_pmd.sh on the remote node
    # cmd='mpirun --hostfile {}'.format(fname) \
    #      + ' -np 1 ./run_pmd.sh {} {}'.format(fparam,str)
    cmd='ssh {0} "cd {1} && ./run_pmd.sh {2} {3}"'.format(node,os.getcwd(),fparam,str)
    procs.append(subprocess.Popen(cmd,shell=True))
for i in range(len(procs)):
    procs[i].wait()

# print " running pmd done."
