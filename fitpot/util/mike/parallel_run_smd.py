#!/usr/local/bin/python

import os,sys,glob,subprocess

usage="""
Usage: python parallel_run_smd.py in.params.NN
       python parallel_run_smd.py in.params.NN smpl_Mg1 smpl_Mg2 ...
"""

if len(sys.argv) < 2:
    print '[Error] Number of arguments was wrong.'
    print usage
    sys.exit()
elif len(sys.argv) == 2:
    fparam= sys.argv[1]
    dirs= glob.glob('smpl_*')
elif len(sys.argv) > 2:
    fparam= sys.argv[1]
    dirs = []
    for iarg,arg in enumerate(sys.argv):
        if iarg < 2:
            continue
        dirs.append(arg)
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

uniqnodes = []
for node in nodes:
    if node not in uniqnodes:
        uniqnodes.append(node)

#...assign to-be-computed directories to each node
dir_per_node= []
ndirs= [ 0 for i in range(len(nodes))]
nrem= len(dirs)%len(nodes)
for i in range(len(nodes)):
    ndirs[i]= len(dirs)/len(nodes)
    if i < nrem:
        ndirs[i] += 1
idir= 0
done= False
for inode in range(len(nodes)):
    arr= []
    for i in range(ndirs[inode]):
        if idir >= len(dirs):
            done= True
            break
        arr.append(dirs[idir])
        idir += 1
    if len(arr) != 0:
        dir_per_node.append(arr)
    if done:
        break

for node in uniqnodes:
    os.system('scp {0} {1}:{2}/'.format(fparam,node,os.getcwd()))

procs= []
for inode in range(len(dir_per_node)):
    node= nodes[inode]
    dir_list= dir_per_node[inode]
    str= ""
    for dir in dir_list:
        str += " "+dir
    #...create node file for smd run
    fname='/tmp/nodefile_{0}'.format(node)
    f= open(fname,'w')
    f.write(node+'\n')
    f.close()
    #...run run_smd.sh on the remote node
    # cmd='mpirun --hostfile {}'.format(fname) \
    #      + ' -np 1 ./run_smd.sh {} {}'.format(fparam,str)
    #os.system('scp {0} {1}:{2}/'.format(fparam,node,os.getcwd()))
    cmd='ssh -q {0} "cd {1} && ./run_smd.sh {2} {3}"'.format(node,os.getcwd(),fparam,str)
    procs.append(subprocess.Popen(cmd,shell=True))
for i in range(len(procs)):
    procs[i].wait()

# print " running smd done."
