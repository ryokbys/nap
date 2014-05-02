#!/usr/local/bin/python
#PBS -N fitpot_pmd
#PBS -o out
#PBS -q batch
#PBS -j oe
#PBS -l nodes=4:ppn=2
#-----------------------------------------------------------------------
# Usage:
#   $ qsub parallel_run_pmd.py
#-----------------------------------------------------------------------

import os,sys,glob,subprocess

if len(sys.argv) != 2:
    print '[Error] Number of arguments was wrong.'
    print ' Usage: ./parallel_run_pmd.py in.params.SW_Si'
    sys.exit()

fparam= sys.argv[1]

#os.system('hostname')
#os.chdir(os.environ.get('PBS_O_WORKDIR'))

dirs= glob.glob('0????')
dirs.sort()
#print dirs

#nodefname= os.environ.get('PBS_NODEFILE')
nodefname= 'nodelist.txt'
nodefile=open(nodefname,'r')
nodes=[]
for line in nodefile.readlines():
    nodes.append(line.split()[0])
nodefile.close()
#print nodes

pmdsrc= os.environ.get('HOME')+'/src/pmd'
expand_pmd= pmdsrc+'/util/expand_pmd.rb'
pmd= pmdsrc+'/src/pmd'
reduce= pmdsrc+'/fitpot/reduce_erg_frc.py'

idir= 0
done= False
while True:
    procs=[]
    nproc=0
    for node in nodes:
        if idir >= len(dirs):
            done= True
            break
        dir= dirs[idir]
        #.....create node file for pmd run
        fname='/tmp/node_for_{}'.format(dir)
        f= open(fname,'w')
        f.write(node+'\n')
        f.close()
        cmd='mpiexec -recvtimeout 100 -machinefile {}'.format(fname) \
             + ' -np 1 ./run_pmd.sh {} {}'.format(dir,fparam)
        arg= cmd.split()
        procs.append(subprocess.Popen(arg))
        nproc += 1
        idir += 1
    for i in range(nproc):
        procs[i].wait()
    if done:
        break
    
    
print " running pmd done."
