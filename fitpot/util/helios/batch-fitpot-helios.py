#!/bin/env python
# -*- coding: utf-8 -*-
#SBATCH -J fitpot-GA       # job name
#SBATCH -A HatD2       # project name
#SBATCH -N 32           # number of nodes
#SBATCH -n 512           # number of tasks
#SBATCH -o out.%j      # stdout filename (%j is jobid)
#SBATCH -e err.%j      # stdout filename (%j is jobid)
#SBATCH -t 2:00:00    # execution time
#SBATCH --mail-type=END
#SBATCH --mail-user=kobayashi.ryo@nitech.ac.jp
#=======================================================================
# Submit this script to SLURM system on helios (Rokkasho supercomputer)
# as,
#   $ sbatch batch-fitpot-helios.py
# 
#=======================================================================

import os,commands,sys
import time,glob

#.....learnig or test
maindir='learning_set'
fitpot=os.environ['HOME']+'/src/nap/fitpot/fitpot.py'
infname='in.fitpot'

def check_if_GA():
    f=open(infname,'r')
    for line in f.readlines():
        if 'fitting_method' in line:
            if line.split()[1] in ('ga','GA','genetic-algorithm'):
                return True
    f.close()
    return False

def divide_nodefile_for_GA(nodefile):
    #.....get node list to be used
    f=open(nodefile,'r')
    nodes=[]
    for lin in f.readlines():
        nodes.append(lin.split()[0])
    f.close()

    #.....get num of GA individuals
    f=open(infname,'r')
    for line in f.readlines():
        if 'ga_num_individuals' in line:
            nindiv= int(line.split()[1])
    f.close()

    #.....get GA directories
    dirs= glob.glob(maindir+"_*")
    if len(dirs) < nindiv:
        print ' [Error] len(dirs) < nindiv !!!'
        print '  len(dirs), nindiv =',len(dirs),nindiv
        sys.exit()

    #.....assign to-be-computed nodes
    ndiv= len(nodes)/nindiv
    nrem= len(nodes)%nindiv
    offset= 0
    for indiv in range(nindiv):
        dirname=maindir+'_{0:05d}'.format(indiv+1)
        nnode= ndiv
        if indiv < nrem:
            nnode += 1
        nodes_for_indiv= nodes[offset:offset+nnode]
        print dirname+':', nodes_for_indiv
        #.....node-info file will be used by parallel_run_pmd.py
        g=open(dirname+'/node-info','w')
        g.write(' offset    {0:10d}\n'.format(offset))
        g.write(' num_tasks {0:10d}\n'.format(nnode))
        for node in nodes_for_indiv:
            g.write(' {0}\n'.format(node))
        g.close()
        offset += nnode

if __name__ == '__main__':
    print "{0:=^72}".format(' FITPOT ')
    t0= time.time()

    os.system("module load intel intelmpi scipy")

    os.system("srun -l /bin/hostname | sort -n | awk '{print $2}' > nodelist.txt")
    nodefile='nodelist.txt'
    os.system("cp {0:} {1:}/".format(nodefile,maindir))
    #.....node-info file will be used by parallel_run_pmd.py
    nnodes= int(commands.getoutput("wc -l {0}".format(nodefile) \
                                   +r" | awk '{print $1}'"))
    f= open(maindir+'/node-info','w')
    f.write(' offset    {0:10d}\n'.format(0))
    f.write(' num_tasks {0:10d}\n'.format(nnodes))
    f.close()
    os.system("cat {0} >> {1}".format(nodefile,maindir+'/node-info'))

    if check_if_GA():
        divide_nodefile_for_GA(nodefile)

    os.system("module load intel intelmpi scipy; " \
                  +"python {0} > out.fitpot 2>&1".format(fitpot))
    
    print '{0:=^72}'.format(' FITPOT finished correctly ')
    print '   Elapsed time = {0:12.2f}'.format(time.time()-t0)
    
