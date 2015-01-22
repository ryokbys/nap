#!/usr/bin/python -u
#
#=========================================================================
# This script must be loaded from fitpot program.
#=========================================================================

import os,sys,glob,subprocess


################################################# FUNCTIONS ############

# def prepare_pmd(dirs,fparam):
#     #.....replace cutoff radius in in.pmd
#     f= open(fparam,'r')
#     rcut= float(f.readline().split()[1])
#     f.close()
#     os.system('copy in.pmd in.pmd.tmp')
#     finpmd=open('in.pmd.tmp','r')
#     finpmdo= open('in.pmd','w')
#     for line in finpmd.readline():
#         if 'cutoff_radius' in line:
#             d= line.split()
#             finpmdo.write(' {0} {1:8.3f}\n'.format(d[0],d[1]))
#         else:
#             finpmdo.write('{0}'.format(line))
#     finpmd.close()
#     finpmdo.close()
# 
#     #.....prepare for pmd run
#     for dir in dirs:
#         os.system('mkdir -p {0}/pmd/0000'.format(dir))
#         os.system('cp {0} {1}/pmd/'.format('in.pmd',dir))
#         os.system('cp {0} {1}/pmd/'.format(fparam,dir))
#         os.system('')

################################################# MAIN ROUTINE #########

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print '[Error] Number of arguments was wrong.'
        print ' Usage: ./parallel_run_pmd.py in.params.SW_Si'
        sys.exit()
    
    fparam= sys.argv[1]
    
    dirs= glob.glob('[0-9]????')
    dirs.sort()
    #print dirs
    
    rankfname= 'ranklist.txt'
    if not os.path.exists(rankfname):
        print ' [Error] {0} does not exist !!!'.format(rankfname)
        sys.exit()
    rankfile=open(rankfname,'r')
    ranks=[]
    for line in rankfile.readlines():
        ranks.append(line.split()[0])
    rankfile.close()
    
    #...assign to-be-computed directories to each rank
    dir_per_rank= []
    ndirs= [ 0 for i in range(len(ranks))]
    nrem= len(dirs)%len(ranks)
    for i in range(len(ranks)):
        ndirs[i]= len(dirs)/len(ranks)
        if i < nrem:
            ndirs[i] += 1
    idir= 0
    done= False
    for irank in range(len(ranks)):
        arr= []
        for i in range(ndirs[irank]):
            if idir >= len(dirs):
                done= True
                break
            arr.append(dirs[idir])
            idir += 1
        if len(arr) != 0:
            dir_per_rank.append(arr)
        if done:
            break
    
    procs= []
    for irank in range(len(ranks)):
        rank= ranks[irank]
        dir_list= dir_per_rank[irank]
        str= ""
        for dir in dir_list:
            str += " "+dir
        #...create rank file for pmd run
        fname='/tmp/rankfile_{0}'.format(rank)
        f= open(fname,'w')
        f.write(rank+'\n')
        f.close()
        #...run run_pmd.sh on the remote rank
        # cmd='mpirun --hostfile {}'.format(fname) \
        #      + ' -np 1 ./run_pmd.sh {} {}'.format(fparam,str)
        #os.system('scp {0} {1}:{2}/'.format(fparam,rank,os.getcwd()))
        #cmd='ssh -q {0} "cd {1} && ./run_pmd.sh {2} {3}"'.format(rank,os.getcwd(),fparam,str)
        cmd='./run_pmd.sh -p {0} {1} {2}'.format(irank,fparam,str)
        procs.append(subprocess.Popen(cmd,shell=True))
    for i in range(len(procs)):
        procs[i].wait()
    
    # print " running pmd done."
