#!/bin/sh
#SBATCH -J neb        # job name
#SBATCH -N 34           # number of nodes
#SBATCH -n 544           # number of tasks
#SBATCH -o out.%j      # stdout filename (%j is jobid)
#SBATCH -e err.%j      # stdout filename (%j is jobid)
#SBATCH -t 2:00:00    # execution time

. /etc/profile.d/00-modules.sh
module load intelmpi # load module

MACHINEFILE="hostfile"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE
NUM_HOSTS=${SLURM_JOB_NUM_NODES}

nnodes=${SLURM_JOB_NUM_NODES}
ntasks=${SLURM_NTASKS}

# #.....mpdboot for intelmpi
# cat $MACHINEFILE | uniq > ./mpd.host
# mpdboot -n ${SLURM_JOB_NUM_NODES} -f ./mpd.host
# #srun -N ${nnodes} "ps ux | grep mpd"

#cwd=$WORKDIR
#echo "cwd= $cwd"

# copy the binary on /tmp on each node
if [ -L ./pmd.local ]; then rm ./pmd.local; fi
srun -N ${nnodes} -n ${ntasks} cp ./pmd_Ito_WHe /tmp/pmd
ln -s /tmp/pmd ./pmd.local

if [ -L ./vasp.local ]; then rm ./vasp.local; fi
srun -N ${nnodes} -n ${ntasks} cp ./vasp /tmp/vasp
ln -s /tmp/vasp ./vasp.local

if [ -L ./qmcl.local ]; then rm ./qmcl.local; fi
srun -N ${nnodes} -n ${ntasks} cp ./qmcl /tmp/qmcl
ln -s /tmp/qmcl ./qmcl.local

if [ -L ./neb.local ]; then rm ./neb.local; fi
srun -N ${nnodes} -n ${ntasks} cp ./neb /tmp/neb
ln -s /tmp/neb ./neb.local

srun -N 1 -n 1 --relative=0 ./neb.local > out.neb

# #.....stop mpd
# mpdallexit
