#!/usr/bin/env python
"""
Manage DFT calculations on some cluster.
- Search directories under given directory whether there are directories
  in which a DFT calculation should be done.
- Compute nodes required to compute those DFT calculations.
- Create job scripts for each submission.

Usage:
  clmgr.py [options] DIRS [DIRS..]

Options:
  -h, --help  Show this message and exit.
  -c CALCULATOR
              Set calculator type. Available types are: vasp
              [default: vasp]
  -d          Dryrun: run clmgr without actually submitting jobs.
              [default: False]
  -j NUM_JOBS
              Number of jobs in one submission.
              In case that there is a limitation to the number of jobs 
              that can run at the same time in the server, you may 
              had better set this larger than 1 to reduce the number of 
              submission by increasing number of jobs in one submission.
              [default: 1]
  -q QUEUE_NAME
              Queue name the jobs are submitted to. [default: default]
"""
from __future__ import print_function

import os
from docopt import docopt
from subprocess import Popen,PIPE
from datetime import datetime
import time
import copy
from logging import getLogger, StreamHandler, FileHandler, DEBUG, WARNING

import nappy
from nappy.clutil.find import find

__author__ = "RYO KOBAYASHI"
__version__ = "170122"

_clmgr_dir = nappy.get_nappy_dir()+'/clmgr'


def pid_exists(pid):
    username = os.environ['USER']
    p = Popen(['ps','-u {0:s}'.format(username)],stdout=PIPE,stderr=PIPE)
    output = p.communicate()[0].splitlines()
    for l in output:
        dat = l.split()
        if not dat[0].isdigit():
            continue
        pid1 = int(dat[1])
        if pid == pid1:
            return True
    return False



class ClusterManager(object):
    
    _stat_file = '.clmgr_stat'
    _pid_file = _clmgr_dir +'/pid'
    _date = datetime.now().strftime("%y%m%d")
    _log_file = _clmgr_dir +'/log_{0:s}_{1:d}'.format(_date,os.getpid())
    _batch_fname = 'batch.clmgr.sh'
    
    def __init__(self,calculator=None,queue=None):
        nappydir = nappy.get_nappy_dir()
        if not os.path.exists(nappydir):
            os.mkdir(nappydir)

        if not os.path.exists(_clmgr_dir):
            os.mkdir(_clmgr_dir)

        self.logger = getLogger(__name__)

        # Set calculator
        self.set_calculator(calculator)

        # Set machine
        try:
            self.machine = Machine()
        except:
            raise
        # Then the scheduler should be already fixed here
        self.sched = self.machine.get_scheduler()

        # Set queue which should come after setting machine
        self.machine.set_queue(queue)
        self.queue = queue

        #...Check if other clmgr is already running,
        #...and if so, stop with error message.
        if self.already_running():
            msg = """clmgr is already running.
Only one clmgr can run at a time.
Please wait until the other clmgr stops or stop it manually.
        """
            raise RuntimeError(msg)

        #...Make pid file to show one clmgr is running with current pid.
        mypid = os.getpid()
        with open(self._pid_file,'w') as f:
            f.write('{0:d}\n'.format(mypid))


    def already_running(self):
        """
        Check whether the clmgr is already running.
        Only one clmgr can run at the same time.
        """
    
        if os.path.exists(self._pid_file):
            with open(self._pid_file,'r') as f:
                file_pid = int(f.readline())
            this_pid = os.getpid()
            if this_pid != file_pid:
                if pid_exists(file_pid):
                    return True
                else:
                    return False
        else:
            return False

    def cleanup(self):
        os.remove(self._pid_file)

    def set_calculator(self,calculator):
        if not calculator:
            raise RuntimeError('No calculator is set.')
        self.calculator = calculator
        if calculator in ('vasp',):
            import nappy.vasp
            self.calc_module = nappy.vasp
            self.Calculator = nappy.vasp.VASP
        else:
            self.logger.debug('? calculator = '+calculator)
            raise ValueError('calculator is not defined, calculator = '+calculator)

    
    def find_dirs_to_work(self,basedirs=None):
        """
        Find directories in which DFT calculations should be done.
        """
        if basedirs is None:
            basedirs = ['.',]

        #...Look for dirs that include POSCAR file.
        #...Among these dirs, look for dirs that include INCAR, KPOINTS, POTCAR as well
        self.dirs_to_work = []
        for basedir in basedirs:
            for path in find(basedir+'/','POSCAR'):
                # logger.debug("find: path = {}".format(path))
                d = os.path.dirname(os.path.abspath(path))
                calc = self.Calculator(d)
                if calc.needs_calc():
                    self.dirs_to_work.append(d)

    def avoid_conflict(self):
        """
        Check if there are already previous process waiting for the DFT calc
        in the given directories.
        And remove directories that are in progress by other clmgrs.
        """
        scheduler = self.machine.scheduler
        if scheduler in ('fx',):
            from nappy.scheduler.fujitsu import get_jobids
            
        elif scheduler in ('torque', 'pbs'):
            from nappy.scheduler.pbs import get_jobids
        elif scheduler in ('sge',):
            from nappy.scheduler.sge import get_jobids
        else:
            raise RuntimeError("No scheduler for {}".format(scheduler))
        
        #...Get jobids currently running on the cluster
        jobids = get_jobids()

        dirs = copy.copy(self.dirs_to_work)
        self.dirs_to_work = []
        for d in dirs:
            if not os.path.exists(d+'/'+self._stat_file):
                self.dirs_to_work.append(d)
                continue
            with open(d+'/'+self._stat_file,'r') as f:
                stat = int(f.readline())
            #...Skip this directory because some job is/will be perfromed here
            if stat in jobids:
                continue
            #...Otherwise, put this directory into new_dirs 
            self.dirs_to_work.append(d)

    def single_job_per_submission(self,dryrun=False):
        """
        Assign all the jobs to corresponding jobscripts in appropriate directories.
        """

        # Prepare mpi command keys.
        # There are some common keys between all the machines
        # such as `npara`, `out` and `exec_path`.
        # And there could be some keys very specific to some machines.
        mpi_command_keys = self.machine.mpi_command_keys()
        mpi_command_dict = {}
        for k in mpi_command_keys:
            if k == 'npara':
                mpi_command_dict[k] = None
            elif k == 'out':
                mpi_command_dict[k] = 'out.'+self.calculator
            elif k == 'exec_path':
                mpi_command_dict[k] = self.calc_module.get_exec_path()
            elif k == 'rankfile':
                mpi_command_dict[k] = None
        
        cwd = os.getcwd()
        pid = os.getpid()
        njobs = 0
        jobs = []
        npn = self.machine.nprocs_per_node
        for i,d in enumerate(self.dirs_to_work):
            os.chdir(d)
            calc = self.Calculator(d)
            nnodes,npn1,npara = calc.estimate_nprocs(max_npn=npn)
            nprocs = nnodes *npn1
            ctime = calc.estimate_calctime(nprocs=nprocs)
            if ctime > self.machine.qattr['limit_sec']:
                logger.info('The estimated calctime at {0:s} seems to '.format(d)
                            +'be longer than the limit_sec:')
                logger.info('  estimated_calctime = {0:d}'.format(ctime))
                logger.info('  limited time       = {0:d}'.format(self.machine.qattr['limit_sec']))
            job_info = {}  # initialize job_info
            job_info['JOB_NAME'] = 'clmgr{0:d}_{1:d}'.format(pid,i)
            job_info['QUEUE'] = queue
            job_info['NNODES'] = nnodes
            job_info['NPROCS_NODE'] = npn1
            job_info['WORKDIR'] = d
            hours = int(ctime / 3600 + 1)
            job_info['WALLTIME'] = '{0:d}:00:00'.format(hours)
            mpi_command_dict['npara'] = npara
            command = self.machine.get_mpi_command(**mpi_command_dict)
            job_info['COMMANDS'] = command
            script = self.sched.script_single(job_info)
            with open(self._batch_fname,'w') as f:
                f.write(script)
    
            if dryrun:
                logger.info(self.sched.get_command('submit')
                            +' '+self._batch_fname)
                jobid = 0
                jobs.append((d,jobid))
            else:
                try:
                    jobid = self.sched.submit(self._batch_fname)
                    njobs += 1
                    jobs.append((d,jobid))
                except Exception, e:
                    logger.warning('There is an error occurred, e = ',e)
                    pass
        os.chdir(cwd)
        return jobs        

    def plural_jobs_per_submission(self,dryrun=False):
        pass

class Machine(object):

    _machine_file = _clmgr_dir +'/machine'
    _default_mpi_command = "mpirun -np {npara} {exec_path} > {out} 2>&1"

    _sample_machine_file = """
scheduler: pbs
queues:
    'fx-debug':  {num_nodes:  32,  default_sec:  3600, limit_sec:   3600}
    'fx-small':  {num_nodes:  16,  default_sec: 86400, limit_sec: 604800}
    'fx-middle': {num_nodes:  96,  default_sec: 86400, limit_sec: 259200}
    'fx-large':  {num_nodes: 192,  default_sec: 86400, limit_sec: 259200}
    'fx-xlarge': {num_nodes: 864,  default_sec: 86400, limit_sec:  86400}
nprocs_per_node: 12
"""

    def __init__(self):
        self.scheduler = None
        self.qattr = None
        self.mpi_command = None
        self.machine_conf = None

        # Check existence of the machine configuration file
        if not os.path.exists(self._machine_file):
            logger.warn("Please create a machine-attribute file at "+_clmgr_dir)
            logger.warn("which is like this (in YAML format):")
            logger.warn(self._sample_machine_file)
            raise RuntimeError('No '+self._machine_file)

        self.load()

        
    def load(self):
        import yaml
        
        #...Read and check scheduler and queue from _machine_file
        with open(self._machine_file,'r') as f:
            self.machine_conf = yaml.load(f)
        
        if not self.machine_conf.has_key('scheduler') \
           or len(self.machine_conf['scheduler']) < 1:
            print('No valid scheduler in '+self._machine_file)
            raise RuntimeError()
        if not self.machine_conf.has_key('queues') \
           or len(self.machine_conf['queues']) < 1:
            print('No valid queues in '+self._machine_file)
            raise RuntimeError()
        if not self.machine_conf.has_key('nprocs_per_node'):
            print('No valid nprocs_per_node in '+self._machine_file)
            raise RuntimeError()

        # Set default MPI command unless something is specified
        if not self.machine_conf.has_key('mpi_command'):
            self.machine_conf['mpi_command'] = self._default_mpi_command
        
        self.scheduler = self.machine_conf['scheduler']
        self.mpi_command = self.machine_conf['mpi_command']
        self.nprocs_per_node = self.machine_conf['nprocs_per_node']

    def set_queue(self,queue):
        if not self.machine_conf['queues'].has_key(queue):
            print('There is no '+queue+' in queues from '+self._machine_file)
            raise RuntimeError()
        self.qattr = self.machine_conf['queues'][queue]

    def get_scheduler(self):
        if self.scheduler in ('pbs','torque'):
            import nappy.scheduler.pbs as sched
        elif self.scheduler in ('fx',):
            import nappy.scheduler.fujitsu as sched
        elif self.scheduler in ('sge',):
            import nappy.scheduler.sge as sched
        else:
            raise NotImplementedError()
        return sched

    def get_mpi_command(self,**kwarg):
        """
        The `mpi_command` in the machine_conf should be like,
        ::
        
          mpi_command: "mpirun -np {npara} {exec_path} > {out}"
          or
          mpi_command: "mpiexec -n {npara} --stdout {out} {exec_path}"
          or
          mpi_command: "mpiexec -n {npara} --vcoordfile {rankfile} --stdout {out} {exec_path}"

        The arguments should correspond to this mpi_command.
        """
        return self.machine_conf['mpi_command'].format(**kwarg)

    def check_mpi_command(self,**kwarg):
        keys = self.mpi_command_keys()
        kwarg_keys = kwarg.keys()
        if not keys == kwarg_keys:
            raise ValueError('Inconsistent keys between mpi_command and kwarg.')

    def mpi_command_keys(self,):
        """
        Parse mpi_command in the machine_conf and return keys in the command.
        """
        import re
        command = self.machine_conf['mpi_command']
        pattern = r"\{(\w+)\}"
        keys = re.findall(pattern,command)
        return keys

    
if __name__ == "__main__":

    start = time.time()

    args = docopt(__doc__)
    calculator = args['-c']
    dry = bool(args['-d'])
    jobs_per_submission = int(args['-j'])
    queue = args['-q']
    
    dirs = args['DIRS']
    if isinstance(dirs,str):
        dirs = [dirs,]
    for i in range(len(dirs)):
        dirs[i] = os.path.abspath(dirs[i])

    logger = getLogger(__name__)
    fHandler = FileHandler(ClusterManager._log_file,'w')
    fHandler.setLevel(DEBUG)
    sHandler = StreamHandler()
    sHandler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(fHandler)
    logger.addHandler(sHandler)
        
    clmgr = ClusterManager(calculator=calculator,queue=queue)
    clmgr.find_dirs_to_work(dirs)
    clmgr.avoid_conflict()
    if jobs_per_submission == 1:
        jobs = clmgr.single_job_per_submission(dryrun=dry)
    else:
        jobs = clmgr.plural_jobs_per_submission(dryrun=dry)

    logger.info("")
    logger.info("Directories treated in this clmgr {0:d} ".format(os.getpid())
                +"and corresponding Job IDs.")
    for job in jobs:
        d = job[0]
        jid = job[1]
        logger.info("{0:s}  {1:d}".format(d,jid))
    logger.info("")
    clmgr.cleanup()

    end = time.time()
    logger.info("Finished correctly with elapsed time = {0:10.2f} sec".format(end-start))

