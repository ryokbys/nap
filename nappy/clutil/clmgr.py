#!/usr/bin/env python
"""
Manage DFT calculations on the cluster computer via some scheduler.
- Search directories under given directory whether there are directories
  in which a DFT calculation should be done.
- Compute nodes required to compute those DFT calculations.
- Create job scripts for each submission.

Usage:
  clmgr.py [options] DIRS...

Options:
  -h, --help  Show this message and exit.
  -c CALCULATOR
              Set calculator type. Available types are: vasp
              [default: vasp]
  -d          Dryrun: run clmgr without actually submitting jobs.
              [default: False]
  -m          Enable multiple jobs in one submission.
              In case that there is a limitation to the number of jobs 
              that can run at the same time in the server, you may 
              had better set this TRUE.
              [default: False]
  -q QUEUE_NAME
              Queue name the jobs are submitted to. [default: default]
  --limit-sec LIMIT_SEC
              Maximum time of computation in second which will replace the
              machine default value. Negative value means to use the machine default.
              [default: -1]
  --template-script PATH
              Path to the template script in which some keywords are replaced by the scheduler.
              [default: None]
"""
import os
from docopt import docopt
from subprocess import Popen,PIPE
from datetime import datetime
import time
import copy
from logging import getLogger, StreamHandler, FileHandler, DEBUG

import nappy
from nappy.clutil.find import find

__author__ = "RYO KOBAYASHI"
__version__ = "200417"

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

def make_rankfile(offset,nnodes,npn,d,logger=None):
    """
    Create rankfile, which is required in Fujitsu FX100 when plural jobs are 
    to be running in one submission, into the directory specified by `d`.
    The rankfile should be like,
    ::

      (0)
      (0)
      ...
      (1)
      (1)
      ...
    
    The number of lines should be the same as the number of processes run by MPI
    for the job. And the `N` in `(N)` is the node index of the submission
    starting from 0.

    Parameters:

    - offset:  Offset of node index `N`.
    - nnodes:  Number of nodes to be used by the job.
    - npn:     Number of processes per node to be used.
    - d:       Directory where the `rankfile` is to be saved.
    """
    logger = logger or getLogger(__name__)
    with open(d+'/rankfile','w') as f:
        for n in range(nnodes):
            for i in range(npn):
                f.write('({0:d})\n'.format(offset+n))
    return None

def sec2hms(sec):
    """
    Convert seconds to hours, minutes and seconds.
    """
    hours = int(sec/3600)
    minutes = int((sec -3600*hours)/60)
    seconds = int(sec -3600*hours -60*minutes)
    return hours,minutes,seconds

class ClusterManager(object):
    
    _stat_file = '.clmgr_stat'
    _pid_file = _clmgr_dir +'/pid'
    _date = datetime.now().strftime("%y%m%d")
    _log_file = _clmgr_dir +'/log_{0:s}_{1:d}'.format(_date,os.getpid())
    _batch_fname = 'batch.clmgr{batch_id}.sh'
    
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
            self.machine = Machine(logger=self.logger)
        except Exception as e:
            raise
        # The scheduler should be already fixed here
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
        elif calculator in ('pmd',):
            import nappy.pmd
            self.calc_module = nappy.pmd
            self.Calculator = nappy.pmd.PMD
        else:
            self.logger.debug('? calculator = '+calculator)
            raise ValueError('calculator is not defined, calculator = '+calculator)

    
    def find_dirs_to_work(self,basedirs=None):
        """
        Find directories in which DFT calculations should be done.
        """
        if basedirs is None:
            basedirs = ['.',]

        if self.calculator in ('vasp',):
            keyfile = 'POSCAR'
        elif self.calculator in ('pmd',):
            keyfile = 'pmdini'

        #...Look for dirs that include the keyfile.
        #...Among these dirs, look for dirs that are ready to run the calculation.
        self.dirs_to_work = []
        for basedir in basedirs:
            for path in find(basedir+'/',keyfile):
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
        if scheduler in ('fx','fujitsu'):
            from nappy.scheduler.fujitsu import get_jobids
        elif scheduler in ('torque', 'pbs'):
            from nappy.scheduler.pbs import get_jobids
        elif scheduler in ('sge',):
            from nappy.scheduler.sge import get_jobids
        else:
            raise RuntimeError("No scheduler for {}".format(scheduler))
        
        #...Get jobids currently running on the cluster
        jobids = get_jobids()
        # self.logger.info('jobids:')
        # for jid in jobids:
        #     self.logger.info('{0:d}'.format(jid))

        dirs = copy.copy(self.dirs_to_work)
        self.dirs_to_work = []
        for d in dirs:
            if not os.path.exists(d+'/'+self._stat_file):
                self.dirs_to_work.append(d)
                continue
            with open(d+'/'+self._stat_file,'r') as f:
                try:
                    stat = int(f.readline())
                except Exception as e:
                    stat = -1
            #...Skip this directory because some job is/will be perfromed here
            if stat in jobids:
                continue
            #...Otherwise, put this directory into new_dirs
            self.dirs_to_work.append(d)

    def make_mpi_command_dict(self):
        """
        Prepare mpi command keys.
        There are some common keys between all the machines
        such as `npara`, `out` and `exec_path`.
        And there could be some keys very specific to some machines.
        """
        
        mpi_command_keys = self.machine.mpi_command_keys()
        self.mpi_command_dict = {}
        for k in mpi_command_keys:
            if k == 'npara':
                self.mpi_command_dict[k] = None
            elif k == 'out':
                self.mpi_command_dict[k] = 'out.'+self.calculator
            elif k == 'exec_path':
                self.mpi_command_dict[k] = self.calc_module.get_exec_path()
            elif k == 'rankfile':
                self.mpi_command_dict[k] = None


    def single_job_per_submission(self,dryrun=False,limit_seconds=-1,
                                  template_script=None):
        """
        Assign all the jobs to corresponding jobscripts in appropriate directories.
        """

        self.make_mpi_command_dict()
        
        limit_sec = self.machine.qattr['limit_sec']
        if limit_seconds > 0:
            limit_sec = min(limit_sec,limit_seconds)
            logger.info('Limit second is modified from machine default value'
                        +' to {0:d}.'.format(limit_sec))

        cwd = os.getcwd()
        pid = os.getpid()
        njobs = 0
        jobs = []
        npn = self.machine.nprocs_per_node
        for i,d in enumerate(self.dirs_to_work):
            os.chdir(d)
            calc = self.Calculator(d)
            nnodes,npn1,npara = calc.estimate_nprocs(max_npn=npn)
            #nprocs = nnodes *npn1
            ctime = calc.estimate_calctime(nprocs=npara)
            if ctime > limit_sec:
                logger.info('Since the estimated calctime {0:s}'.format(d)
                            +' seems to be longer than the limit_sec,'
                            +' the maximum calucation time of the queue is applied.')
                ctime = limit_sec
            # ctime = max(ctime,3600)
            job_info = {}  # initialize job_info
            job_info['JOB_NAME'] = 'clmgr{0:d}_{1:d}'.format(pid,i)
            job_info['QUEUE'] = queue
            job_info['NNODES'] = nnodes
            job_info['NPROCS_NODE'] = npn1
            job_info['WORKDIR'] = d
            hours,minutes,seconds = sec2hms(ctime)
            # hours = int(ctime / 3600 + 1)
            job_info['WALLTIME'] = '{0:d}:{1:02d}:{2:02d}'.format(hours,
                                                                  minutes,
                                                                  seconds)
            self.mpi_command_dict['npara'] = npara
            job_info['NPROCS'] = npara
            # if self.mpi_command_dict.has_key('rankfile'):
            #     make_rankfile(0,nnodes,npn1,d)
            #     self.mpi_command_dict['rankfile'] = './rankfile'
            command = self.machine.get_mpi_command(**self.mpi_command_dict)
            job_info['COMMANDS'] = command
            script = self.sched.script_single(job_info,template_script)
            batch_fname = self._batch_fname.format(batch_id=i)
            with open(batch_fname,'w') as f:
                f.write(script)
    
            if dryrun:
                logger.info(self.sched.get_command('submit')
                            +' '+batch_fname
                            +' at '+d)
                jobid = 0
                jobs.append((d,jobid))
            else:
                try:
                    jobid = self.sched.submit(batch_fname)
                    njobs += 1
                    jobs.append((d,jobid))
                except Exception as e:
                    print('Error: ',e)
                    logger.warning('There is an error occurred, e = ',e)
                    pass
            with open(self._stat_file,'w') as f:
                f.write('{0:d}\n'.format(jobid))
        os.chdir(cwd)
        return jobs        

    def plural_jobs_per_submission(self,dryrun=False,limit_seconds=-1,
                                   template_script=None):
        """
        Assign all the jobs by grouping some jobs to one submission.
        Run groupped submission at ~/.nappy/clmgr/tmpdir/ without specific
        output.
        """
        self.make_mpi_command_dict()

        tmpdir = _clmgr_dir+'/tmpdir'
        os.system('mkdir -p '+tmpdir)

        cwd = os.getcwd()
        pid = os.getpid()
        njobs = 0
        jobs = []
        #...Num procs per node
        npn = self.machine.nprocs_per_node
        #...Num nodes per submission
        limit_nodes = self.machine.qattr['num_nodes']
        limit_sec = self.machine.qattr['limit_sec']
        if limit_seconds > 0:
            limit_sec = min(limit_sec,limit_seconds)
            logger.info('Limit second is modified from machine default value'
                        +' to {0:d}.'.format(limit_sec))
        max_ctime = 0.0
        #...Initialize job_info
        job_info = {}
        #...Initial submission-ID
        isub = 1  
        job_info['WORKDIR'] = tmpdir
        job_info['QUEUE'] = queue
        #...Use full procs per node
        job_info['NPROCS_NODE'] = npn
        dirs = []
        sum_nodes = 0
        commands = ""
        for i,d in enumerate(self.dirs_to_work):
            os.chdir(d)
            calc = self.Calculator(d)
            nnodes,npn1,npara = calc.estimate_nprocs(max_npn=npn)
            if sum_nodes+nnodes > limit_nodes:
                #...Register job_info and dirs to jobs
                isub = len(jobs)+1
                job_info['NNODES'] = sum_nodes
                job_info['NPROCS'] = npn *sum_nodes
                job_info['JOB_NAME'] = 'clmgr{0:d}_{1:d}'.format(pid,isub)
                job_info['COMMANDS'] = commands
                # hours = int(max_ctime /3600 +1)
                hours,minutes,seconds = sec2hms(max_ctime)
                job_info['WALLTIME'] = '{0:d}:{1:02d}:{2:02d}'.format(hours,
                                                                      minutes,
                                                                      seconds)
                job = {}
                job['info'] = copy.copy(job_info)
                job['dirs'] = copy.copy(dirs)
                jobs.append(job)
                #...Initialize some
                sum_nodes = 0
                commands = ""
                dirs = []
                max_ctime = 0.0
            make_rankfile(sum_nodes,nnodes,npn1,d)
            sum_nodes += nnodes
            dirs.append(d)
            ctime = calc.estimate_calctime(nprocs=npara)
            max_ctime = min(max(max_ctime,ctime),limit_sec)
            #...max_ctime = max(max_ctime,3600)
            self.mpi_command_dict['npara'] = npara
            self.mpi_command_dict['rankfile'] = './rankfile'
            commands += "cd {0:s}\n".format(d) \
                        +self.machine.get_mpi_command(**self.mpi_command_dict) \
                        +" &\n"
        if dirs:
            # Register rest of works to jobs
            isub = len(jobs) +1
            job_info['JOB_NAME'] = 'clmgr{0:d}_{1:d}'.format(pid,isub)
            job_info['NNODES'] = sum_nodes
            job_info['NPROCS'] = npn *sum_nodes
            job_info['COMMANDS'] = commands
            # hours = max(int(max_ctime /3600),1)
            hours,minutes,seconds = sec2hms(max_ctime)
            # job_info['WALLTIME'] = '{0:d}:00:00'.format(hours)
            job_info['WALLTIME'] = '{0:d}:{1:02d}:{2:02d}'.format(hours,
                                                                  minutes,
                                                                  seconds)
            job = {}
            job['info'] = copy.copy(job_info)
            job['dirs'] = copy.copy(dirs)
            jobs.append(job)

        #...Submit jobs to the scheduler
        # print('len(jobs) = ',len(jobs))
        os.chdir(tmpdir)
        jobs2 = []
        for i, job in enumerate(jobs):
            jobnum = i+1
            job_info = job['info']
            dirs = job['dirs']
            script = self.sched.script_plural(job_info,template_script)
            batch_fname = self._batch_fname.format(batch_id=
                                                   '{0:d}_{1:d}'.format(pid,jobnum))
            with open(batch_fname,"w") as f:
                f.write(script)

            if dryrun:
                logger.info(self.sched.get_command('submit')
                            +' '+batch_fname+' at '+os.getcwd())
                jobid = 0
            else:
                try:
                    jobid = self.sched.submit(batch_fname)
                    njobs += 1
                except Exception as e:
                    logger.warn('There is an error occurred, e = ',e)
                    pass
            for d in dirs:
                jobs2.append((d,jobid))
                with open(d+'/'+self._stat_file,'w') as f:
                    f.write('{0:d}\n'.format(jobid))
        os.chdir(cwd)
        return jobs2


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
mpi_command: 'mpirun -np {npara} {exec_path} > {out} 2>&1'
"""

    def __init__(self,logger=None):
        self.scheduler = None
        self.qattr = None
        self.mpi_command = None
        self.machine_conf = None

        self.logger = logger or getLogger(__name__)
        
        self.load()

        
    def load(self):
        # Check existence of the machine configuration file
        if os.path.exists(self._machine_file+'.yaml'):
            import yaml
            try:
                with open(self._machine_file+'.yaml','r') as f:
                    self.machine_conf = yaml.safe_load(f)
            except Exception as e:
                raise
        elif os.path.exists(self._machine_file) or \
             os.path.exists(self._machine_file+'.json'):
            import json
            try:
                with open(self._machine_file,'r') as f:
                    self.machine_conf = json.load(f)
            except Exception as e:
                raise
        else:
            self.logger.warn("Please create a machine-attribute file 'machine'"
                             +" or 'machine.yaml' at "+_clmgr_dir)
            self.logger.warn("which is like this (in JSON format):")
            self.logger.warn(self._sample_machine_file)
            raise RuntimeError('No '+self._machine_file)

        
        if 'scheduler' not in self.machine_conf \
           or len(self.machine_conf['scheduler']) < 1:
            self.logger.warn('No valid scheduler in '+self._machine_file)
            raise RuntimeError()
        if 'queues' not in self.machine_conf \
           or len(self.machine_conf['queues']) < 1:
            self.logger.warn('No valid queues in '+self._machine_file)
            raise RuntimeError()
        if 'nprocs_per_node' not in self.machine_conf:
            self.logger.warn('No valid nprocs_per_node in '+self._machine_file)
            raise RuntimeError()

        # Set default MPI command unless something is specified
        if 'mpi_command' not in self.machine_conf:
            self.machine_conf['mpi_command'] = self._default_mpi_command
        
        self.scheduler = self.machine_conf['scheduler']
        self.mpi_command = self.machine_conf['mpi_command']
        self.nprocs_per_node = self.machine_conf['nprocs_per_node']

    def set_queue(self,queue):
        if queue not in self.machine_conf['queues']:
            self.logger.warn('There is no '+queue+' in queues from '+self._machine_file)
            raise RuntimeError()
        self.qattr = self.machine_conf['queues'][queue]

    def get_scheduler(self):
        if self.scheduler in ('pbs','torque'):
            import nappy.scheduler.pbs as sched
        elif self.scheduler in ('fx','fujitsu'):
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
        command = self.machine_conf['mpi_command'].format(**kwarg)
        if '--vcoordfile None' in command:
            command = command.replace('--vcoordfile None','')
        return command

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
    multiple_jobs_per_submission = bool(args['-m'])
    queue = args['-q']
    limit_sec = int(args['--limit-sec'])
    template_path = args['--template-script']

    if template_path == 'None':
        template_path = None
    if not os.path.exists(template_path):
        raise IOError('Template script does not exist at ',template_path)
    else:
        with open(template_path,'r') as f:
            lines = f.readlines()
            template_script = ""
            for l in lines:
                template_script += l
    
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
    if multiple_jobs_per_submission:
        jobs = clmgr.plural_jobs_per_submission(dryrun=dry,
                                                limit_seconds=limit_sec,
                                                template_script=template_script)
    else:
        jobs = clmgr.single_job_per_submission(dryrun=dry,
                                               limit_seconds=limit_sec,
                                               template_script=template_script)

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

