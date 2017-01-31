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
  -d          Run clmgr without actually submitting jobs. [default: False]
  -n NPN
              Number of processers per node. [default: 12]
  -j NUM_JOBS
              Number of jobs in one submission.
              In case that there is a limitation to the number of jobs 
              that can run at the same time in the server, you may 
              had better set this larger than 1 to reduce the number of 
              submission by increasing number of jobs in one submission.
              [default: 1]
  -m MAX_NNODES_PER_SUBMISSION
              Maximum number of nodes per submission that can be used in the
              queue specified by the following QUEUE_NAME.
              [default: 16]
  -q QUEUE_NAME
              Queue name the jobs are submitted to. [default: default]
"""
from __future__ import print_function

import os
from docopt import docopt
from subprocess import Popen,PIPE
from datetime import datetime
import time
from logging import getLogger, StreamHandler, FileHandler, DEBUG, WARNING

import nappy
from nappy.clutil.find import find

__author__ = "RYO KOBAYASHI"
__version__ = "170122"

_stat_file = '.clmgr_stat'
_clmgr_dir = nappy.get_nappy_dir()+'/clmgr'
_pid_file = _clmgr_dir +'/pid'
_machine_file = _clmgr_dir +'/machine'
_date = datetime.now().strftime("%y%m%d")
_log_file = _clmgr_dir +'/log_{0:s}_{1:d}'.format(_date,os.getpid())
_batch_fname = 'batch.clmgr.sh'

def already_running():
    """
    Check whether the clmgr is already running.
    Only one clmgr can run at the same time.
    """

    if os.path.exists(_pid_file):
        with open(_pid_file,'r') as f:
            file_pid = int(f.readline())
        this_pid = os.getpid()
        if this_pid != file_pid:
            if pid_exists(file_pid):
                return True
            else:
                return False
    else:
        return False

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


def find_dirs_to_work(calculator,basedirs=None):
    """
    Find directories in which DFT calculations should be done.
    """
    logger = getLogger(__name__)
    
    if calculator in ('vasp',):
        from nappy.vasp import VASP
        Calculator = VASP
    else:
        raise ValueError('calculator is not defined, calculator =',
                         calculator)
    
    if basedirs is None:
        basedirs = ['.',]

    #...Look for dirs that include POSCAR file.
    #...Among these dirs, look for dirs that include INCAR, KPOINTS, POTCAR as well
    dirs_to_work = []
    for basedir in basedirs:
        for path in find(basedir+'/','POSCAR'):
            # logger.debug("find: path = {}".format(path))
            d = os.path.dirname(os.path.abspath(path))
            calc = Calculator(d)
            if calc.needs_calc():
                dirs_to_work.append(d)
    return dirs_to_work

def avoid_conflict(dirs,scheduler):
    """
    Check if there are already previous process waiting for the DFT calc
    in the given directories.
    And remove directories that are in progress by other clmgrs.
    """
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
    
    new_dirs = []
    for d in dirs:
        if not os.path.exists(d+'/'+_stat_file):
            new_dirs.append(d)
            continue
        with open(d+'/'+_stat_file,'r') as f:
            stat = int(f.readline())
        #...Skip this directory because some job is/will be perfromed here
        if stat in jobids:
            continue
        #...Otherwise, put this directory into new_dirs 
        new_dirs.append(d)

    return new_dirs



def one_job_one_submission(scheduler,calculator,calc_dirs,
                           queue,qattr,npn,
                           dryrun=False,logger=None):
    """
    Assign all the jobs to corresponding jobscripts in appropriate directories.
    """
    logger = logger or getLogger(__name__)
    
    if scheduler in ('pbs',):
        import nappy.scheduler.pbs as sched
    else:
        raise NotImplementedError()

    if calculator in ('vasp',):
        from nappy.vasp import VASP
        Calculator = VASP
    else:
        raise ValueError('calculator is not defined, calculator =',
                         calculator)
    cwd = os.getcwd()
    pid = os.getpid()
    njobs = 0
    jobs = []
    for i,d in enumerate(calc_dirs):
        os.chdir(d)
        calc = Calculator(d)
        nnodes,npn1,npara = calc.estimate_nprocs(max_npn=npn)
        nprocs = nnodes *npn1
        ctime = calc.estimate_calctime(nprocs=nprocs)
        if ctime > qattr['limit_sec']:
            logger.info('The estimated calctime at {0:s} seems to '.format(d)
                        +'be longer than the limit_sec:')
            logger.info('  estimated_calctime = {0:d}'.format(ctime))
            logger.info('  limited time       = {0:d}'.format(qattr['limit_sec']))
        job_info = {}  # initialize job_info
        job_info['JOB_NAME'] = 'clmgr{0:d}_{1:d}'.format(pid,i)
        job_info['QUEUE'] = queue
        job_info['NNODES'] = nnodes
        job_info['NPROCS_NODE'] = npn1
        job_info['WORKDIR'] = d
        hours = int(ctime / 3600 + 1)
        job_info['WALLTIME'] = '{0:d}:00:00'.format(hours)
        job_info['COMMANDS'] = calc.command_text().format(NPARA=npara)
        script = sched.script_for_one_job(job_info)
        with open(_batch_fname,'w') as f:
            f.write(script)

        if dryrun:
            logger.info(sched.get_command('submit')+' '+_batch_fname)
            jobid = 0
            jobs.append((d,jobid))
        else:
            try:
                jobid = sched.submit(_batch_fname)
                njobs += 1
                jobs.append((d,jobid))
            except Exception, e:
                logger.warning('There is an error occurred, e = ',e)
                pass
    os.chdir(cwd)
    return jobs

def initialize():
    """
    Initializer of the clmgr.
    """
    nappydir = nappy.get_nappy_dir()
    if not os.path.exists(nappydir):
        os.mkdir(nappydir)

    if not os.path.exists(_clmgr_dir):
        os.mkdir(_clmgr_dir)

    logger = getLogger(__name__)
    fHandler = FileHandler(_log_file,'w')
    fHandler.setLevel(DEBUG)
    sHandler = StreamHandler()
    sHandler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(fHandler)
    logger.addHandler(sHandler)

    if not os.path.exists(_machine_file):
        logger.warn("Please create a machine-attribute file at "+_clmgr_dir)
        logger.warn("which is like this (in YAML format):")
        msg = """
scheduler : pbs
queues :
    'fx-debug': { num_nodes : 32,  default_sec : 3600, limit_sec :  3600}
    'fx-small': { num_nodes : 16,  default_sec :86400, limit_sec :604800}
    'fx-middle':{ num_nodes : 96,  default_sec :86400, limit_sec :259200}
    'fx-large': { num_nodes :192,  default_sec :86400, limit_sec :259200}
    'fx-xlarge':{ num_nodes :864,  default_sec :86400, limit_sec : 86400}
"""
        logger.warn(msg)
        raise RuntimeError('No '+_machine_file)
        

    #...Check if other clmgr is already running,
    #...and if so, stop with error message.
    if already_running():
        msg = """clmgr is already running.
Only one clmgr can run at a time.
Please wait until the other clmgr stops or stop it manually.
        """
        raise RuntimeError(msg)

    #...Make pid file to show one clmgr is running with current pid.
    mypid = os.getpid()
    with open(_pid_file,'w') as f:
        f.write('{0:d}\n'.format(mypid))

    logger.info("initialize done")

def finalize(logger=None):
    logger = logger or getLogger(__name__)
    os.remove(_pid_file)
    logger.info("finalize done")

def main(basedirs,calculator,npn,jobs_per_submission,queue,dry,logger=None):
    """
    Search under given directories and check which directories
    are waiting for calculations, and then assign nodes and processes
    according to the estimated work load of these calculations.
    """
    import yaml

    logger = logger or getLogger(__name__)
    if calculator in ('vasp',):
        from nappy.vasp import VASP
        Calculator = VASP
    else:
        logger.debug('? calculator = '+calculator)
        raise ValueError('calculator is not defined, calculator = '+calculator)

    #...Read and check scheduler and queue from _machine_file
    with open(_machine_file,'r') as f:
        mdata = yaml.load(f)
    if not mdata.has_key('scheduler') or len(mdata['scheduler']) < 1:
        print('No valid scheduler in '+_machine_file)
        raise RuntimeError()
    if not mdata.has_key('queues') or len(mdata['queues']) < 1:
        print('No valid queues in '+_machine_file)
        raise RuntimeError()
    if not mdata['queues'].has_key(queue):
        print('There is no '+queue+' in queues from '+_machine_file)
        raise RuntimeError()

    scheduler = mdata['scheduler']
    qattr = mdata['queues'][queue]
    
    dirs0 = find_dirs_to_work(calculator,basedirs)
    calc_dirs = avoid_conflict(dirs0,scheduler)

    if jobs_per_submission == 1:
        jobs = one_job_one_submission(scheduler,calculator,calc_dirs,
                                      queue,qattr,npn,dryrun=dry)
    else:
        pass

    logger.info("")
    logger.info("Directories treated in this clmgr {0:d} ".format(os.getpid())
                +"and corresponding Job IDs.")
    for job in jobs:
        d = job[0]
        jid = job[1]
        logger.info("{0:s}  {1:d}".format(d,jid))
    logger.info("")
    

if __name__ == "__main__":

    start = time.time()

    args = docopt(__doc__)
    calculator = args['-c']
    dry = bool(args['-d'])
    npn = int(args['-n'])
    jobs_per_submission = int(args['-j'])
    queue = args['-q']
    
    dirs = args['DIRS']
    if isinstance(dirs,str):
        dirs = [dirs,]
    for i in range(len(dirs)):
        dirs[i] = os.path.abspath(dirs[i])

    initialize()
    logger = getLogger(__name__)
    main(dirs,calculator,npn,jobs_per_submission,queue,dry)
    finalize()

    end = time.time()
    logger.info("Finished correctly with elapsed time = {0:10.2f} sec".format(end-start))
