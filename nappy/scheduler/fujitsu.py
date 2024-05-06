# -*- coding: utf-8 -*-
"""
The scheduler used in Fujitsu PRIMEHPC FX100.
Since it seems to be no specific name for the scheduler, actually the scheduler
is a part of the "Technical Computing Suite (TCS)", this module was named as fujitsu.py.
The schedulers used in FX10 and K computer are supposed to be the same as this,
but it has not been tested for those two systems.
"""
import re
from logging import getLogger

__author__  = "Ryo KOBAYASHI"
__version__ = "200417"
__email__   = "kobayashi.ryo@nitech.ac.jp"

# This maps Fujitsu-FX state codes to our own status list

## List of states from the man page of pjstat
## ACC: Accepted job submission
## RJT: Rejected job submission
## QUE: Waiting for job execution
## SIN: During stage-in [FX10]
## RDY: Waiting for beginning of job execution [FX10]
## RNA: Acquiring resources required for job execution
## RNP: Executing prologue
## RUN: Executing job
## RNE: Executing epilogue
## RNO: Waiting for completion of job termination processing
## SOT: During stage-out [FX10]
## SWO: Swap-out in progress [FX100]
## SWD: Swap-out done [FX100]
## SWI: Swap-in in progress [FX100]
## EXT: Exited job end execution
## CCL: Exited job execution by interruption
## HLD: In fixed state due to users
## ERR: In fixed state due to an error

# _map_status_fx = {
#     'ACC':job_states.QUEUED,   ##  Accepted job submission
#     'RJT':job_states.DONE,     ##  Rejected job submission
#     'QUE':job_states.QUEUED,   ##  Waiting for job execution
#     'SIN':job_states.RUNNING,  ##  During stage-in [FX10]
#     'RDY':job_states.QUEUED,   ##  Waiting for beginning of job execution [FX10]
#     'RNA':job_states.QUEUED,   ##  Acquiring resources required for job execution
#     'RNP':job_states.RUNNING,  ##  Executing prologue
#     'RUN':job_states.RUNNING,  ##  Executing job
#     'RNE':job_states.RUNNING,  ##  Executing epilogue
#     'RNO':job_states.RUNNING,  ##  Waiting for completion of job termination processing
#     'SOT':job_states.RUNNING,  ##  During stage-out [FX10]
#     'SWO':job_states.RUNNING,  ##  Swap-out in progress [FX100]
#     'SWD':job_states.RUNNING,  ##  Swap-out done [FX100]
#     'SWI':job_states.RUNNING,  ##  Swap-in in progress [FX100]
#     'EXT':job_states.DONE,     ##  Exited job end execution
#     'CCL':job_states.DONE,     ##  Exited job execution by interruption
#     'HLD':job_states.QUEUED_HELD,     ##  In fixed state due to users
#     'ERR':job_states.DONE,     ##  In fixed state due to an error
#     }

# A posssible line in Fujitsu-FX is:
#   pjsub: [INFO] PJM 0000 pjsub Job 388645 submitted.
_fx_submitted_regexp = re.compile(
    r'.*(\[INFO\] PJM)\s.*\s(pjsub Job)\s+(?P<jobid>\d+)\s+(submitted.)'
)

# FX time format
_fx_time_format = '%Y/%m/%d %H:%M:%S'

# Separator between fields in the output of squeue
_field_separator = ","

# Fields to query or to parse
# The name like job_id should not be changed like reserved names.
_fields = [
    ("jid",'job_id'),     # job id
    ("jnam",'job_name'),  # job name
    ("st",'state_raw'),   # state
    ("usr",'username'),   # username
    ("nnumr",'number_nodes'),  # number of nodes requested
    ("cnumr",'number_cpus'),   # number of cpus
    ("elpl",'time_limit'),     # time limit in "DD hh:mm:ss"
    ("elp",'time_used'),       # time used by the job in "DD hh:mm:ss"
    ("sdt",'dispatch_time'),   # actual or expected dispatch time (start time)
    ("adt",'submission_time'),
    ("cmt",'annotation'),  # command in FX
    ("rscg",'partition'),  # resource group in FX
    ("nidlu",'allocated_machines'),
]

# Fields that are returned by pjstat command
_pjstat_fields = (
    'JOB ID',
    'JOB NAME',
    'STATE',
    'USER',
    'NODE NUM',
    'CPU NUM',
    'ELAPSE TIME (LIMIT)',
    'ELAPSE TIME (USE)',
    'RESOURCE GROUP',
)

_commands = {
    'submit': 'pjsub',
    'delete': 'pjdel',
    'status': 'pjstat',
}

# Total number of cores in a node, 32 for FX100
_num_cores_in_node = 32

_script_template_single = """#!/bin/bash
#PJM -N {JOB_NAME}
#PJM -L "rscgrp={QUEUE}"
#PJM -L "node={NNODES}"
#PJM --mpi "proc={NPROCS}"
#PJM --mpi "rank-map-bychip"
#PJM -L "elapse={WALLTIME}"
#PJM -o out_{JOB_NAME}
#PJM -j
#PJM -X
#PJM -S

cd {WORKDIR}

echo "started at `date`"
{COMMANDS}
echo "ended at `date`"
"""

#...Script template for plural jobs per one submission, where only one cycle of
#...plural jobs will be done and wait all the jobs done.
_script_template_plural = """#!/bin/bash
#PJM -N {JOB_NAME}
#PJM -L "rscgrp={QUEUE}"
#PJM -L "node={NNODES}"
#PJM --mpi "proc={NPROCS}"
#PJM --mpi "rank-map-bychip"
#PJM -L "elapse={WALLTIME}"
#PJM -o out_{JOB_NAME}
#PJM -j
#PJM -X
#PJM -S

export OMP_NUM_THREADS=1
export PARALLEL=1

echo "started at `date`"
cwd=`pwd`
{COMMANDS}
wait
echo "ended at `date`"
"""

def get_command(command_type):
    if command_type in _commands:
        return _commands[command_type]
    else:
        return None

def get_joblist_command():
    command = ["pjstat","-S",
               "--choose {}".format(_field_separator.join(
                   f[0] for f in _fields)) ]
    comm = ' '.join(command)
    return comm
    
def parse_jobdata(command_out):
    """
    Parse job data from output of `pjstat -S` command.
    `command_out` should be obtained like following,
    .. code:: 

       from subprocess import Popen, PIPE
       cmd = get_joblist_command()
       p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE)
       command_out,err = p.communicate()
    
    """

    output = [ o for o in command_out.splitlines() ]

    jobdata = []
    reading_job = False
    job = {}
    for line in output:
        if type(line) == bytes:
            line = line.decode('utf-8')
        if ':' not in line:  # job entry not started
            if reading_job:
                jobdata.append(job)
            reading_job = False
            job = {}  # reset job
        if 'JOB ID' in line:  # new job entry starts
            reading_job = True
            job = {}
        for field in _pjstat_fields:
            if field in line:
                job[field] = line.split(':')[1].strip()
        #jobdata.append(':'.join(line.split(':')[1:]).strip())
    if reading_job:
        jobdata.append(job)
        reading_job = False
    
    return jobdata

def get_jobids():
    from subprocess import Popen, PIPE
    
    cmd = get_joblist_command()
    p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE)
    command_out,err = p.communicate()
    jobdata = parse_jobdata(command_out.decode('utf-8'))

    jobids = []
    for job in jobdata:
        jobids.append(int(job['JOB ID']))
    return jobids

def script_single(job_info,template=None):
    """
    Create job script context using given `job_info`.
    The `job_info` has to contain the following keys:

    - JOB_NAME
    - QUEUE
    - NNODES
    - NPROCS
    - WORKDIR
    - COMMANDS

    This makes a script for only one calculation.
    """

    if template is not None:
        script = template
    else:
        script = _script_template_single

    return script.format(**job_info)

def script_plural(job_info,template=None):
    """
    Create job script context using given `job_info`.
    The `job_info` has to contain the following keys:

    - JOB_NAME
    - QUEUE
    - NNODES
    - NPROCS
    - WORKDIR
    - COMMANDS

    This makes a script for plural calculations.
    """

    if template is not None:
        script = template
    else:
        script = _script_template_plural

    return script.format(**job_info)


def submit(script_path,logger=None):
    """
    Submit a job with a given `script_path` and return its `JOB ID`.
    Output from pjsub command is something like,
    ::

      [INFO] PJM 0000 pjsub Job 543330 submitted.

    In case of error, something like,
    ::

      [ERR.] PJM 0020 pjsub Cannot open the file (rscgrp=fx-debug) : 2.

    """
    from subprocess import Popen,PIPE

    logger = logger or getLogger(__name__)

    command = '{0} {1}'.format('pjsub',script_path)
    try:
        p = Popen(command,shell=True,stdout=PIPE,stderr=PIPE)
        out,err = p.communicate()
        if type(out) == bytes:
            out = out.decode('utf-8')
        data = out.split(' ')
        if 'ERR' in data[0]:
            logger.warn('Error in submission, {0}'.format(script_path))
            #raise ValueError('Error in submission.')
        else:
            # JOB ID is at the 5th position of output (counting from 0)
            jobid = int(data[5])
    except Exception as e:
        raise
        
    return jobid
