# -*- coding: utf-8 -*-
"""
Functions for torque/pbs scheduler.
"""
from __future__ import division
import re

__author__  = "Ryo KOBAYASHI"
__version__ = "170123"
__email__   = "kobayashi.ryo@nitech.ac.jp"

_fields = (
    'Job Id', 'Job_Name', 'Job_Owner', 'job_state',
    'queue', 'server', 'Checkpoint', 'ctime', 
    'Error_Path', 'exec_host', 'Hold_Types'
    'Join_Path', 'Keep_Files', 'Mail_Points'
    'mtime', 'Output_Path', 'Priority', 'qtime', 'Rerunable',
    'Resource_List.cput', 'Resource_List.nodect', 
    'Resource_List.nodes', 'session_id',
    'Variable_List', 'etime'
    'submit_args', 'start_time', 'start_count',
)

_commands = {
    'submit' : 'qsub',
    'delete' : 'qdel',
    'status' : 'qstat',
    'full-status' : 'qstat -f',
}

_script_template_single = """#!/bin/bash
#PBS -N {JOB_NAME}
#PBS -o out
#PBS -q {QUEUE}
#PBS -j oe
#PBS -l nodes={NNODES}:ppn={NPROCS_NODE}
#PBS -l walltime={WALLTIME}

cd {WORKDIR}

echo 'started at ' `date`
{COMMANDS}
echo 'ended at ' `date`
"""

def get_command(command_type):
    if _commands.has_key(command_type):
        return _commands[command_type]
    else:
        return None

def get_jobids():
    from subprocess import Popen, PIPE

    cmd = get_command('full-status')
    p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE)
    command_out,err = p.communicate()
    jobdata = parse_jobdata(command_out)

    jobids = []
    for job in jobdata:
        jobids.append(int(job['Job Id']))
    return jobids

def parse_jobdata(command_out):
    """
    Parse job data from output of `qstat -f` command.
    `command_out` should be passed to this function as argument.
    """

    output = [ o for o in command_out.splitlines() ]

    jobdata = []
    reading_job = False
    job = {}
    for line in output:
        if 'Job Id' in line:
            if len(job) != 0:
                jobdata.append(job)
                job = {}  ## refresh job
            jobid_str = line.split()[2]
            job['Job Id'] = jobid_str.split('.')[0]
        else:
            for fld in _fields:
                if fld in line:
                    entry = line.split()
                    value = entry[1]
                    if value == '=':
                        value = entry[2]
                    job[fld] = value
    if len(job) > 0:
        jobdata.append(job)
    return jobdata
    

def script_single(job_info):
    """
    Create job script context using given `job_info`.
    The `job_info` has to contain the following keys:

    - JOB_NAME
    - QUEUE
    - NNODES
    - NPROCS_NODE
    - WORKDIR
    - COMMANDS

    This only makes a script for only one calculation.
    The job script that can run more than 1 calculation is now under construction...
    """

    return _script_template_single.format(**job_info)

def submit(script_path):
    from subprocess import Popen,PIPE

    command = '{0} {1}'.format(get_command('submit'),script_path)
    try:
        p = Popen(command,shell=True,stdout=PIPE,stderr=PIPE)
        out,err = p.communicate()
        jobid = int(out.split('.')[0])
    except:
        raise
        
    return jobid

def delete(jobid):
    from subprocess import Popen, PIPE

    command = '{0} {1:d}'.format(get_command('delete'),jobid)
    try:
        p = Popen(command,shell=True,stdout=PIPE,stderr=PIPE)
        out,err = p.communicate()
    except:
        raise
    return out


