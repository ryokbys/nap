# Cluster Manager utility

## What's `clmgr.py`?

This python utility automates the process of submitting jobs on clusters and supercomputers such as

* to make job scripts
* to submit jobs with specifying appropriate resources
* to combine some jobs to submit in one job script to reduce the actual number of submissions


### What is not provided by `clmgr.py`?

* to make input files for calculations
* to send input files from local computer to remote clusters
* to gather results from remote clusters
* to parse and check calculation results

There are some very sophisticated programs that can do these processes as well as the processes that `clmgr.py` provides, such as AiiDA by MARVEL and fireworks by materials project, but to remember how to use these very sophisticated program is also not so easy. Compared to those programs, `clmgr.py` is thought to be much easier to setup and use.


### Why is it included in nappy?

Basically this `clmgr.py` can be provided as a separated module, but it requires some parsers such as VASP parser which is already implemented in nappy and it is just easier for me to make this by utilizing the existing package. And also this will make installation process easier as well.


## Install nappy

```
export PYTHONPATH=/path/to/nap:$PYTHONPATH
```

Then you can test if the installation went well by
```
$ python -c "import nappy; print nappy.__file__"
```
If you get the path to nappy directory, the installation is OK.

## Setup on remote cluster machine

Following files are required to run `/path/to/nappy/clutil/clmgr.py`,
probably on a remote cluster or supercomputer.

* `~/.nappy/clmgr`
* `~/.nappy/clmgr/machine`
* `~/.nappy/vasp.conf`


### `machine` file

The machine file should be written in YAML format and contain information about:

* scheduler-type to be used
* list of queues and corresponding resources
* number of processors/cores in one node
* MPI command to be used

An example of `machine` file for Fujitsu FX100 is as follows,
```
scheduler: fujitsu
queues:
    'fx-debug':  {'num_nodes':  32, 'default_sec':  3600, 'limit_sec':   3600}
    'fx-small':  {'num_nodes':  16, 'default_sec': 86400, 'limit_sec': 604800}
    'fx-middle': {'num_nodes':  96, 'default_sec': 86400, 'limit_sec': 259200}
    'fx-large':  {'num_nodes': 192, 'default_sec': 86400, 'limit_sec': 259200}
    'fx-xlarge': {'num_nodes': 864, 'default_sec': 86400, 'limit_sec':  86400}
nprocs_per_node: 32
mpi_command: "mpiexec --vcoordfile {rankfile} -n {npara} --stdout {out} {exec_path}"
```

#### scheduler

`scheduler` can take values either one of the following,

* `fujitsu`
* `pbs` or `torque`


#### queue

`queues` should be a list of queue names and each queue name should have a dictionary with entries, `num_nodes`, `default_sec`, and `limit_sec`.

#### nprocs_per_node

This value is very machine specific and you must provide it.

#### mpi_command

The `mpi_command` entry specifies how the MPI command is used in the machine.
Currently following parameters can be used in the command,

* `{npara}`: number of parallel processes for the MPI run
* `{out}`: output file path
* `{exec_path}`: executable file path
* `{rankfile}`: path to the rankfile which is used in Fujitsu FX100


### `vasp.conf` file

This file should be written in YAML format and should contains the path to the *VASP* executable file like,

```
exec_path: '~/bin/vasp541fx'
```


## Usage

The usage can be shown by typing `clmgr.py -h`,
```
Usage:
  clmgr.py [options] DIRS [DIRS..]

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
```

So, for example, in Linux cluster without any limitation of running jobs per user, one can run `clmgr.py` like,
```
$ /path/to/clmgr.py -q queue_name calc_dir
```

Or in some supercomputer which limits the number of running jobs at a time,
```
$ /path/to/clmgr.py -q queue_name -m calc_dir
```

### log file

`clmgr.py` creates a log file at `~/.nappy/clmgr/` with the name including *PID* of the process that `clmgr.py` was running. In the file, you can find which directories are treated as directories where jobs will be running and to which job id the jobs are assigned.

