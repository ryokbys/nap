"""
Python modules and classes for PMD.
"""
from __future__ import print_function

import os
import copy
import math

import nappy

import json


def get_conf_path():
    return nappy.get_nappy_dir()+'/pmd.conf'

def get_exec_path():
    conf_path = get_conf_path()
    with open(conf_path,'r') as f:
        config = json.load(f)
    if not 'exec_path' in config:
        msg = """
Error: self.exec_path has not been set yet.
You should write a path to the PMD executable in {0}.
It should be in JSON format like,
::

  {
    "exec_path":  "/home/username/bin/pmd"
  }


""".format(get_conf_path())
        raise RuntimeError('config does not have exec_path.')
    return config['exec_path']


class PMD:
    
    _input_files = (
        'in.pmd',
        'pmdini',
    )

    _output_files = (
        'out.pmd',
        'out.erg'
    )

    def __init__(self, path):
        if not os.path.exists(path):
            raise RuntimeError('path does not exist.')

        self.path = path
        self.input_files = {}
        for f in self._input_files:
            if os.path.exists(self.path+'/'+f):
                self.input_files[f] = os.stat(self.path+'/'+f).st_mtime
            else:
                self.input_files[f] = None
        self.output_files = {}
        for f in self._output_files:
            if os.path.exists(self.path+'/'+f):
                self.output_files[f] = os.stat(self.path+'/'+f).st_mtime
            else:
                self.output_files[f] = None

        self.config = {}
        self.load_config()

        
    def needs_calc(self):
        """
        Check whether the PMD calc is needed.
        """
        #...check if necessary for PMD files exist
        for f,t in self.input_files.items():
            if not t:
                return False
        #...if there is not OUTCAR in the directory, VASP should be performed
        for f,t in self.output_files.items():
            if not t:
                return True
        
        #...compare dates between out.pmd and infiles
        #...to determine if the PMD calc should be performed
        outpmd_mtime = self.output_files['out.pmd']
        for fin,tin in self.input_files.items():
            if tin > outpmd_mtime:
                return True

        #...Check whether the calculation done correctly.
        if not self.calc_done():
            return True

        #...otherwise no need to perform VASP
        return False

    def calc_done(self):
        if not os.path.exists('out.pmd'):
            return False
        with open('out.pmd','r') as f:
            lines = f.readlines()
        if 'finished' in lines[-1]:
            return True
        else:
            return False

    
    def estimate_calctime(self,nprocs=1):
        """
        Estimate the computation time.
        It is extremely hard to estimate comp time since it depends on
        which force-fields are used, how many atoms, how long is the MD...
        Just return 1 min = 60 sec.
        """
        return 60

    
    def estimate_nprocs(self,max_npn=16,limit_npn=None):
        """
        Estimate the nodes to be used.
        Currently, just return 1.
        """
        nnodes = 1
        npn = 1
        npara = 1

        return nnodes, npn, npara
        

    def get_exec_command(self):
        """
        Make command text to run PMD using mpirun.
        """
        if not 'exec_path' in self.config:
            msg = """
Error: self.exec_path has not been set yet.
You should write a path to the PMD executable in {0}.
It should be in JSON format like,
::

  {
    exec_path:  /home/username/bin/pmd
  }


""".format(get_conf_path())
            raise RuntimeError(msg)

        text = self.config['exec_path']
        # text = 'mpirun -np {{NPARA}} {path}'.format(path=self.config['exec_path']) \
        #         +' > out.pmd 2>&1'
        return text

    def set_exec_path(self,path):
        self.config['exec_path'] = path
        return None

    def get_exec_path(self):
        if not 'exec_path' in self.config:
            msg = """
Error: self.exec_path has not been set yet.
You should write a path to the PMD executable in {0}.
It should be in JSON format like,
::

  {
    "exec_path":  "/home/username/bin/pmd"
  }


""".format(get_conf_path())
            raise RuntimeError(msg)
        return self.config['exec_path']


    def load_config(self):
        """
        Load config from `~/.nappy/pmd.conf` file.
        """
        with open(get_conf_path(),'r') as f:
            self.config = json.load(f)

        return None
        
    def save_config(self):
        try:
            self.config
        except:
            raise RuntimeError('self.config has not been set.')

        #...Read all the current config from file
        dic = {}
        with open(get_conf_path(),'r') as f:
            lines = f.readlines()
        for line in lines:
            dat = line.split()
            dic[dat[0]] = dat[1]

        for k,v in self.config.items():
            dic[k] = v
        #...Override the config file
        with open(get_conf_path(),'w') as f:
            for k,v in dic.items():
                f.write('{0:s}  {1:s}\n').format(k,v)
        return None
