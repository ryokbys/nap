"""
Python classes and functions useful for using DFT code VASP.
"""
from __future__ import print_function

import os
import copy
import math

from nappy import vasp

# def needs_calc(path):
#     """
#     Check whether the VASP calculation is needed on the given path.
#     """
#     import os

#     infiles = ('POSCAR','INCAR','POTCAR','KPOINTS')
#     #...check if necessary for VASP files exist
#     for f in infiles:
#         if not os.path.exists(path+'/'+f):
#             return False

#     #...if there is not OUTCAR in the directory, VASP should be performed
#     outcar = 'OUTCAR'
#     if not os.path.exists(path+'/'+outcar):
#         return True
    
#     #...compare dates between OUTCAR and infiles
#     #...to determine if the VASP calc should be performed
#     outcar_mtime = os.stat(path+'/'+outcar).st_mtime
#     for f in infiles:
#         infile_mtime = os.stat(path+'/'+f).st_mtime
#         if infile_mtime > outcar_mtime:
#             return True

#     #...otherwise no need to perform VASP
#     return False


class VASP(object):
    
    _input_files = (
        'POSCAR',
        'INCAR',
        'POTCAR',
        'KPOINTS'
    )

    _output_files = (
        'OUTCAR',
        'vasprun.xml'
    )

    def __init__(self, path):
        if not os.path.exists(path):
            raise RuntimeError('path does not exist.')

        self.path = path
        self.input_files = {}
        for f in self._input_files:
            if os.path.exists(self.path+'/'+f):
                self.input_files['f'] = os.stat(self.path+'/'+f).st_mtime
            else:
                self.input_files['f'] = None
        self.output_files = {}
        for f in self._output_files:
            if os.path.exists(self.path+'/'+f):
                self.output_files['f'] = os.stat(self.path+'/'+f).st_mtime
            else:
                self.output_files['f'] = None

        self.config = {}
        self.load_config()
        
    def needs_calc(self):
        """
        Check whether the VASP calc is needed.
        """
        #...check if necessary for VASP files exist
        for f,t in self.input_files.items():
            if not t:
                return False
        #...if there is not OUTCAR in the directory, VASP should be performed
        for f,t in self.output_files.items():
            if not t:
                return True
        
        #...compare dates between OUTCAR and infiles
        #...to determine if the VASP calc should be performed
        outcar_mtime = self.output_files['OUTCAR']
        for fin,tin in self.input_files:
            if tin > outcar_mtime:
                return True

        #...otherwise no need to perform VASP
        return False

    def get_num_valence(self):
        """
        Get number of valence electrons by reading POSCAR and POTCAR.
        Since POTCAR is large, reading POTCAR may take a while.
        """
        #...Read POSCAR file and get number of atoms of each species
        try:
            self.poscar
        except NameError:
            poscar = vasp.poscar.POSCAR(self.path+'/POSCAR')
        poscar.read()
        num_atoms = copy.copy(poscar.num_atoms)

        #...Read POTCAR file and get number of electrons of each species
        try:
            self.potcar
        except NameError:
            potcar = vasp.potcar.read_POTCAR(self.path+'/POTCAR')
        valences = potcar['valence']

        #...Total num of valence electrons
        nel = 0
        for i,n in enumerate(num_atoms):
            nel += n *valences[i]
        return nel

    
    def parse_KPOINTS(fname='KPOINTS'):
        """
        Parse KPOINTS file and get k-points.
        Note: The k-points read from file are not actually the number of 
        k-points used in the calculation, because it can be reduced by 
        using symmetry of the system considered.
        Estimation of the number of actual number of k-points is not easy 
        at the moment.
        """
        with open(fname,'r') as f:
            lines = f.readlines()

        comment = lines[0]
        auto = int(lines[1].split()[0])
        method = lines[2].split()[0]
        kpt = [ int(x) for x in lines[3].split()[0:3] ]
        if auto == 0:
            nkpt = kpt[0]*kpt[1]*kpt[2]
        else:
            nkpt = auto
        return nkpt

    
    def estimate_calctime(self,nprocs=1):
        """
        Estimate the computation time wrt num of valence electrons,
        num of bands, num of k-points and num of processes used.
        ::

          Time = 0.005 sec *sqrt(nprocs) /(nel *nband *nkpt)

        It is a very rough estimation of computation time.
        """
        try:
            self.incar
        except NameError:
            incar = vasp.incar.parse_INCAR()
            
        nel = self.get_num_valence()
        nkpt = self.parse_KPOINTS()
        
        nband = incar['NBAND']
        estime = 0.005 *math.sqrt(nprocs) /(nel *nband *nkpt)
        return estime

    
    def estimate_nprocs(self,procs_per_node=16):
        """
        Estimate the computation time wrt num of valence electrons,
        num of bands, num of k-points and procs_per_node.
        Make nprocs multiples of procs_per_node.
        """
        try:
            self.incar
        except NameError:
            self.incar = vasp.incar.parse_INCAR()

        nel = self.get_num_valence()
        nkpt = self.parse_KPOINTS()
        nband = self.incar['NBAND']
        if self.incar.has_key('NPAR') and self.incar.has_key('NCORE'):
            nprocs = self.incar['NPAR'] *self.incar['NCORE']
        elif self.incar.has_key('NPAR'):
            nprocs = self.incar['NPAR']
        elif self.incar.has_key('NCORE'):
            nprocs = self.incar['NCORE']
        #...Make nprocs to be multiples of procs_per_node
        nprocs = max(nprocs,procs_per_node)
        if nprocs % procs_per_node != 0:
            n = int(nprocs/procs_per_node)
            nprocs = (n+1) *procs_per_node

        return nprocs
        

    def command_text(self):
        """
        Make command text to run VASP using mpirun.
        """
        import nappy
        if not self.config.has_key('exec_path'):
            msg = """
Error: self.exec_path has not been set yet.
You should write a path to the VASP executable in ~/{0}/vasp.conf .
Like,
::

  exec_path  /home/username/bin/vasp535-openmpi


""".format(nappy._nappy_dir)
            raise RuntimeError(msg)
        
        text = 'mpirun -np {{txt}} {path}'.format(txt='NPROCS',
                                                  path=self.config['exec_path']) \
               +' > out.vasp 2>&1'
        return text
    
    def set_exec_path(self,path):
        self.config['exec_path'] = path
        return None

    def load_config(self):
        """
        Load config from `~/.nappy/vasp.conf` file.
        """
        import nappy
        with open(nappy._nappy_dir+'/vasp.conf','r') as f:
            lines = f.readlines()

        for line in lines:
            entry = line.split()
            self.config[entry[0]] = entry[1]
        return None
        
    def save_config(self):
        import nappy
        try:
            self.config
        except NameError:
            raise RuntimeError('self.config has not been set.')

        #...Read all the current config from file
        dic = {}
        with open(nappy._nappy_dir+'/vasp.conf','r') as f:
            lines = f.readlines()
        for line in lines:
            dat = line.split()
            dic[dat[0]] = dat[1]

        for k,v in self.config.items():
            dic[k] = v
        #...Override the config file
        with open(nappy._nappy_dir+'/vasp.conf','w') as f:
            for k,v in dic.items():
                f.write('{0:s}  {1:s}\n').format(k,v)
        return None
