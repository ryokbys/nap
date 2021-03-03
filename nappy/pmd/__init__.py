"""
PMD module.
"""
from __future__ import print_function

import os
import sys
import copy
try:
    import nappy.pmd.pmd_wrapper as pw
except:
    pass
import nappy
import numpy as np

__all__ = ['inpmd','pmd2phonopy']

from . import *

__author__ = "RYO KOBAYASHI"
__version__ = "210302"

def str2char(string,nlen):
    slen = len(string)
    if slen > nlen:
        raise ValueError('slen > nlen !')
    c = np.empty(nlen,dtype='c')
    c[:] = ' '
    c[0:slen] = [ string[i] for i in range(slen) ]
    c[slen:] = [ ' ' for i in range(nlen-slen) ]
    return c

class PMD:

    def __init__(self, nsys=None):
        if 'nappy.pmd.pmd_wrapper' not in sys.modules:
            raise ImportError('pmd_wrapper is not loaded.\n'
                              +'Probably you need to compile it at nap/nappy/pmd/.')
        self.params = nappy.pmd.inpmd.get_default()
        self.params['naux'] = 0
        if nsys is not None:
            self.nsys = nsys
            self.params['specorder'] = nsys.specorder
        else:
            self.nsys = None
        return None

    def run(self, nstp=0, dt=1.0 ):
        """
        Run pmd.
        """
        if self.nsys == None:
            raise ValueError('nsys must be set beofre calling run().')
        self.params['num_iteration'] = nstp
        self.params['time_interval'] = dt
        self.update_params()
        rtot = self.nsys.get_scaled_positions()
        vtot = np.zeros(rtot.shape)
        naux = self.params['naux']
        hmat = np.zeros((3,3,2))
        hmat[0:3,0:3,0] = self.nsys.get_hmat()
        ispcs = self.nsys.atoms.sid.values
        self.result = pw.run(rtot.T,vtot.T,naux,hmat,ispcs)
        return None

    def update_params(self):
        """
        Update pmdvars-module variables.
        """
        naux = 0
        iprint = 0
        rc = 6.0
        nstp = 0
        laux = 6
        
        keys = self.params.keys()
        if 'specorder' in keys:
            nspmax = 9
            cspcs = np.empty((nspmax,3),dtype='c')
            cspcs[:] = 'x  '
            specorder = self.params['specorder']
            for i in range(len(specorder)):
                cspcs[i] = str2char(specorder[i],3)
        if 'force_type' in keys:
            forces = self.params['force_type']
            nfrcs = len(forces)
            cfrcs = np.empty((nfrcs,128),dtype='c')
            for i in range(nfrcs):
                cfrcs[i] = str2char(forces[i],128)
        if 'cutoff_radius' in keys:
            rc = self.params['cutoff_radius']
        if 'print_level' in keys:
            iprint = self.params['print_level']
        if 'naux' in keys:
            naux = self.params['naux']
        if 'num_iteration' in keys:
            nstp = self.params['num_iteration']
        if 'Coulomb' in self.params['force_type']:
            naux = max(naux,2)
            self.params['naux'] = naux
        cauxarr = np.empty((naux,laux),dtype='c')
        cauxarr[:] = '      '
        if naux >= 2:
            cauxarr[0] = 'chg   '
            cauxarr[1] = 'chi   '
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        fcomm = comm.py2f()
        size = comm.Get_size()
        rank = comm.Get_rank()
        pw.set_pmdvars(cspcs,cfrcs,rc,fcomm,rank,size,iprint,nstp,cauxarr)
        return None

    def load_inpmd(self,):
        if not os.path.exists('in.pmd'):
            raise FileNotFoundError('in.pmd does not exist.')
        from .inpmd import read_inpmd
        inputs = read_inpmd('in.pmd')
        self.params = copy.copy(inputs)
        return None
        
    def set_system(self,nsys):
        self.nsys = nsys
        self.params['specorder'] = nsys.specorder
        return None

    def set_potential(self, potential):
        if type(potential) is str:
            potential = potential.split()
        if type(potential) not in (list,tuple):
            raise TypeError('potential should be given as string or list of strings.')
        if len(potential) == 0:
            raise ValueError('potential should be given.')
        forces = self.params['force_type']
        if forces is None:
            forces = []
        for p in potential:
            if p not in forces:
                forces.append(p)
        self.params['force_type'] = forces
        return None

    def get_kinetic_energy(self):
        if not hasattr(self,'result'):
            return None
        return self.result[8]

    def get_potential_energy(self):
        if not hasattr(self,'result'):
            return None
        return self.result[9]

    def get_stress_tensor(self):
        if not hasattr(self,'result'):
            return None
        return self.result[10]

    def get_system(self):
        if not hasattr(self,'result'):
            raise ValueError('No MD result.')
        hmat = self.result[7]
        self.nsys.set_hmat(hmat[:,:,0])
        self.nsys.set_scaled_positions(self.result[0].T)
        self.nsys.set_scaled_velocities(self.result[1].T)
        self.nsys.atoms['epi'] = [ x for x in self.result[5] ]
        self.nsys.set_kinetic_energy(self.result[8])
        self.nsys.set_potential_energy(self.result[9])
        self.nsys.set_stress_tensor(self.result[10])
        return self.nsys
