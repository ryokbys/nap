"""
PMD module.
"""
from __future__ import print_function

import sys
try:
    import nappy.pmd.pmd_wrapper as pw
except:
    pass
import nappy
import numpy as np

__all__ = ['inpmd','pmd2phonopy']

from . import *

__author__ = "RYO KOBAYASHI"
__version__ = ""

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
        self.pmdparams = nappy.pmd.inpmd.get_default()
        self.pmdparams['naux'] = 0
        if nsys is not None:
            self.nsys = nsys
            self.pmdparams['specorder'] = nsys.specorder
        else:
            self.nsys = None
        return None

    def set_system(self,nsys):
        self.nsys = nsys
        self.pmdparams['specorder'] = nsys.specorder
        return None

    def run(self, nstp=0, dt=1.0, ):
        """
        Run pmd.
        """
        if self.nsys == None:
            raise ValueError('nsys must be set beofre calling run().')
        self.pmdparams['num_iteration'] = nstp
        self.pmdparams['time_interval'] = dt
        self.set_params()
        rtot = self.nsys.get_scaled_positions()
        vtot = np.zeros(rtot.shape)
        naux = self.pmdparams['naux']
        hmat = np.zeros((3,3,2))
        hmat[0:3,0:3,0] = self.nsys.get_hmat()
        ispcs = self.nsys.atoms.sid.values
        self.res = pw.run(rtot.T,vtot.T,naux,hmat,ispcs)
        return None

    def set_params(self):
        """
        Set pmdvars-module variables.
        """
        naux = 0
        iprint = 5
        rc = 6.0
        nstp = 0
        laux = 6
        
        keys = self.pmdparams.keys()
        if 'specorder' in keys:
            nspmax = 9
            cspcs = np.empty((nspmax,3),dtype='c')
            cspcs[:] = 'x  '
            specorder = self.pmdparams['specorder']
            for i in range(len(specorder)):
                cspcs[i] = str2char(specorder[i],3)
        if 'force_type' in keys:
            forces = self.pmdparams['force_type']
            nfrcs = len(forces)
            cfrcs = np.empty((nfrcs,128),dtype='c')
            for i in range(nfrcs):
                cfrcs[i] = str2char(forces[i],128)
        if 'cutoff_radius' in keys:
            rc = self.pmdparams['cutoff_radius']
        if 'print_level' in keys:
            iprint = self.pmdparams['print_level']
        if 'naux' in keys:
            naux = self.pmdparams['naux']
        if 'num_iteration' in keys:
            nstp = self.pmdparams['num_iteration']
        if 'Coulomb' in self.pmdparams['force_type']:
            naux = max(naux,2)
            self.pmdparams['naux'] = naux
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

    def set_potential(self,potential):
        if type(potential) is str:
            potential = potential.split()
        if type(potential) not in (list,tuple):
            raise TypeError('potential should be given as string or list of strings.')
        if len(potential) == 0:
            raise ValueError('potential should be given.')
        forces = self.pmdparams['force_type']
        if forces is None:
            forces = []
        for p in potential:
            if p not in forces:
                forces.append(p)
        self.pmdparams['force_type'] = forces
        return None
