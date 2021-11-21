"""
PMD module.
"""
from __future__ import print_function

import os
import sys
import copy
import nappy
import numpy as np
try:
    import nappy.pmd.pmd_wrapper as pw
except:
    pass


__all__ = ['inpmd','pmd2phonopy','pairlist']

from . import *

__author__ = "RYO KOBAYASHI"
__version__ = "210930"

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
        if not ('nappy.pmd.pmd_wrapper' in sys.modules or 'pw' in sys.modules):
            raise ImportError('pmd_wrapper is not loaded.\n'
                              +'Probably you need to compile it at nap/nappy/pmd/.')
        self.params = nappy.pmd.inpmd.get_default()
        self.params['naux'] = 0
        if nsys is not None:
            self.nsys = copy.deepcopy(nsys)
            self.params['specorder'] = nsys.specorder
        else:
            self.nsys = None
        return None

    def run(self, nstp=0, dt=1.0, ifdmp=0, dmp=0.99, conveps=1e-5, convnum=3,
            initialize=True ):
        """
        Call pmd and store the result to self.result.
        """
        if self.nsys == None:
            raise ValueError('nsys must be set beofre calling run().')
        self.set_params(num_iteration=nstp, time_interval=dt, flag_damping=ifdmp,
                        damping_coeff=dmp, converge_eps=conveps, converge_num=convnum)
        # self.params['num_iteration'] = nstp
        # self.params['time_interval'] = dt
        self.update_params()
        self.update_mpivars()
        rtot = self.nsys.get_scaled_positions()
        vtot = np.zeros(rtot.shape)
        naux = self.params['naux']
        hmat = np.zeros((3,3,2))
        hmat[0:3,0:3,0] = self.nsys.get_hmat()
        ispcs = self.nsys.atoms.sid.values
        res = pw.run(rtot.T,vtot.T,naux,hmat,ispcs,initialize)
        self.result = {}
        self.result['rtot'] = res[0]
        self.result['vtot'] = res[1]
        self.result['atot'] = res[2]
        self.result['stot'] = res[3]
        self.result['ekitot'] = res[4]
        self.result['epitot'] = res[5]
        self.result['auxtot'] = res[6]
        self.result['hmat'] = res[7]
        self.result['ekin'] = res[8]
        self.result['epot'] = res[9]
        self.result['stnsr'] = res[10]
        return None

    def set_params(self,**kwargs):
        """
        Set parameters for pmdvars-module stored in self.params dictionary.
        """
        unknown_keys = []
        for k,v in kwargs.items():
            if k in self.params:
                self.params[k] = v
            else:
                unknown_keys.append(k)
        if len(unknown_keys) > 0:
            print(' Some of given parameters are unknown: ',unknown_keys)
        return None

    def update_params(self):
        """
        Update pmdvars-module variables.
        """
        
        keys = self.params.keys()
        if 'specorder' in keys:
            nspmax = 9
            cspcs = np.empty((nspmax,3),dtype='c')
            cspcs[:] = 'x  '
            specorder = self.params['specorder']
            for i in range(len(specorder)):
                cspcs[i] = str2char(specorder[i],3)
        else:
            raise KeyError('PMD().params has no specorder key.')
        if 'force_type' in keys:
            forces = self.params['force_type']
            nfrcs = len(forces)
            cfrcs = np.empty((nfrcs,128),dtype='c')
            for i in range(nfrcs):
                cfrcs[i] = str2char(forces[i],128)

        #...Extract variables from self.params with initial guesses
        rc = self.param2var('cutoff_radius',6.0)
        rbuf = self.param2var('cutoff_buffer',0.3)
        iprint = self.param2var('print_level',0)
        nstp = self.param2var('num_iteration',0)
        dt = self.param2var('time_interval',1.0)
        naux = self.param2var('naux',0)
        ifdmp = self.param2var('flag_damping',0)
        dmpcoeff = self.param2var('damping_coeff',0.99)
        conveps = self.param2var('converge_eps',1e-5)
        convnum = self.param2var('converge_num',3)
        tinit = self.param2var('initial_temperature',300.0)
        # tfin  = self.param2var('final_temperature',-10.0)
        # tcontrol = self.param2var('temperature_control','none')
        # ttgt = self.param2var('temperature_target',[300.0,100.0])
        # trlx = self.param2var('temperature_relax_time',50.0)
        sctrl = self.param2var('stress_control','vc-Berendsen')
        ptgt = self.param2var('pressure_target',0.0)
        stgt = np.array(self.param2var('stress_target',
                                       [[0.1, 0.0, 0.0],
                                        [0.0, 0.1, 0.0],
                                        [0.0, 0.0, 0.1]]))
        srlx = self.param2var('stress_relax_time',50.0)
        ifpmd = self.param2var('flag_out_pmd',2)
        npmd = self.param2var('num_out_pmd',0)
        nerg = self.param2var('num_out_energy',100)
        nnmax = self.param2var('max_num_neighbors',200)

        if 'Coulomb' in self.params['force_type']:
            naux = max(naux,2)
            self.params['naux'] = naux

        laux = 6
        cauxarr = np.empty((naux,laux),dtype='c')
        cauxarr[:] = '      '
        if naux >= 2:
            cauxarr[0] = 'chg   '
            cauxarr[1] = 'chi   '

        cpctrl = np.empty(20,dtype='c')
        cpctrl = str2char(sctrl,20)

        pw.set_pmdvars(cspcs,cfrcs,rc,rbuf,iprint,nstp,dt,cauxarr,
                       ifdmp,dmpcoeff,conveps,convnum,
                       cpctrl,ptgt,stgt.T,srlx,
                       ifpmd,npmd,nerg,nnmax)
        return None

    def param2var(self,s,v0):
        """
        Convert parameter name to its variable if exists in self.params, otherwise use v0 as the initial value.
        """
        if s in self.params.keys():
            return self.params[s]
        else:
            return v0

    def update_mpivars(self):
        """
        Set MPI-related variables in the pmdvars-module.
        """
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        fcomm = comm.py2f()
        size = comm.Get_size()
        rank = comm.Get_rank()
        #print('fcomm,size,rank=',fcomm,size,rank)
        pw.set_mpivars(fcomm,size,rank)
        return None

    def load_inpmd(self,):
        if not os.path.exists('in.pmd'):
            raise FileNotFoundError('in.pmd does not exist.')
        from .inpmd import read_inpmd
        inputs = read_inpmd('in.pmd')
        self.params = copy.deepcopy(inputs)
        return None
        
    def set_system(self,nsys):
        self.nsys = copy.deepcopy(nsys)
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
        return self.result['ekin']

    def get_potential_energy(self):
        if not hasattr(self,'result'):
            return None
        return self.result['epot']

    def get_stress_tensor(self):
        if not hasattr(self,'result'):
            return None
        return self.result['stnsr']

    def get_system(self):
        if not hasattr(self,'result'):
            raise ValueError('No MD result.')
        hmat = self.result['hmat']
        self.nsys.set_hmat(hmat[:,:,0])
        self.nsys.set_scaled_positions(self.result['rtot'].T)
        self.nsys.set_scaled_velocities(self.result['vtot'].T)
        self.nsys.atoms['epi'] = [ x for x in self.result['epitot'] ]
        self.nsys.set_kinetic_energy(self.result['ekin'])
        self.nsys.set_potential_energy(self.result['epot'])
        self.nsys.set_stress_tensor(self.result['stnsr'])
        return self.nsys
