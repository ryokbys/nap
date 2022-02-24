#!/usr/bin/env python
"""
Tree-based Parzen Estimator (TPE).

Usage:
  {0:s} [options]

Options:
  -h, --help  Show this message and exit.
"""
from __future__ import print_function

import os, sys
from docopt import docopt
import numpy as np
import random
import copy
from multiprocessing import Process, Pool
from time import time

__author__ = "RYO KOBAYASHI"
__version__ = "rev220123"

def test_func(var, vranges, **kwargs):
    x,y= var
    res= x**2 +y**2 +100.0*exp(-x**2 -y**2)*sin(2.0*(x+y))*cos(2*(x-y)) \
         +80.0*exp(-(x-1)**2 -(y-1)**2)*cos(x+4*y)*sin(2*x-y) \
         +200.0*sin(x+y)*exp(-(x-3)**2-(y-1)**2)
    return res

def test_write_func(vs,vrs,fname,**kwargs):
    with open(fname,'w') as f:
        for i,v in enumerate(vs):
            vr = vrs[i]
            f.write(' {0:10.3f}  {1:10.3f}  {2:10.3f}\n'.format(v,*vr))
    return None

def random_variables(vranges):
    """
    Create random variables within the given ranges.
    """
    rvs = np.empty(len(vranges))
    for i in range(len(vranges)):
        vmin, vmax = vranges[i]
        rvs[i] = random.random()*(vmax -vmin) +vmin
        # print(' i,vmin,vmax,v=',i,vmin,vmax,v)
    return rvs

def gauss_kernel(x):
    return np.exp(-0.5*x*x)/np.sqrt(2.0*np.pi)

class Sample:
    """
    A sample point that has parameter vector and loss value.
    """
    def __init__(self, iid, ndim, loss_func):
        self.iid = iid
        self.ndim = ndim
        self.loss_func = loss_func
        self.variables = np.zeros(self.ndim)
        self.val = None
        return None

    def set_variables(self,variables):
        if len(variables) != len(self.variables):
            raise ValueError('len(variables) != len(self.variables)')

        self.variables[:] = variables[:]
        self.val = None
        return None

    def calc_loss_func(self, vlogs, **kwargs):
        """
        Compute loss function value using self.loss_func function given in the constructor.
        """
        vec = copy.copy(self.variables)
        # If vlog != None, some of variables may be expressed in log domain,
        # they must be transformed back to non-log domain.
        for i in range(len(vec)):
            if vlogs[i]:
                vec[i] = np.exp(vec[i])
        val = self.loss_func(vec, **kwargs)
        return val,kwargs['index']

class TPE:
    """
    Class for Tree-based Parzen Estimator (TPE) or Weighted Parzen Estimator (WPE).
    
    Refs:
      1. Bergstra, J., Bardenet, R., Bengio, Y. & Kégl, B.  in Proc. NIPS-24th, 2546–2554
    """

    def __init__(self, nbatch, variables, vranges, vlimits, loss_func,
                 write_func, **kwargs):
        """
        Conctructor of TPE class.

        Parameters:
          nbatch : int
            Number of samples in a batch which is the same as num of processes.
          variables : 1-D np.array
            np.array of variables to be optimized.
          vranges : 2-D np.array
            Current lower and upper limit of variables, which are used only at random sampling processes.
          vlimits : 2-D np.array
            Hard limits of variables, which are fixed during TPE iterations. 
          loss_func : function
            Loss function to be minimized with variables and **kwargs.
          write_func : function
            Function for outputing some info.
        """
        if nbatch < 1:
            raise ValueError('nbatch must be > 0.')
        self.nbatch = nbatch
        self.ndim = len(variables)
        self.vars0 = variables
        self.vranges = vranges
        self.vlimits = vlimits
        self.vlogs = [ False for i in range(len(self.vars0))]
        if 'vlogs' in kwargs.keys():
            self.vlogs = kwargs['vlogs']
        self.loss_func = loss_func
        self.write_func = write_func
        self.kwargs = kwargs
        self.best_pnt = None
        self.print_level = 0
        self.nsmpl_prior = 100
        self.ntrial = 100
        self.method = kwargs['fitting_method']
        self.fname_smpl = 'out.{0:s}.samples'.format(self.method)

        self.gamma = 0.15
        #...Change default values if specified
        if 'print_level' in kwargs.keys():
            self.print_level = int(kwargs['print_level'])
        if 'tpe_nsmpl_prior' in kwargs.keys():
            self.nsmpl_prior = int(kwargs['tpe_nsmpl_prior'])
        if 'tpe_ntrial' in kwargs.keys():
            self.ntrial = int(kwargs['tpe_ntrial'])
        if 'tpe_gamma' in kwargs.keys():
            self.gamma = float(kwargs['tpe_gamma'])

        if self.gamma < 0.0 or self.gamma > 1.0:
            raise ValueError('gamma must be within 0. and 1., whereas gamma = ',self.gamma)

        #...Write info
        print('')
        if self.method in ('wpe','WPE'):
            print('   {0:s} infomation:'.format(self.method))
            print(f'     Num of prior samples = {self.nsmpl_prior:d}')
            print(f'     Num of top samples used for density estimation = {self.ntrial:d}')
        elif self.method in ('tpe','TPE'):
            print('   {0:s} infomation:'.format(self.method))
            print(f'     Num of prior samples = {self.nsmpl_prior:d}')
            print(f'     Num of trials for sampling = {self.ntrial:d}')
            print(f'     Gamma for dividing high and low = {self.gamma:4.2f}')
        print('')
    
        #...Change vrange if log domain
        for i in range(self.ndim):
            if self.vlogs[i]:
                self.vars0[i] = np.log(self.vars0[i])
                for l in range(2):
                    self.vranges[i,l] = np.log(self.vranges[i,l])
                    self.vlimits[i,l] = np.log(self.vlimits[i,l])

        #...Initialize sample history
        self.history = []  # History of all samples
        self.iidmax = 0
        for i in range(self.nbatch):
            self.iidmax += 1
            smpl = Sample(self.iidmax, self.ndim, self.loss_func)
            if i == 0:
                smpl.set_variables(self.vars0)
            else:
                smpl.set_variables(random_variables(self.vranges))
            self.history.append(smpl)

        return None

    def keep_best(self):
        vals = []
        for i,si in enumerate(self.history):
            if si.val == None:
                raise ValueError('Something went wrong.')
            vals.append(si.val)

        minval = min(vals)
        if self.bestsmpl == None or minval < self.bestsmpl.val:
            idx = vals.index(minval)
            self.bestsmpl = copy.deepcopy(self.history[idx])
        return None

    def run(self,maxstp=0):
        """
        Perfom TPE.
        """

        starttime = time()

        fsmpl = open(self.fname_smpl,'w')
        #...Headers
        fsmpl.write('# {0:>7s}  {1:>12s}'.format('iid', 'loss'))
        for i in range(len(self.history[0].variables)):
            fsmpl.write(' {0:>8d}-th'.format(i+1))
        fsmpl.write('\n')

        #...Create pool before going into maxstp-loop,
        #...since creating pool inside could cause "Too many files" error.
        pool = Pool(processes=self.nbatch)

        #...Evaluate sample losses of initial sets
        prcs = []
        for i,ci in enumerate(self.history):
            kwtmp = copy.copy(self.kwargs)
            kwtmp['index'] = i
            kwtmp['iid'] = ci.iid
            prcs.append(pool.apply_async(ci.calc_loss_func, (self.vlogs,), kwtmp))

        results = [ res.get() for res in prcs ]
        for res in results:
            val, i = res
            self.history[i].val = val
        #self.history.extend(candidates) # not need
        
        #...Check best
        self.bestsmpl = self.history[0]
        for si in self.history[1:]:
            if si.val < self.bestsmpl.val:
                self.bestsmpl = si
        fname = 'in.vars.fitpot.{0:d}'.format(self.bestsmpl.iid)
        self.write_variables(self.bestsmpl,
                             fname=fname,
                             **self.kwargs)
        os.system('cp -f {0:s} in.vars.fitpot.best'.format(fname))
        
        #...Write sample data
        for i,si in enumerate(self.history):
            self._write_smpl_data(fsmpl,si)

        if self.print_level > 0:
            self._write_step_info(0,starttime)

        #...TPE loop starts
        for istp in range(1,maxstp):
            #...Create candidates by either random or TPE
            if len(self.history) <= self.nsmpl_prior:
                #...Create random candidates
                candidates = []
                for i in range(self.nbatch):
                    self.iidmax += 1
                    newsmpl = Sample(self.iidmax, self.ndim, self.loss_func)
                    newsmpl.set_variables(random_variables(self.vranges))
                    candidates.append(newsmpl)
            else:
                if self.method in ('wpe,''WPE'):
                    #...Create candidates by WPE
                    candidates = self._candidates_by_WPE()
                else:
                    #...Create candidates by TPE
                    candidates = self._candidates_by_TPE()
                    

            #...Evaluate sample losses of initial sets
            prcs = []
            for i,ci in enumerate(candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = i
                kwtmp['iid'] = ci.iid
                prcs.append(pool.apply_async(ci.calc_loss_func, (self.vlogs,), kwtmp))

            results = [ res.get() for res in prcs ]
            for res in results:
                val, i = res
                candidates[i].val = val
            self.history.extend(candidates)

            #...Check best
            best_updated = False
            for si in self.history[-self.nbatch:]:
                if si.val < self.bestsmpl.val:
                    self.bestsmpl = si
                    best_updated = True
            if best_updated:
                fname = 'in.vars.fitpot.{0:d}'.format(self.bestsmpl.iid)
                self.write_variables(self.bestsmpl,
                                     fname=fname,
                                     **self.kwargs)
                os.system('cp -f {0:s} in.vars.fitpot.best'.format(fname))

            for i,si in enumerate(self.history[-self.nbatch:]):
                self._write_smpl_data(fsmpl,si)

            #...Write info
            if self.print_level > 0:
                self._write_step_info(istp,starttime)
            
        fsmpl.close()
        pool.close()
        return None

    def _write_smpl_data(self,f,smpl):
        f.write(' {0:8d}  {1:12.4e}'.format(smpl.iid, smpl.val))
        for j,vj in enumerate(smpl.variables):
            f.write(' {0:11.3e}'.format(vj))
        f.write('\n')
        return None

    def _write_step_info(self,istp,starttime):
        print(' step,time,best,vars='
              +' {0:6d} {1:8.1f}  {2:8.4f}'.format(istp, time()-starttime,
                                                   self.bestsmpl.val),end="")
        for i in range(min(16,self.ndim)):
            print(' {0:6.3f}'.format(self.bestsmpl.variables[i]),end="")
        print('', flush=True)
        

    def _candidates_by_TPE(self,):
        """
        Create candidates by using TPE.
        """
        vals = np.zeros(len(self.history))
        for i,si in enumerate(self.history):
            vals[i] = si.val
        nlow = int(self.gamma *len(self.history))
        nhigh = len(self.history) -nlow
        argpart = np.argpartition(vals,nlow)
        Xlow = np.zeros((nlow,self.ndim))
        Xhigh = np.zeros((nhigh,self.ndim))
        for i in range(nlow):
            iid = argpart[i]
            Xlow[i,:] = self.history[iid].variables[:]
        for i in range(nhigh):
            iid = argpart[nlow+i]
            Xhigh[i,:] = self.history[iid].variables[:]

        #...Sampling variable candidates
        xcandidates = np.empty((self.nbatch,self.ndim))
        ntrial = max(self.ntrial,self.nbatch)
        for idim in range(self.ndim):
            xhmin = self.vlimits[idim,0]
            xhmax = self.vlimits[idim,1]
            xlowsrt = np.sort(Xlow[:,idim])
            npnt = len(xlowsrt)
            #...Determine smoothness parameter, h, by Silverman's method
            q75, q25 = np.percentile(xlowsrt, [75,25])
            std = np.std(xlowsrt)
            sgm = min(std, (q75-q25)/1.34)
            h = 1.06 *sgm /np.power(npnt, 1.0/5)
            #...Prepare for g(x)
            xhighsrt = Xhigh[:,idim]
            q75, q25 = np.percentile(xhighsrt, [75,25])
            sgmh = min(np.std(xhighsrt), (q75-q25)/1.34)
            hh = 1.06 *sgmh /np.power(len(xhighsrt),1.0/5)
            #...Several trials for selection
            aquisition = np.zeros(ntrial)
            xs = np.empty(ntrial)
            for itry in range(ntrial):
                ipnt = int(random.random()*npnt)
                xi = xlowsrt[ipnt]
                r = h *np.sqrt(-2.0*np.log(random.random()))
                th = 2.0 *np.pi *random.random()
                x = xi + r*np.cos(th)
                #...Wrap by vranges
                x = min(max(x,xhmin),xhmax)
                xs[itry] = x
                #...Compute l(x) and g(x)
                lx = 0.0
                for j in range(npnt):
                    z = (x-xlowsrt[j])/h
                    lx += np.exp(-0.5*z*z)
                lx /= npnt*h *np.sqrt(2.0*np.pi)
                gx = 0.0
                for j in range(len(xhighsrt)):
                    z = (x-xhighsrt[j])/hh
                    gx += np.exp(-0.5*z*z)
                gx /= len(xhighsrt)*hh *np.sqrt(2.0*np.pi)
                aquisition[itry] = gx/lx
            #...Pick nbatch of minimum aquisition points
            idxsort = np.argsort(aquisition)
            xcandidates[:,idim] = xs[idxsort[0:self.nbatch]]
        #...Create sample with xcandidate as variables
        candidates = []
        for ib in range(self.nbatch):
            self.iidmax += 1
            smpl = Sample(self.iidmax, self.ndim, self.loss_func)
            smpl.set_variables(xcandidates[ib,:])
            candidates.append(smpl)
        return candidates

    def _candidates_by_WPE(self,):
        """
        Create candidates by using WPE.
        """
        vals = np.zeros(len(self.history))
        for i,si in enumerate(self.history):
            vals[i] = si.val
        if len(self.history) > self.ntrial:
            iargs = np.argpartition(vals,self.ntrial)
            tmpsmpls = [ self.history[i] for i in iargs[:self.ntrial] ]
            vals = vals[ iargs[:self.ntrial] ]
        else:
            tmpsmpls = copy.copy(self.history)
        vmin = vals.min()
        wgts = np.array([ np.exp(-(v-vmin)/vmin) for v in vals ])
        xtmps = np.zeros((len(tmpsmpls),self.ndim))
        for i in range(len(tmpsmpls)):
            xtmps[i,:] = tmpsmpls[i].variables[:]
        #...Sampling variable candidates
        xcandidates = np.empty((self.nbatch,self.ndim))
        for idim in range(self.ndim):
            xhmin = self.vlimits[idim,0]
            xhmax = self.vlimits[idim,1]
            xsrt = np.sort(xtmps[:,idim])
            #...Determine smoothness parameter, h, by Silverman's method
            q75, q25 = np.percentile(xsrt, [75,25])
            sgm = min(np.std(xsrt), (q75-q25)/1.34)
            h = 1.06 *sgm /np.power(len(xsrt), 1.0/5)
            xs = np.empty(self.nbatch)
            for ib in range(self.nbatch):
                ipnt = random.choices([j for j in range(len(wgts))],weights=wgts)[0]
                xi = xtmps[ipnt,idim]
                r = h *np.sqrt(-2.0*np.log(random.random()))
                th = 2.0 *np.pi *random.random()
                x = xi + r*np.cos(th)
                #...Wrap by vranges
                x = min(max(x,xhmin),xhmax)
                xs[ib] = x
            xcandidates[:,idim] = xs[:]
        
        #...Create sample with xcandidate as variables
        candidates = []
        for ib in range(self.nbatch):
            self.iidmax += 1
            smpl = Sample(self.iidmax, self.ndim, self.loss_func)
            smpl.set_variables(xcandidates[ib,:])
            candidates.append(smpl)
        
        return candidates

    def write_variables(self,smpl,fname='in.vars.fitpot',**kwargs):
        vs = smpl.variables
        vrs = self.vranges
        self.write_func(vs,vrs,fname,**kwargs)
        return None

def main():
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])))
    return None

if __name__ == "__main__":

    main()
