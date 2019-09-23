#!/usr/bin/env python
"""
Cuckoo search.

Usage:
  cs.py [options]

Options:
  -h, --help  Show this message and exit.
  -n N        Number of generations in CS. [default: 20]
  --print-level LEVEL
              Print verbose level. [default: 1]
"""
from __future__ import print_function

import os,sys
from docopt import docopt
import numpy as np
from numpy import exp, sin, cos
import random
import copy
from multiprocessing import Process, Queue
from time import time
from scipy.special import gamma

__author__ = "RYO KOBAYASHI"
__version__ = "rev190920"

_fname_gen = 'out.cs.generations'
_fname_ind = 'out.cs.individuals'

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

def wrap(vs,vrs):
    vsnew = copy.copy(vs)
    for i,v in enumerate(vsnew):
        vmin, vmax = vrs[i]
        vsnew[i] = min(max(v,vmin),vmax)
    return vsnew

class Individual:
    """
    Individual class that consists of variables as vector elements.
    """
    def __init__(self, iid, ndim, vranges, loss_func):
        self.iid = iid
        self.ndim = ndim
        self.loss_func = loss_func
        self.vector = np.zeros(self.ndim)
        self.vranges = vranges
        self.val = None

    def set_variable(self,variables):
        if len(variables) != len(self.vector):
            raise ValueError()

        self.vector = variables
        # print('iid, v before wrap,vrs =',self.iid,self.vector,self.vranges)
        self.wrap_range()
        # print('iid, v after wrap,vrs  =',self.iid,self.vector,self.vranges)
        self.val = None
        return None

    def init_random(self):
        for i in range(self.ndim):
            vmin, vmax = self.vranges[i]
            # vmin = self.vranges[i,0]
            # vmax = self.vranges[i,1]
            v = random.random()*(vmax -vmin) +vmin
            self.vector[i] = v
            # print(' i,vmin,vmax,v=',i,vmin,vmax,v)
        self.wrap_range()
        self.val = None
        return None

    def wrap_range(self):
        self.vector = wrap(self.vector, self.vranges)

    def calc_loss_func(self,kwargs,q):
        """
        Compute loss function value using self.loss_func function given in the constructor.
        In order to return a result in multiprocessing.Process, it also takes an argument q.
        """
        # print('type(kwargs)=',type(kwargs))
        val = self.loss_func(self.vector, self.vranges, **kwargs)
        # print(' iid,v,val=',self.iid,self.vector,val)
        q.put(val)
        return None

class CS:
    """
    Cuckoo search class.
    """

    def __init__(self, N, F, variables, vranges, loss_func, write_func, **kwargs):
        """
        Conctructor of CS class.

        N:  Number of individuals.
        F:  Fraction of worse individuals to be abondoned.
        loss_func:
            Loss function to be minimized with variables and **kwargs.
        """
        if N < 2:
            raise ValueError('N must be greater than 1 in CS!')
        self.N = N   # Number of individuals in a generation
        self.F = F   # Fraction of worse individuals to be abondoned
        self.ndim = len(variables)
        self.vs = variables
        self.vrs = vranges
        self.vws = np.zeros(self.ndim)
        for i in range(self.ndim):
            self.vws[i] = max(self.vrs[i,1] -self.vrs[i,0], 0.0)
        # print('original variables=',self.vs,self.vrs)
        self.loss_func = loss_func
        self.write_func = write_func
        self.kwargs = kwargs
        self.bestind = None
        self.print_level = 0
        if 'print_level' in kwargs.keys():
            self.print_level = kwargs['print_level']

        self.beta = 1.5
        self.betai = 1.0 /self.beta
        self.usgm = (gamma(1+self.beta)*np.sin(np.pi*self.beta/2)/ \
                     gamma((1+self.beta)/2)*self.beta*2.0**((self.beta-1)/2))**self.betai
        self.vsgm = 1.0

        #...initialize population
        self.population = []
        self.iidmax = 0
        for i in range(N):
            self.iidmax += 1
            ind = Individual(self.iidmax, self.ndim, self.vrs, self.loss_func)
            if i == 0:
                ind.set_variable(self.vs)
            else:
                ind.init_random()
            self.population.append(ind)

        #...Evaluate loss function values
        qs = [ Queue() for i in range(self.N) ]
        prcs = []
        for ip,pi in enumerate(self.population):
            kwtmp = copy.copy(self.kwargs)
            kwtmp['index'] = ip
            kwtmp['iid'] = pi.iid
            prcs.append(Process(target=pi.calc_loss_func, args=(kwtmp,qs[ip])))
        for p in prcs:
            p.start()
        for p in prcs:
            p.join()
        for ip,pi in enumerate(self.population):
            pi.val = qs[ip].get()
            # print('ip,val,vec=',ip,pi.val,pi.vector)
        
        self.keep_best()
        if self.print_level > 2:
            for pi in self.population:
                self.write_variables(pi,
                                     fname='in.vars.fitpot.{0:d}'.format(pi.iid),
                                     **self.kwargs)
        else:
            self.write_variables(self.bestind,
                                 fname='in.vars.fitpot.{0:d}'.format(self.bestind.iid),
                                 **self.kwargs)
        return None

    def keep_best(self):
        vals = []
        for i,pi in enumerate(self.population):
            # print('i,val,vec=',i,pi.val,pi.vector)
            if pi.val == None:
                raise ValueError('Something went wrong.')
            vals.append(pi.val)

        minval = min(vals)
        if self.bestind == None or minval < self.bestind.val:
            idx = vals.index(minval)
            self.bestind = copy.deepcopy(self.population[idx])
        return None

    def sort_individuals(self):

        jtop = self.N
        for i in range(self.N):
            jtop -= 1
            for j in range(jtop):
                pj = self.population[j]
                pjp = self.population[j+1]
                if pj.val > pjp.val:
                    self.population[j] = pjp
                    self.population[j+1] = pj
        

    def run(self,maxiter=100):
        """
        Perfom CS.
        """

        if 'start' in self.kwargs.keys():
            start = self.kwargs['start']
        else:
            start = time()
        fgen = open(_fname_gen,'w')
        find = open(_fname_ind,'w')
        for i,ind in enumerate(self.population):
            fgen.write('     0  {0:8d}  {1:12.4e}\n'.format(ind.iid, ind.val))
            find.write(' {0:8d}  {1:12.4e}'.format(ind.iid, ind.val))
            for j,vj in enumerate(ind.vector):
                find.write(' {0:11.3e}'.format(vj))
            find.write('\n')

        if self.print_level > 0:
            print(' step,time,best,vars= {0:6d} {1:8.1f}  {2:8.4f}'.format(0, time()-start,
                                                                           self.bestind.val),end="")
            for i in range(min(16,self.ndim)):
                print(' {0:6.3f}'.format(self.bestind.vector[i]),end="")
            print('', flush=True)
            
        for it in range(maxiter):

            candidates = []

            self.sort_individuals()
            #...Create candidates by Levy flight
            vbest = self.bestind.vector
            for ip,pi in enumerate(self.population):
                vi = pi.vector
                vnew =np.array(vi)
                for iv in range(self.ndim):
                    u = np.random.normal()*self.usgm
                    v = abs(np.random.normal()*self.vsgm)
                    v = max(v,1.0e-8)
                    w = u/v**self.betai
                    zeta = self.vws[iv] *0.01 *w
                    # zeta = self.vws[iv]*0.01 *w *(vi[iv] -vbest[iv])
                    # if ip == 0:
                    #     zeta = self.vws[iv] *0.001 *w
                    # else:
                    #     zeta = 0.01 *w *(vi[iv] -vbest[iv])
                    vnew[iv] = vnew[iv] +zeta*np.random.normal()
                #...create new individual for trial
                # print('ip,vi,vnew=',ip,vi,vnew)
                self.iidmax += 1
                newind = Individual(self.iidmax, self.ndim, self.vrs, self.loss_func)
                newind.set_variable(vnew)
                candidates.append(newind)

            #...Evaluate loss function values
            qs = [ Queue() for i in range(len(candidates)) ]
            prcs = []
            for ic,ci in enumerate(candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = ic
                kwtmp['iid'] = ci.iid
                prcs.append(Process(target=ci.calc_loss_func, args=(kwtmp,qs[ic])))
            for p in prcs:
                p.start()
            for p in prcs:
                p.join()
            for ic,ci in enumerate(candidates):
                ci.val = qs[ic].get()

            #...Pick j that is to be compared with i
            js = random.sample(range(self.N),k=self.N)
            #...Decide whether or not to adopt new one
            for jc,jv in enumerate(js):
                pj = self.population[jv]
                cj = candidates[jc]
                dval = cj.val -pj.val
                if dval < 0.0: # replace with new individual
                    self.population[jv] = cj
                    find.write(' {0:8d}  {1:12.4e}'.format(cj.iid, cj.val))
                    for k,vk in enumerate(cj.vector):
                        find.write(' {0:11.3e}'.format(vk))
                    find.write('\n')
                else:
                    pass

            #...Rank individuals
            self.sort_individuals()
            
            #...Abandon bad ones and replace with random ones
            iab = int((1.0 -self.F)*self.N)
            candidates = []
            for iv in range(iab,self.N):
                self.iidmax += 1
                newind = Individual(self.iidmax, self.ndim, self.vrs, self.loss_func)
                newind.init_random()
                candidates.append(newind)
            
            #...Evaluate loss function values of new random ones
            qs = [ Queue() for i in range(len(candidates)) ]
            prcs = []
            for ic,ci in enumerate(candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = ic
                kwtmp['iid'] = ci.iid
                prcs.append(Process(target=ci.calc_loss_func, args=(kwtmp,qs[ic])))
            for p in prcs:
                p.start()
            for p in prcs:
                p.join()
            for ic,ci in enumerate(candidates):
                ci.val = qs[ic].get()

            #...Replace them with old ones
            ic = 0
            for iv in range(iab,self.N):
                ci = candidates[ic]
                ic += 1
                self.population[iv] = ci

            #...Check best
            for ic,ci in enumerate(self.population):
                if ci.val < self.bestind.val:
                    self.bestind = ci
                    self.write_variables(ci,
                                         fname='in.vars.fitpot.{0:d}'.format(ci.iid),
                                         **self.kwargs)
            
            if self.print_level > 0:
                print(' step,time,best,vars= {0:6d} {1:8.1f}  {2:8.4f}'.format(it+1, time()-start,
                                                                               self.bestind.val),end="")
                for i in range(min(16,self.ndim)):
                    print(' {0:6.3f}'.format(self.bestind.vector[i]),end="")
                print('', flush=True)

            for i,ind in enumerate(self.population):
                fgen.write(' {0:5d}  {1:8d}  {2:12.4e}\n'.format(it+1, ind.iid, ind.val))
        fgen.close()
        find.close()
        #...Finaly write out the best one
        self.write_variables(self.bestind,fname='in.vars.fitpot.best',**self.kwargs)
        return None

    def write_variables(self,ind,fname='in.vars.fitpot',**kwargs):
        
        vs = ind.vector
        vrs = ind.vranges
        self.write_func(vs,vrs,fname,**kwargs)
        return None


if __name__ == "__main__":

    args = docopt(__doc__)
    n = int(args['-n'])
    kwargs = {}
    kwargs['print_level'] = int(args['--print-level'])


    vs = np.array([1.0, -0.5])
    vrs = np.array([[-1.0, 2.0],[-1.0, 1.0]])
    
    cs = CS(10, 0.25, vs, vrs, test_func, test_write_func, **kwargs)
    cs.run(n)
    
