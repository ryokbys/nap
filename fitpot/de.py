#!/usr/bin/env python
"""
Differential evolution test program.

Usage:
  de.py [options]

Options:
  -h, --help  Show this message and exit.
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

__author__ = "RYO KOBAYASHI"
__version__ = "rev190904"

_fname_gen = 'out.de.generations'
_fname_ind = 'out.de.individuals'

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

class DE:
    """
    Differential evolution class.
    """

    def __init__(self, N, F, CR, T, variables, vranges, loss_func, write_func, **kwargs):
        """
        Conctructor of DE class.

        loss_func:
            Loss function to be minimized with variables and **kwargs.
        """
        if N < 4:
            raise ValueError('N must be greater than 3 in DE!')
        self.N = N   # Number of individuals in a generation
        self.F = F   # Fraction of mixing in DE
        self.CR = CR # Cross-over rate
        self.T = T   # Temperature (kT) to compute adoption probability
        if self.T < 1e-10:
            raise ValueError('T is too small.')
        self.ndim = len(variables)
        self.vs = variables
        self.vrs = vranges
        # print('original variables=',self.vs,self.vrs)
        self.loss_func = loss_func
        self.write_func = write_func
        self.kwargs = kwargs
        self.bestind = None
        self.print_level = 0
        if 'print_level' in kwargs.keys():
            self.print_level = kwargs['print_level']

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
                
    def run(self,maxiter=100):
        """
        Perfom DE.
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
            print(' step,time,best,vars= {0:8d} {1:8.1f}  {2:9.4f}'.format(0, time()-start,
                                                                           self.bestind.val),end="")
            for i in range(min(10,self.ndim)):
                print(' {0:11.3e}'.format(self.bestind.vector[i]),end="")
            print('', flush=True)

        
        for it in range(maxiter):

            candidates = []

            #...Create candidates
            for ip,pi in enumerate(self.population):
                vi = pi.vector
                #...pick other 3 individuals
                indices= [ j for j in range(self.N) if j != i ]
                irand = int(random.random()*len(indices))
                i1 = indices.pop(irand)
                irand = int(random.random()*len(indices))
                i2 = indices.pop(irand)
                irand = int(random.random()*len(indices))
                i3 = indices.pop(irand)
                # print('i,i1,i2,i3=',i,i1,i2,i3)
                ind1 = self.population[i1]
                ind2 = self.population[i2]
                ind3 = self.population[i3]
                v1 = ind1.vector
                v2 = ind2.vector
                v3 = ind3.vector
                vd = v1 +self.F *(v2 -v3)
                #...cross over
                vnew = np.array(vd)
                for k in range(len(vi)):
                    r = random.random()
                    if r > self.CR:
                        vnew[k] = vi[k]
                #...create new individual for trial
                self.iidmax += 1
                newind = Individual(self.iidmax, self.ndim, self.vrs, self.loss_func)
                newind.set_variable(vnew)
                candidates.append(newind)

            #...Evaluate loss func values of candidates
            #...This block can be parallelize and it makes the program much faster
            # for ic,ci in enumerate(candidates):
            #     self.kwargs['index'] = ic
            #     ci.calc_loss_func(self.kwargs)
            #...Evaluate loss function values
            qs = [ Queue() for i in range(self.N) ]
            prcs = []
            for ic,ci in enumerate(candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = ic
                prcs.append(Process(target=ci.calc_loss_func, args=(kwtmp,qs[ic])))
            for p in prcs:
                p.start()
            for p in prcs:
                p.join()
            for ic,ci in enumerate(candidates):
                ci.val = qs[ic].get()

            #...Check best
            for ic,ci in enumerate(candidates):
                if ci.val < self.bestind.val:
                    self.bestind = ci
                    self.write_variables(ci,
                                         fname='in.vars.fitpot.{0:d}'.format(ci.iid),
                                         **self.kwargs)
            if self.print_level > 2:
                for ci in candidates:
                    self.write_variables(ci,
                                         fname='in.vars.fitpot.{0:d}'.format(ci.iid),
                                         **self.kwargs)

            #...Decide whether or not to adopt new one
            for ic,ci in enumerate(candidates):
                pi = self.population[ic]
                # #...Skip if pi is the current best
                # if pi.iid == self.bestind.iid:
                #     continue
                #...adoption probability
                prob = min(1.0, np.exp(-(ci.val -pi.val)/self.T))
                r = random.random()
                if r < prob:  # replace with new individual
                    self.population[ic] = ci
                    find.write(' {0:8d}  {1:12.4e}'.format(ci.iid, ci.val))
                    for k,vk in enumerate(ci.vector):
                        find.write(' {0:11.3e}'.format(vk))
                    find.write('\n')
                else:
                    pass

            if self.print_level > 0:
                print(' step,time,best,vars= {0:8d} {1:8.1f}  {2:9.4f}'.format(it+1, time()-start,
                                                                               self.bestind.val),end="")
                for i in range(min(10,self.ndim)):
                    print(' {0:11.3e}'.format(self.bestind.vector[i]),end="")
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
    kwargs = {}
    kwargs['print_level'] = int(args['--print-level'])

    vs = np.array([1.0, -0.5])
    vrs = np.array([[-1.0, 2.0],[-1.0, 1.0]])
    
    de = DE(10, 0.8, 0.5, 1.0, vs, vrs, test_func, test_write_func, **kwargs)
    de.run(20)
    
