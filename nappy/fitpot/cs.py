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

import os
import sys
from docopt import docopt
import numpy as np
from numpy import exp, sin, cos
import random
import copy
from multiprocessing import Process, Pool
from time import time
from scipy.special import gamma

__author__ = "RYO KOBAYASHI"
__version__ = "190920"

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

def update_vrange(vrs,vrsh,all_indivisuals):
    """
    Update variable ranges adaptively using all the individuals information.
    """
    #...Extract top NTOPS individuals from all
    ntops = 100
    tops = []
    # print('len(all_indivisuals)=',len(all_indivisuals))
    for i,ind in enumerate(all_indivisuals):
        if len(tops) < ntops:  # add the individual
            # print(' i (< ntops)=',i)
            for it,t in enumerate(tops):
                if ind.val < t.val:
                    tops.insert(it,ind)
                    break
            if not ind in tops:
                tops.append(ind)
        else: # insert the individual and pop out the worst one
            # print(' i (>=ntops)=',i)
            for it,t in enumerate(tops):
                if ind.val < t.val:
                    tops.insert(it,ind)
                    break
            if len(tops) > ntops:
                del tops[ntops:len(tops)]

    # print('len(tops)=',len(tops))
    # print('iids= ',[t.iid for t in tops])
    #...Get new ranges
    new_vrs = np.array(vrs)
    vss = np.zeros((len(tops),len(vrs)))
    for i,ind in enumerate(tops):
        vi = ind.vector
        vss[i,:] = vi[:]
        # print('i,vi=',i,vi)
    for j in range(len(new_vrs)):
        # print('j,min,max=',i,min(vss[:,j]),max(vss[:,j]))
        new_vrs[j,0] = max(new_vrs[j,0],min(vss[:,j]))
        new_vrs[j,1] = min(new_vrs[j,1],max(vss[:,j]))

    #...Set best variables center in the ranges
    fbest = tops[0].val
    vbest = tops[0].vector
    for j in range(len(vbest)):
        vjmin = new_vrs[j,0]
        vjmax = new_vrs[j,1]
        wmax = max(abs(vjmin-vbest[j]),abs(vjmax-vbest[j]))
        new_vrs[j,0] = max(min(vjmin,vbest[j]-wmax),vrsh[j,0])
        new_vrs[j,1] = min(max(vjmax,vbest[j]+wmax),vrsh[j,1])
    
    return new_vrs

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

    def calc_loss_func(self,kwargs):
        """
        Compute loss function value using self.loss_func function given in the constructor.
        """
        # print('type(kwargs)=',type(kwargs))
        val = self.loss_func(self.vector, **kwargs)
        # print(' iid,v,val=',self.iid,self.vector,val)
        return val,kwargs['index']

class CS:
    """
    Cuckoo search class.
    """

    def __init__(self, N, F, variables, vranges, vhardlimit, loss_func, write_func,
                 nproc=0,**kwargs):
        """
        Conctructor of CS class.

        N:  Number of individuals.
        F:  Fraction of worse individuals to be abondoned.
        loss_func:
            Loss function to be minimized with variables and **kwargs.
        nproc:
            Number of processes used to run N individuals.
        """
        if N < 2:
            raise ValueError('N must be greater than 1 in CS!')
        self.N = N   # Number of individuals in a generation
        self.F = F   # Fraction of worse individuals to be abondoned
        self.nproc = nproc
        self.ndim = len(variables)
        self.vs = variables
        self.vrs0 = vranges
        self.vrs = copy.copy(self.vrs0)
        self.vrsh = vhardlimit
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
            self.print_level = int(kwargs['print_level'])
        if 'update_vrange' in kwargs.keys():
            self.update_vrs_per = kwargs['update_vrange']

        self.beta = 1.5
        self.betai = 1.0 /self.beta
        self.usgm = (gamma(1+self.beta)*np.sin(np.pi*self.beta/2)/ \
                     gamma((1+self.beta)/2)*self.beta*2.0**((self.beta-1)/2))**self.betai
        self.vsgm = 1.0

        #...initialize population
        self.population = []
        self.all_indivisuals = []
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
        prcs = []
        if self.nproc > 0 :  # use specified number of cores by nproc
            pool = Pool(processes=self.nproc)
        else:
            pool = Pool()
            
        for ip,pi in enumerate(self.population):
            kwtmp = copy.copy(self.kwargs)
            kwtmp['index'] = ip
            kwtmp['iid'] = pi.iid
            #prcs.append(pool.apply_async(pi.calc_loss_func, (kwtmp,qs[ip])))
            prcs.append(pool.apply_async(pi.calc_loss_func, (kwtmp,)))
        results = [ res.get() for res in prcs ]
        for res in results:
            val,ip = res
            self.population[ip].val = val

        pool.close()
        
        self.keep_best()
        self.all_indivisuals.extend(self.population)
        if self.print_level > 2:
            for pi in self.population:
                fname = 'in.vars.fitpot.{0:d}'.format(pi.iid)
                self.write_variables(pi,
                                     fname=fname,
                                     **self.kwargs)
        else:
            fname = 'in.vars.fitpot.{0:d}'.format(self.bestind.iid)
            self.write_variables(self.bestind,
                                 fname=fname,
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
        """
        Sort individuals in the population in the ascending order.
        """
        jtop = self.N
        for i in range(self.N):
            jtop -= 1
            for j in range(jtop):
                pj = self.population[j]
                pjp = self.population[j+1]
                if pj.val > pjp.val:
                    self.population[j] = pjp
                    self.population[j+1] = pj
        return None
        
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
        #...Headers
        fgen.write('# {0:>4s}  {1:>8s}  {2:12s}\n'.format('gen','iid','loss'))
        find.write('# {0:>7s}  {1:>12s}'.format('iid', 'loss'))
        for i in range(len(self.population[0].vector)):
            find.write(' {0:>8d}-th'.format(i+1))
        find.write('\n')
        
        for i,ind in enumerate(self.population):
            fgen.write('     0  {0:8d}  {1:12.4e}\n'.format(ind.iid, ind.val))
            find.write(' {0:8d}  {1:12.4e}'.format(ind.iid, ind.val))
            for j,vj in enumerate(ind.vector):
                find.write(' {0:11.4e}'.format(vj))
            find.write('\n')

        if self.print_level > 0:
            print(' step,time,best,vars= {0:6d} {1:8.1f}  {2:8.4f}'.format(0, time()-start,
                                                                           self.bestind.val),end="")
            for i in range(min(16,self.ndim)):
                print(' {0:6.3f}'.format(self.bestind.vector[i]),end="")
            print('', flush=True)

        #...Create pool before going into maxiter-loop,
        #...since creating pool inside could cause "Too many files" error.
        if self.nproc > 0 :  # use specified number of cores by nproc
            pool = Pool(processes=self.nproc)
        else:
            pool = Pool()
            
        for it in range(maxiter):
            self.sort_individuals()
            #...Create candidates from current population using Levy flight
            candidates = []
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

            #...Create new completely random candidates
            iab = int((1.0 -self.F)*self.N)
            rnd_candidates = []
            for iv in range(iab,self.N):
                self.iidmax += 1
                newind = Individual(self.iidmax, self.ndim, self.vrs, self.loss_func)
                newind.init_random()
                rnd_candidates.append(newind)

            #...Evaluate loss function values of updated candidates and new random ones
            prcs = []
            for ic,ci in enumerate(candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = ic
                kwtmp['iid'] = ci.iid
                # prcs.append(Process(target=ci.calc_loss_func, args=(kwtmp,qs[ic])))
                prcs.append(pool.apply_async(ci.calc_loss_func, (kwtmp,)))
            rnd_prcs = []
            for ic,ci in enumerate(rnd_candidates):
                kwtmp = copy.copy(self.kwargs)
                kwtmp['index'] = len(candidates) +ic
                kwtmp['iid'] = ci.iid
                # prcs.append(Process(target=ci.calc_loss_func, args=(kwtmp,qs[ic])))
                rnd_prcs.append(pool.apply_async(ci.calc_loss_func, (kwtmp,)))
            
            results = [ res.get() for res in prcs ]
            rnd_results = [ res.get() for res in rnd_prcs ]

            for res in results:
                val,ic = res
                candidates[ic].val = val
            self.all_indivisuals.extend(candidates)

            for res in rnd_results:
                val,ic_rnd = res
                ic = ic_rnd -len(candidates)
                rnd_candidates[ic].val = val
            self.all_indivisuals.extend(rnd_candidates)

            #...Pick j that is to be compared with i
            js = random.sample(range(self.N),k=self.N)
            #...Decide whether or not to adopt new one
            for jc,jv in enumerate(js):
                pj = self.population[jv]
                cj = candidates[jc]
                dval = cj.val -pj.val
                if dval < 0.0:  # replace with new individual
                    self.population[jv] = cj
                    find.write(' {0:8d}  {1:12.4e}'.format(cj.iid, cj.val))
                    for k,vk in enumerate(cj.vector):
                        find.write(' {0:11.4e}'.format(vk))
                    find.write('\n')
                    find.flush()
                else:
                    pass

            #...Rank individuals
            self.sort_individuals()
            
            #...Replace to-be-abandoned ones with new random ones
            ic = 0
            for iv in range(iab,self.N):
                ci = rnd_candidates[ic]
                ic += 1
                self.population[iv] = ci

            #...Check best
            best_updated = False
            for ic,ci in enumerate(self.population):
                if ci.val < self.bestind.val:
                    self.bestind = ci
                    best_updated = True
            if best_updated:
                fname = 'in.vars.fitpot.{0:d}'.format(self.bestind.iid)
                self.write_variables(self.bestind,
                                     fname=fname,
                                     **self.kwargs)
                os.system('cp -f {0:s} in.vars.fitpot.best'.format(fname))

            #...Update variable ranges if needed
            if self.update_vrs_per > 0 and (it+1) % self.update_vrs_per == 0:
                self.vrs = update_vrange(self.vrs,self.vrsh,self.all_indivisuals)
                print(' Update variable ranges')
                for i in range(len(self.vrs)):
                    print(' {0:2d}:  {1:7.3f}  {2:7.3f}'.format(i+1,self.vrs[i,0],self.vrs[i,1]))
                #...Set variable ranges of all individuals in the population
                for iv in range(len(self.population)):
                    self.population[iv].vranges = self.vrs
            
            if self.print_level > 0:
                print(' step,time,best,vars= {0:6d} {1:8.1f}  {2:8.4f}'.format(it+1, time()-start,
                                                                               self.bestind.val),end="")
                for i in range(min(16,self.ndim)):
                    print(' {0:6.3f}'.format(self.bestind.vector[i]),end="")
                print('', flush=True)

            for i,ind in enumerate(self.population):
                fgen.write(' {0:5d}  {1:8d}  {2:12.4e}\n'.format(it+1, ind.iid, ind.val))
                fgen.flush()
        fgen.close()
        find.close()
        pool.close()
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
    
