u""" Genetic Algrithm.
"""

import numpy as np
from random import random
import copy
from multiprocessing import Process, Queue
from math import exp,cos,sin,log

large= 1e+10
tiny = 1e-10
init_murate= 0.06
maxid= 0

def test_rosen(var,*args):
    import scipy.optimize as opt
    return opt.rosen(var)

def test_func(var,*args):
    x,y= var
    res= x**2 +y**2 +100*exp(-x**2 -y**2)*sin(2.0*(x+y))*cos(2*(x-y)) \
        +80*exp(-(x-1)**2 -(y-1)**2)*cos(x+4*y)*sin(2*x-y) \
        +200*sin(x+y)*exp(-(x-3)**2-(y-1)**2)
    return res

def fitfunc1(val):
    return exp(-val)

def fitfunc2(val):
    return log(1.0/val +1.0)

def dec_to_bin(dec,nbitlen):
    bin= np.zeros(nbitlen,dtype=int)
    for i in range(nbitlen):
        bin[i]= dec % 2
        dec /= 2
    return bin

def bin_to_dec(bin,nbitlen):
    dec= 0
    for i in range(nbitlen):
        dec += bin[i]*2**(i)
    return dec

def crossover(ind1,ind2):
    u"""
    Homogeneous crossover of two individuals to create a offspring
    which has some similarities to the parents.
    """
    ind= copy.deepcopy(ind1)
    nbitlen= ind.genes[0].nbitlen
    for i in range(len(ind.genes)):
        g1= ind1.genes[i]
        g2= ind1.genes[i]
        for ib in range(nbitlen):
            if g1.brep[ib] != g2.brep[ib] and random() < 0.5:
                ind.genes[i].brep[ib]= g2.brep[ib]
    return ind

def make_pairs(num):
    u"""makes random pairs from num elements.
    """
    arr= range(num)
    pairs=[]
    while len(arr) >= 2:
        i= int(random()*len(arr))
        ival= arr.pop(i)
        j= int(random()*len(arr))
        jval= arr.pop(j)
        pairs.append((ival,jval))
    return pairs

class Gene:
    u"""Gene made of *nbitlen* bits.
    """

    def __init__(self,nbitlen,var,min=-large,max=large):
        self.nbitlen= nbitlen
        self.brep = np.zeros(nbitlen,dtype=int)
        self.set_range(min,max)
        self.set_var(var)

    def set_range(self,min,max):
        self.min= float(min)
        self.max= float(max)

    def set_var(self,var):
        dec= int((var-self.min)/(self.max-self.min) *(2**self.nbitlen-1))
        self.brep= dec_to_bin(dec,self.nbitlen)

    def get_var(self):
        dec= bin_to_dec(self.brep,self.nbitlen)
        return self.min +dec*(self.max-self.min)/(2**self.nbitlen-1)

    def mutate(self,rate):
        for i in range(self.nbitlen):
            if random() < rate:
                self.brep[i] = (self.brep[i]+1) % 2



class Individual:
    u"""Individual made of some genes which should return evaluation value..
    """

    def __init__(self,id,ngene,murate,func,*args):
        self.id= id
        self.ngene= ngene
        self.murate= murate
        self.func= func
        self.args= args
        self.value= 0.0

    def set_genes(self,genes):
        if len(genes) != self.ngene:
            print "{:*>20}: len(genes) != ngene !!!".format(' Error')
            exit()
        self.genes= genes

    def calc_func_value(self,q):
        u"""
        calculates the value of given function.
        """
        vars= np.zeros(len(self.genes))
        for i in range(len(self.genes)):
            vars[i]= self.genes[i].get_var()
        val= self.func(vars,self.args)
        q.put(val)
        #self.value= self.func(vars,self.args)
        print ' ID{:05d}: value= {:15.7f}'.format(self.id,val)

    def mutate(self):
        for gene in self.genes:
            gene.mutate(self.murate)

    def set_mutation_rate(self,rate):
        self.murate= rate

    def get_variables(self):
        vars= []
        for gene in self.genes:
            vars.append(gene.get_var())
        return vars


class GA:
    u""" Genetic Algorithm class.
    """

    def __init__(self,nindv,ngene,nbitlen,func,vars,vranges,fitfunc,*args):
        u"""Constructor of GA class.

        func
          function to be evaluated with arguments vars and *args.

        fitfunc
          function for fitness evaluation with using function value obtained above.
        """
        self.nindv= nindv
        self.ngene= ngene
        self.nbitlen= nbitlen
        self.func= func
        self.vars= vars
        self.vranges= vranges
        self.fitfunc= fitfunc
        self.args= args
        self.create_population()

    def keep_best_individual(self):
        vals= []
        for i in range(len(self.population)):
            vals.append(self.population[i].value)
        idx= vals.index(min(vals))
        self.best_individual= copy.deepcopy(self.population[idx])

    def create_population(self):
        u"""creates *nindv* individuals around the initial guess."""
        global maxid
        self.population= []
        #.....0th individual is the initial guess if there is
        ind= Individual(0,self.ngene,init_murate,self.func,self.args)
        genes=[]
        for ig in range(self.ngene):
            g= Gene(self.nbitlen,self.vars[ig]
                    ,min=self.vranges[ig,0],max=self.vranges[ig,1])
            genes.append(g)
        ind.set_genes(genes)
        self.population.append(ind)
        #.....other individuals whose genes are randomly distributed
        for i in range(self.nindv-1):
            ind= Individual(i+1,self.ngene,init_murate,self.func,self.args)
            maxid= i+1
            genes= []
            for ig in range(self.ngene):
                g= Gene(self.nbitlen,self.vars[ig]
                        ,min=self.vranges[ig,0],max=self.vranges[ig,1])
                #.....randomize by mutating with high rate
                g.mutate(0.25)
                genes.append(g)
            ind.set_genes(genes)
            self.population.append(ind)
        
    def roulette_selection(self):
        u"""selects *nindv* individuals according to their fitnesses
        by means of roulette.
        """
        #.....calc all the probabilities
        prob= []
        for ind in self.population:
            prob.append(self.fitfunc(ind.value))
        print prob

        self.keep_best_individual()

        istore=[]
        for i in range(len(self.population)):
            istore.append(0)
        print istore

        for i in range(self.nindv):
            ptot= 0.0
            for ii in range(len(self.population)):
                if istore[ii] == 1: continue
                ptot += prob[ii]
            prnd= random()*ptot
            ptot= 0.0
            for ii in range(len(self.population)):
                if istore[ii] == 1: continue
                ptot= ptot +prob[ii]
                #print ii,prnd,ptot
                if prnd < ptot:
                    istore[ii]= 1
                    break
        print istore

        while istore.count(0) > 0:
            idx= istore.index(0)
            del self.population[idx]
            del istore[idx]

        if len(self.population) != self.nindv:
            print "{:*>20}: len(self.population != self.nindv) !!!".format(' Error')
            print len(self.population), self.nindv
            exit()

    def run(self,maxiter=100):
        u"""main loop of GA.
        """
        global maxid
        #.....parallel processes of function evaluations
        prcs= []
        qs= []
        for i in range(self.nindv):
            qs.append(Queue())
            prcs.append(Process(target=self.population[i].calc_func_value
                                ,args=(qs[i],)))
        for p in prcs:
            p.start()
        for p in prcs:
            p.join()
        for i in range(self.nindv):
            self.population[i].value= qs[i].get()

        for it in range(maxiter):
            print ' step= {:8d}'.format(it+1)
            #.....give birth to some offsprings by crossover
            pairs= make_pairs(self.nindv)
            for pair in pairs:
                new_ind= crossover(self.population[pair[0]],
                                   self.population[pair[1]])
                maxid += 1
                new_ind.id= maxid
                self.population.append(new_ind)
            #.....mutation of new born offsprings
            for i in range(self.nindv,len(self.population)):
                self.population[i].mutate()
            #.....evaluate function values of new born offsprings
            prcs= []
            qs= []
            j=0
            for i in range(self.nindv,len(self.population)):
                qs.append(Queue())
                prcs.append(Process(target=self.population[i].calc_func_value,
                                    args=(qs[j],)))
                j += 1
            for p in prcs:
                p.start()
            for p in prcs:
                p.join()
            j=0
            for i in range(self.nindv,len(self.population)):
                self.population[i].value= qs[j].get()
                j += 1
            #.....selection
            self.roulette_selection()
            #.....output the current best if needed
            self.out_current_best()
        print ' Best record:'
        best= self.best_individual
        print '  ID= {:05d}'.format(best.id)
        print '  value= {:15.7f}'.format(best.value)
        print '  variables= ',best.get_variables()

    def out_current_best(self):
        best= self.best_individual
        print ' current best ID{:05d}, '.format(best.id) \
            +'value= {:15.7f}'.format(best.value)
        f=open('out.current_best','w')
        f.write('ID= {:05d}\n'.format(best.id))
        f.write('value= {:15.7f}\n'.format(best.value))
        f.write('variables:\n')
        vars= best.get_variables()
        for var in vars:
            f.write('{:15.7e}\n'.format(var))
        f.close()

if __name__ == '__main__':
    vars= np.array([0.1,0.2])
    vranges= np.array([[-5.0,5.0],[-5.0,5.0]])
    ga= GA(10,2,16,test_func,vars,vranges,fitfunc1)
    ga.run(20)
