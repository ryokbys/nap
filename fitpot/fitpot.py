u"""Fit parameters of a certain potential to DFT data.

The potential must be specified in pmd input file, in.pmd.
"""

import os
import time
import glob
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

#.....import from local modules
from parse_input import read_input
from MD_system import MD_system
from ga import GA

#.....global variables
samples= []
sample_dirs= []
ergrefs= []
frcrefs= []
ergpmds= []
frcpmds= []
inputs={}
params=[]
prange=[]
nprms= 0
rcut= 0.0
vars= []
vranges=[]

#.....input parameters
nsmpl= 1
niter= 1
fmethod= 'test'
maindir= 'learning_set'
parfile= 'in.params.SW_Si'
runmode= 'serial'
xtol= 1e-5
gtol= 1e-5
ftol= 1e-5
eps = 1e-8
#.....GA parameters
ga_nindv= 10
ga_nbit= 16
ga_temp= 1.0

#.....constants
large= 1.0e+30
tiny = 1.0e-8

def get_sample_dirs():
    if not os.path.exists(maindir):
        print "{:*>20}: {} does not exist !!!".format(' Error',maindir)
        exit()
    lst= glob.glob(maindir+'/0*')
    for i in range(len(lst)):
        lst[i]= lst[i][len(maindir):]
    return lst

def read_pos():
    u"""read position data from learning_set/#####/pos.
    """
    global sample_dirs
    global samples
    #.....for each directory written in dir-list.txt file...
    for dir in sample_dirs:
        sys= MD_system()
        sys.read_pmd(maindir+'/'+dir+'/pos')
        sys.set_id(dir)
        samples.append(sys)

def read_params(fname):
    u"""Read potential parameters from in.params.???? file.
    The in.params.???? file should have the num of parameters and
    cutoff radius at the 1st line.
    The in.params.???? can have the range of each parameter.
    And for the parameters to be fixed can be also set.
    """
    global nprms,rcut,params,prange
    f=open(fname,'r')
    data=f.readline().split()
    nprms= int(data[0])
    rcut= float(data[1])
    params= np.zeros(nprms)
    prange= np.zeros([nprms,2])
    for i in range(nprms):
        data= f.readline().split()
        #.....depending the num of columns, the range or constraint
        #.....  of the parameter is set
        params[i]= float(data[0])
        if len(data) == 1: # no contraint for the range
            prange[i,0]= -large
            prange[i,1]=  large
        elif len(data) == 2:
            prange[i,0]=  float(data[1])
            prange[i,1]=  float(data[1])
        elif len(data) == 3:
            prange[i,0]=  float(data[1])
            prange[i,1]=  float(data[2])
    f.close()
    params_to_vars(params,prange)
    print ' number of variables to be fitted=',len(vars)

def write_params(fname,x):
    global nprms,rcut
    vars_to_params(x)
    f=open(fname,'w')
    f.write(' {:6d} {:10.4f}\n'.format(nprms,rcut))
    for i in range(nprms):
        if abs(prange[i,0]) >= large and abs(prange[i,1]) >= large:
            f.write('{:14.6e} \n'.format(params[i]))
        elif prange[i,0] == prange[i,1]:
            f.write('{:14.6e} {:14.6e}\n'.format(params[i],prange[i,0]))
        else:
            f.write('{:14.6e} {:14.6e} {:14.6e}\n'.format(params[i],
                                                          prange[i,0],
                                                          prange[i,1]))
    f.close()

def params_to_vars(x,xrange):
    global vars
    nvars= 0
    for i in range(nprms):
        if xrange[i,0] != xrange[i,1]:
            nvars += 1
    vars= np.zeros(nvars)
    vranges np.zeros(nvars,2)
    j=0
    for i in range(nprms):
        if xrange[i,0] != xrange[i,1]:
            vars[j]= x[i]
            vranges[j,0]= xrange[i,0]
            vranges[j,1]= xrange[i,1]
            j += 1

def vars_to_params(x):
    global params
    j=0
    for i in range(nprms):
        if prange[i,0] != prange[i,1]:
            params[i]= x[j]
            j += 1

def set_input_params(dict):
    global nsmpl,niter,fmethod,maindir,parfile,runmode,eps,xtol,gtol,ftol
    if 'num_samples' in dict:
        nsmpl= dict['num_samples']
    if 'num_iteration' in dict:
        niter= dict['num_iteration']
    if 'fitting_method' in dict:
        fmethod= dict['fitting_method']
    if 'main_directory' in dict:
        maindir= dict['main_directory']
    if 'param_file' in dict:
        parfile= dict['param_file']
    if 'run_mode' in dict:
        runmode= dict['run_mode']
    if 'eps' in dict:
        eps= dict['eps']
    if 'xtol' in dict:
        xtol= dict['xtol']
    if 'gtol' in dict:
        gtol= dict['gtol']
    if 'ftol' in dict:
        ftol= dict['ftol']
    #.....GA parameters
    if 'ga_num_individuals' in dict:
        ga_nindv= dict['ga_num_individuals']
    if 'ga_num_bit' in dict:
        ga_nbit= dict['ga_num_bit']
    if 'ga_temperature' in dict:
        ga_temp= dict['ga_temperature']

    
def show_input_params(input_params):
    print '>>>>> input parameters:'
    for key,value in input_params.items():
        print ' {:>20}: '.format(key), value

def gather_pmd_data():
    global ergpmds,frcpmds
    #.....initialize variables
    ergpmds=np.zeros(len(samples))
    for smpl in samples:
        frcpmds.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        #.....force
        ff=open(maindir+'/'+dir+'/frc.pmd','r')
        natm= int(ff.readline().split()[0])
        #.....energy
        f=open(maindir+'/'+dir+'/erg.pmd','r')
        ergpmds[i]= float(f.readline().split()[0])/natm
        f.close()
        for j in range(natm):
            data= ff.readline().split()
            for k in range(3):
                frcpmds[i][j,k]= float(data[k])
        ff.close()

def gather_ref_data():
    global ergrefs,frcrefs
    #.....initialize variables
    ergrefs=np.zeros(len(samples))
    for smpl in samples:
        frcrefs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        #.....force
        ff=open(maindir+'/'+dir+'/frc.ref','r')
        natm= int(ff.readline().split()[0])
        #.....energy
        f=open(maindir+'/'+dir+'/erg.ref','r')
        ergrefs[i]= float(f.readline().split()[0])/natm
        f.close()
        #.....read forces
        for j in range(natm):
            data= ff.readline().split()
            for k in range(3):
                frcrefs[i][j,k]= float(data[k])
        ff.close()

#============================================= function evaluation
def func(x,*args):
    u"""evaluate function L=sum_{samples}[E(pmd)-E(ref)]^2.
    
    This will be called from scipy.optimize.fmin_cg().
    The 1st argument x should be 1-D array of variables.
    """
    #.....write parameters to in.params.????? file
    dir= args[0]
    write_params(dir+'/'+parfile,x)
    
    #.....run pmd in all sample directories
    os.chdir(dir)
    if runmode in {'serial','Serial','SERIAL'}:
        os.system('./serial_run_pmd.sh '+parfile)
    elif runmode in {'parallel','Parallel','PARALLEL'}:
        os.system('./parallel_run_pmd.py '+parfile)
    else:
        print "{:*>20}: no such run_mode !!!".format(' Error', runmode)
        exit()
    os.chdir(cwd)

    #.....gather pmd results
    gather_pmd_data()

    #.....calc function value of L
    val= 0.0
    nval= 0
    for i in range(len(samples)):
        val += ((ergpmds[i]-ergrefs[i])/ergrefs[i])**2
        nval += 1
        natm= samples[i].natm
        for j in range(natm):
            for k in range(3):
                if abs(frcrefs[i][j,k]) > tiny:
                    val += ((frcpmds[i][j,k]-frcrefs[i][j,k])
                            /frcrefs[i][j,k])**2
                else:
                    val += ((frcpmds[i][j,k]-frcrefs[i][j,k]))**2
                nval += 1
    val /= nval
    print ' x, val=',x,val
    return val

#================================================== output data
def output_energy_relation(fname='out.pmd-vs-dft'):
    f= open(fname,'w')
    for i in range(len(ergrefs)):
        f.write(' {:15.7e} {:15.7e}\n'.format(ergrefs[i],ergpmds[i]))
    f.close()

def plot_energy_relation(fname='graph.eps'):
    x=ergrefs
    y=ergpmds
    plt.scatter(x,y)
    xmin,xmax= plt.xlim()
    ymin,ymax= plt.ylim()
    xmax= max(xmax,ymax)
    xmin= min(xmin,ymin)
    dat= np.arange(xmin,xmax,50.0)
    plt.plot(dat,dat,'b--')
    plt.xlabel('DFT energy (eV)')
    plt.ylabel('pmd energy (eV)')
    plt.title('DFT and pmd energy relation')
    plt.xlim(xmin,xmax)
    plt.ylim(xmin,xmax)
    plt.savefig(fname)
    #plt.show()

#================================================== GA wrapper
def fitfunc(val):
    return exp(-val/ga_temp)

def ga_wrapper():
    ga_check_range()
    ga= GA(ga_nindv,ga_nbit,vars,vranges,fitfunc,args=(maindir,))
    return ga.run(niter)

def ga_check_range():
    wrong=False
    for i in range(len(vars)):
        min= vranges[i,0]
        max= vranges[i,1]
        if abs(max-min) > 2.0**ga_nbit:
            wrong= True
            print ' A range seems to be too wide [{},{}]'.format(min,max)
    if wrong:
        print '{:*>20}: Some ranges are too wide.'.format(' Error')
        print '  Hoping you know what you are doing...'
        exit()
    
#============================================= main routine hereafter
if __name__ == '__main__':
    print "{:=^72}".format(' FITPOT ')
    t0= time.time()
    cwd= os.getcwd()
    #.....inputs: parameters in in.fitpot as a dictionary
    inputs= read_input('in.fitpot')
    set_input_params(inputs)
    show_input_params(inputs)
    #.....params: parameters in in.params.?????
    read_params(maindir+'/'+parfile)
    write_params(maindir+'/'+parfile
                 +'.{:03d}'.format(0),vars)
    
    #.....get samples from ##### directories
    sample_dirs= get_sample_dirs()
    if nsmpl != len(sample_dirs):
        print '{:*>20}: num_samples in in.fitpot is wrong.'.format(' Error')
        exit()
    read_pos()

    #.....initial data
    gather_ref_data()
    gather_pmd_data()

    output_energy_relation(fname='out.pmd-vs-dft.ini')

    if fmethod in {'cg','CG','conjugate-gradient'}:
        print '>>>>> conjugate-gradient was selected.'
        solution= opt.fmin_cg(func,vars,args=(maindir,),maxiter=niter,disp=True
                              ,epsilon=eps,gtol=gtol)
        print ' CG solution:',solution
    elif fmethod in {'qn','quasi-Newtown','QN','bfgs','BFGS'}:
        print '>>>>> quasi-Newton was selected.'
        solution= opt.fmin_bfgs(func,vars,args=(maindir,),maxiter=niter,disp=True
                                ,epsilon=eps,gtol=gtol)
        print ' QN solution:',solution
    elif fmethod in {'NM','Nelder-Mead','downhill-simplex'}:
        print '>>>>> Nelder-Mead was selected.'
        solution= opt.fmin(func,vars,args=(maindir,),maxiter=niter,disp=True)
        print ' NM solution:',solution
    elif fmethod in {'ga','GA','genetic-algorithm'}:
        print '>>>>> genetic algorithm was selected.'
        solution= ga_wrapper()
    elif fmethod in {'test','TEST'}:
        print '>>>>> TEST was selected.'
        func(vars,(maindir))
        solution= vars

    write_params(maindir+'/'+parfile+'.fin',solution)
    
    output_energy_relation(fname='out.pmd-vs-dft.fin')

    print '{:=^72}'.format(' FITPOT finished correctly ')
    print '   Elapsed time = {:12.2f}'.format(time.time()-t0)
