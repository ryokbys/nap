#!/usr/local/bin/python
# -*- coding: utf-8 -*-

"""Fit parameters of a certain potential to DFT data.

The potential must be specified in pmd input file, in.pmd.
"""

import os,sys
import time
import glob
import numpy as np
import scipy.optimize as opt
# import matplotlib.pyplot as plt
import math
import multiprocessing as mp

#.....import from local modules
from MD_system import MD_system
import NN
from pga import GA

#.....constants
_large= 1.0e+30
_tiny = 1.0e-8
max_species= 10
_valmax= 1.0e+30


#.....global variables
samples= []
sample_dirs= []
ergrefs= []
frcrefs= []
ergpmds= []
frcpmds= []
inputs={}
params=[]
pranges=[]
nprms= 0
rcut= 0.0
vars= []
vranges=[]
l1st_call= True
bases= []  # bases in linreg that can be recycled during fitting
bmax= []
_valmin= _valmax
_init_time= 0.0

#.....input parameters
nsmpl= 1
niter= 1
fmethod= 'test'
maindir= 'learning_set'
parfile= 'in.params.SW_Si'
runmode= 'serial'
potential= 'none'
gradient= 'none'
grad_scale= False
xtol= 1e-5
gtol= 1e-5
ftol= 1e-5
eps = 1e-8
eatom= np.zeros(max_species)
fmatch= True
regularize= False
penalty= 'no'
pweight= 1.0
lswgt= False
swgt= []
swbeta= 1.0
nprcs= 1
#.....GA parameters
ga_nindv= 10
ga_nbit= 16
ga_temp= 1.0
ga_murate= 0.01

def get_sample_dirs():
    if not os.path.exists(maindir):
        print "{0:*>20}: {1} does not exist !!!".format(' Error',maindir)
        exit()
    lst= glob.glob(maindir+'/[0-9]*')
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
    global nprms,rcut,params,pranges,vars,vranges
    f=open(fname,'r')
    data=f.readline().split()
    nprms= int(data[0])
    rcut= float(data[1])
    params= np.zeros(nprms)
    pranges= np.zeros([nprms,2])
    for i in range(nprms):
        data= f.readline().split()
        #.....depending the num of columns, the range or constraint
        #.....  of the parameter is set
        params[i]= float(data[0])
        if len(data) == 1: # no contraint for the range
            pranges[i,0]= -_large
            pranges[i,1]=  _large
        elif len(data) == 2:
            pranges[i,0]=  float(data[1])
            pranges[i,1]=  float(data[1])
        elif len(data) == 3:
            pranges[i,0]=  float(data[1])
            pranges[i,1]=  float(data[2])
    f.close()
    vars,vranges= params_to_vars(params,pranges)
    print ' number of variables to be fitted=',len(vars)

def write_params(fname,x):
    global params,pranges
    
    params,pranges= vars_to_params(x,vranges,params,pranges)
    f=open(fname,'w')
    f.write(' {0:6d} {1:10.4f}\n'.format(len(params),rcut))
    #if potential in ('linreg') and not fmethod in ('test','TEST'):
    if potential in ('linreg') and regularize:
        print '>>>>> writing params by multiplying bmax...'
        print ' potential =',potential
        print ' regularize=',regularize
        for i in range(nprms):
            if abs(pranges[i,0]) >= _large and abs(pranges[i,1]) >= _large:
                f.write('{0:22.14e} \n'.format(params[i]/bmax[i]))
            elif pranges[i,0] == pranges[i,1]:
                f.write('{0:22.14e} {1:22.14e}\n'.format(params[i]/bmax[i],pranges[i,0]))
            else:
                f.write('{0:22.14e} {1:22.14e} {2:22.14e}\n'.format(params[i]/bmax[i],
                                                                    pranges[i,0],
                                                                    pranges[i,1]))
    else:
        for i in range(nprms):
            if abs(pranges[i,0]) >= _large and abs(pranges[i,1]) >= _large:
                f.write('{0:22.14e} \n'.format(params[i]))
            elif pranges[i,0] == pranges[i,1]:
                f.write('{0:22.14e} {1:22.14e}\n'.format(params[i],pranges[i,0]))
            else:
                f.write('{0:22.14e} {1:22.14e} {2:22.14e}\n'.format(params[i],
                                                                    pranges[i,0],
                                                                    pranges[i,1]))
    f.close()

def params_to_vars(params,pranges):
    """
    Converts params to vars.
    Params and pranges have richer information than vars and vranges.
    Thus the information will be reduced to vars.
    """
    nvars= 0
    for i in range(len(params)):
        if pranges[i,0] != pranges[i,1]:
            nvars += 1
    vars= np.zeros(nvars)
    vranges= np.zeros([nvars,2])
    j=0
    for i in range(len(params)):
        if pranges[i,0] != pranges[i,1]:
            vars[j]= params[i]
            vranges[j,0]= pranges[i,0]
            vranges[j,1]= pranges[i,1]
            j += 1
    return (vars,vranges)

def vars_to_params(vars,vranges,params,pranges):
    """
    Converts vars to params.
    Length of vars should be less or equal to length of params.
    Thus vars and vranges do not have enough information of params and pranges.
    """
    j=0
    for i in range(len(params)):
        if pranges[i,0] != pranges[i,1]:
            params[i]= vars[j]
            j += 1
    return (params,pranges)

def read_input(fname='in.fitpot'):
    global _conf
    global nsmpl,niter,fmethod,maindir,parfile,runmode,eps \
           ,xtol,gtol,ftol,fmatch,penalty,pweight,gradient,potential \
           ,regularize
    global lswgt,swbeta,nprcs
    global eatom,ga_nindv,ga_nbit,ga_temp

    f= open(fname,'r')
    for line in f.readlines():
        data= line.split()
        # skip if the line is empty or comment line
        if len(data)==0 or \
                line[0]=='#' or \
                line[0]=='!' or \
                line[0]=='%':
            continue
        else:
            if data[0] == 'num_samples':
                nsmpl= int(data[1])
            elif data[0] == 'num_iteration':
                niter= int(data[1])
            elif data[0] == 'fitting_method':
                fmethod= data[1]
            elif data[0] == 'main_directory':
                maindir= data[1]
            elif data[0] == 'param_file':
                parfile= data[1]
            elif data[0] == 'run_mode':
                runmode= data[1]
            elif data[0] == 'eps':
                eps= float(data[1])
            elif data[0] == 'xtol':
                xto= float(data[1])
            elif data[0] == 'gtol':
                gtol= float(data[1])
            elif data[0] == 'ftol':
                ftol= float(data[1])
            elif data[0] == 'atom_energy':
                eatom[int(data[1])]= float(data[2])
            elif data[0] == 'force_match':
                if data[1] in ('false','False','no','No','0'):
                    fmatch= False
            elif data[0] == 'penalty':
                penalty= data[1]
            elif data[0] == 'penalty_weight':
                pweight= float(data[1])
            elif data[0] == 'potential':
                potential= data[1]
            elif data[0] == 'gradient':
                gradient= data[1]
            elif data[0] == 'grad_scale':
                if data[1] in ('true','True','yes','YES','y','Y','1'):
                    grad_scale= True
            elif data[0] == 'regularize':
                if data[1] in ('true','True','yes','YES','y','Y','1'):
                    regularize= True
            elif data[0] == 'sample_weight':
                if data[1] in ('true','True','yes','YES','y','Y','1'):
                    lswgt= True
            elif data[0] == 'sample_weight_beta':
                swbeta= float(data[1])
            elif data[0] == 'num_multiprocess':
                nprcs= int(data[1])
            #.....GA parameters
            elif data[0] == 'ga_num_individuals':
                ga_nindv= int(data[1])
            elif data[0] == 'ga_num_bit':
                ga_nbit= int(data[1])
            elif data[0] == 'ga_temperature':
                ga_temp= float(data[1])
            elif data[0] == 'ga_murate':
                ga_murate= float(data[1])
            else:
                print '{0:*^72}'.format('*')
                print '* There is no keyword like: ',data[0]
                print '* Please confirm the keyword again...'
                print '{0:*^72}'.format('*')
    f.close()

def show_inputs(fname='in.fitpot'):
    print '>>>>> configuration:'
    # for key,value in input_params.items():
    #     print ' {0:>20}: '.format(key), value
    f= open(fname,'r')
    for line in f.readlines():
        data= line.split()
        # skip if the line is empty or comment line
        if len(data)==0 or \
                line[0]=='#' or \
                line[0]=='!' or \
                line[0]=='%':
            continue
        else:
            print '  '+line.rstrip()
    f.close()

def gather_pmd_data(basedir):
    global samples,sample_dirs
    # print ' basedir=',basedir
    # print ' len(samples)=',len(samples)
    # print ' len(sample_dirs)=',len(sample_dirs)
    #.....initialize variables
    ergs=np.zeros(len(samples))
    frcs= []
    for smpl in samples:
        frcs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        #.....force
        ff=open(basedir+'/'+dir+'/frc.pmd','r')
        natm= int(ff.readline().split()[0])
        #.....energy
        f=open(basedir+'/'+dir+'/erg.pmd','r')
        ergs[i]= float(f.readline().split()[0])
        f.close()
        for j in range(natm):
            data= ff.readline().split()
            for k in range(3):
                frcs[i][j,k]= float(data[k])
        ff.close()
    # print 'ergs:',ergs
    # print 'frcs:',frcs
    return (ergs,frcs)

def gather_smd_data(basedir):
    global samples,sample_dirs
    # print ' basedir=',basedir
    # print ' len(samples)=',len(samples)
    # print ' len(sample_dirs)=',len(sample_dirs)
    #.....initialize variables
    ergs=np.zeros(len(samples))
    frcs= []
    for smpl in samples:
        frcs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        #.....force
        ff=open(basedir+'/'+dir+'/frc.smd','r')
        natm= int(ff.readline().split()[0])
        #.....energy
        f=open(basedir+'/'+dir+'/erg.smd','r')
        ergs[i]= float(f.readline().split()[0])
        f.close()
        for j in range(natm):
            data= ff.readline().split()
            for k in range(3):
                frcs[i][j,k]= float(data[k])
        ff.close()
    # print 'ergs:',ergs
    # print 'frcs:',frcs
    return (ergs,frcs)

def gather_ref_data(basedir):
    global ergrefs,frcrefs
    global swgt
    
    #.....initialize variables
    ergrefs=np.zeros(len(samples))
    for smpl in samples:
        frcrefs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        #print dir
        #.....force
        ff=open(basedir+'/'+dir+'/frc.ref','r')
        natm= int(ff.readline().split()[0])
        #print 'ismpl,natm=',i,natm
        #.....energy
        f=open(basedir+'/'+dir+'/erg.ref','r')
        ergrefs[i]= float(f.readline().split()[0])
        #.....need to subtract atomic energies from total energy
        #.....to get the cohesive energy
        for isp in range(1,max_species):
            num= smpl.num_of_species(isp)
            if num != 0:
                ergrefs[i] -= eatom[isp]*num
        f.close()
        #.....read forces
        for j in range(natm):
            data= ff.readline().split()
            for k in range(3):
                frcrefs[i][j,k]= float(data[k])
        ff.close()

    #.....calc sample weight if required
    swgt= np.array([ 1.0 for i in range(len(samples))])
    if lswgt:
        #.....get minimum energy
        emin= 1.0e+30
        for ismpl in range(len(samples)):
            natm= samples[ismpl].natm
            emin= min(emin,ergrefs[ismpl]/natm)
        #.....compute sample weight
        for ismpl in range(len(samples)): 
            swgt[ismpl]= math.exp(-swbeta*(ergrefs[ismpl]/natm \
                                           -emin))

#============================================= function evaluation
def func(x,*args):
    """evaluate function L=sum_{samples}[E(pmd)-E(ref)]^2.
    
    This will be called from scipy.optimize.fmin_cg().
    The 1st argument x should be 1-D array of variables.
    """
    global _valmin
    
    t0= time.time()
    #.....write parameters to in.params.????? file
    dir= args[0]

#    if fmethod in ('test','TEST') or \
    if fmethod in ('test','TEST','check_grad') or \
            not potential in ('linreg','NN'):
        #.....store original file
        os.system('cp '+dir+'/'+parfile+' '+dir+'/'+parfile+'.tmp')
        write_params(dir+'/'+parfile,x)
        #.....run smd in all sample directories
        os.chdir(dir)
        #print os.getcwd(),dir
        if runmode in ('serial','Serial','SERIAL','sequential','single'):
            os.system('./serial_run_smd.sh '+parfile)
        elif runmode in ('parallel','Parallel','PARALLEL'):
            os.system('python ./parallel_run_smd.py '+parfile)
        else:
            print "{0:*>20}: no such run_mode !!!".format(' Error', runmode)
            exit()
        os.chdir(cwd)
        #.....restore original file
        os.system('cp '+dir+'/'+parfile+' '+dir+'/'+parfile+'.current')
        os.system('cp '+dir+'/'+parfile+'.tmp'+' '+dir+'/'+parfile)
        #.....gather smd results
        ergs,frcs=gather_smd_data(dir)
    elif potential in ('linreg'):
        #.....calc ergs and frcs from bases data and x (variables)
        read_bases(dir)
        ergs,frcs=calc_ef_from_bases(x,*args)
    elif potential in ('NN'):
        #.....now it is possible to compute only from bases
        ergs,frcs= NN.calc_ef_from_bases(x,*args)

    #.....calc function value of L
    val= eval_L(ergs,frcs,ergrefs,frcrefs,samples)
    #.....output temporal results
    output_energy_relation(ergs,ergrefs,samples,sample_dirs, \
                               fname='out.erg.pmd-vs-dft.tmp')
    output_force_relation(frcs,frcrefs,samples,sample_dirs, \
                              fname='out.frc.pmd-vs-dft.tmp')

    print
    print ' L value=',val

    if penalty in ('ridge','Ridge','RIDGE') and potential in ('linreg'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            p += math.sqrt(x[n]**2)
        print ' penalty value=',p*pweight
        val += p*pweight
        print ' total L value=',val

    elif penalty in ('lasso','LASSO') and potential in ('linreg'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            p += abs(x[n])
        print ' penalty value=',p*pweight
        val += p*pweight
        print ' total L value=',val
    sys.stdout.flush()

    #.....if L value is minimum ever, store this parameter file
    if val < _valmin:
        _valmin= val
        if potential in ('linreg','NN'):
            write_params(dir+'/'+parfile+'.min',x)
        else:
            os.system('cp '+dir+'/'+parfile+'.current' \
                          +' '+dir+'/'+parfile+'.min')
        
    print ' ===> time func: {0:12.3f} sec'.format(time.time()-t0) \
          +', {0:12.3f} sec'.format(time.time()-_init_time)
    return val

def eval_L(cergs,cfrcs,rergs,rfrcs,samples):
    val= 0.0
    for i in range(len(samples)):
        natm= samples[i].natm
        sw= 1.0
        if lswgt:
            sw= swgt[i]
        vi= (cergs[i]-rergs[i])**2 /natm
        vi *= sw
        val += vi
        if not fmatch:
            continue
        for j in range(natm):
            for k in range(3):
                val += (cfrcs[i][j,k]-rfrcs[i][j,k])**2  \
                       /(3*natm) *sw
    return val

def calc_ef_from_bases(x,*args):
    """Calculate energies and forces of every samples using bases data.
    """

    #.....initialize variables
    es=np.zeros(len(samples))
    fs= []
    for smpl in samples:
        fs.append(np.zeros((smpl.natm,3)))

    p= mp.Pool(nprcs)
    
    if potential in ('linreg'):
        # #.....calc energies
        # for ismpl in range(len(samples)):
        #     smpl= samples[ismpl]
        #     for ia in range(smpl.natm):
        #         for iprm in range(len(params)):
        #             es[ismpl] += x[iprm] *bases[ismpl][ia,iprm]
        # 
        # #.....calc forces
        # if fmatch:
        #     for ismpl in range(len(samples)):
        #         smpl= samples[ismpl]
        #         for ia in range(smpl.natm):
        #             for iprm in range(len(params)):
        #                 dbs=bases[len(samples)+ismpl][ia,iprm]
        #                 fs[ismpl][ia,0] -= x[iprm] *dbs[0]
        #                 fs[ismpl][ia,1] -= x[iprm] *dbs[1]
        #                 fs[ismpl][ia,2] -= x[iprm] *dbs[2]
        if nprcs == 1:
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                est,fst= calc_ef_linreg(ismpl,x,*args)
                es[ismpl]= est
                for ia in range(smpl.natm):
                    fs[ismpl][ia,0] += fst[ia,0]
                    fs[ismpl][ia,1] += fst[ia,1]
                    fs[ismpl][ia,2] += fst[ia,2]
        else:
            func_args=[]
            for ismpl in range(len(samples)):
                func_args.append( (calc_ef_linreg,ismpl,x) )
            results= p.map(arg_wrapper,func_args)
            p.close()
            p.join()
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                est,fst= results[ismpl]
                es[ismpl]= est
                for ia in range(smpl.natm):
                    fs[ismpl][ia,0] += fst[ia,0]
                    fs[ismpl][ia,1] += fst[ia,1]
                    fs[ismpl][ia,2] += fst[ia,2]

    return (es,fs)

def arg_wrapper(args):
    return args[0](*args[1:])

def read_bases(dir):
    global bases,l1st_call,_gsf,_hl1,_aml

    #.....read bases from pmd directories
    if l1st_call:
        if potential in ('linreg'):
            bases= gather_basis_linreg(dir)
            if regularize:
                regularize_bases_linreg(bases)
        l1st_call= False

def scale_vars(x,fac):
    if len(x) != len(fac):
        print ' [Error] len(x) != len(fac) !!!'
        exit()
    newx= np.zeros(len(x))
    for i in range(len(x)):
        newx[i]= x[i]*fac[i]
    return newx

#==================================================== linreg
def calc_ef_linreg(ismpl,x,*args):
    smpl= samples[ismpl]
    es= 0.0
    fs= np.zeros((smpl.natm,3))
    #.....calc energy
    for ia in range(smpl.natm):
        for iprm in range(len(params)):
            es += x[iprm] *bases[ismpl][ia,iprm]
    #.....calc forces
    if fmatch:
        basis= bases[len(samples)+ismpl]
        for ia in range(smpl.natm):
            for iprm in range(len(params)):
                dbs=basis[ia,iprm]
                fs[ia,0] -= x[iprm] *dbs[0]
                fs[ia,1] -= x[iprm] *dbs[1]
                fs[ia,2] -= x[iprm] *dbs[2]
    return (es,fs)

def grad_linreg(x,*args):
    global bases
    t0= time.time()
    dir= args[0]

    read_bases(dir)
    #bases= gather_basis_linreg(dir)
    #.....gather pmd results
    ergs,frcs= calc_ef_from_bases(x,dir)
    # ergs,frcs= gather_pmd_data(dir)

    p=mp.Pool(nprcs)
    
    grad= np.zeros(len(params))
    # for ismpl in range(len(samples)):
    #     smpl= samples[ismpl]
    #     sw= 1.0
    #     if lswgt:
    #         sw= swgt[ismpl]
    #     ediff= (ergs[ismpl]-ergrefs[ismpl])/smpl.natm *sw
    #     for iprm in range(len(params)):
    #         bs= 0.0
    #         for ia in range(smpl.natm):
    #             bs += bases[ismpl][ia,iprm]
    #         grad[iprm] += 2.0*ediff*bs/smpl.natm
    #         
    # if fmatch:
    #     for ismpl in range(len(samples)):
    #         smpl= samples[ismpl]
    #         sw= 1.0
    #         if lswgt:
    #             sw= swgt[ismpl]
    #         for iprm in range(len(params)):
    #             dlprm= 0.0
    #             dbs= np.zeros(3)
    #             fdiff= np.zeros(3)
    #             for ia in range(smpl.natm):
    #                 fdiff= frcs[ismpl][ia] -frcrefs[ismpl][ia]
    #                 dbs= bases[len(samples)+ismpl][ia,iprm]
    #                 dlprm -= 2.0*( fdiff[0]*dbs[0] \
    #                                    +fdiff[1]*dbs[1] \
    #                                    +fdiff[2]*dbs[2] ) \
    #                                    /(smpl.natm*3) *sw
    #             grad[iprm] += dlprm
    if nprcs == 1:
        for ismpl in range(len(samples)):
            gs= grad_linreg_core(ismpl,ergs,frcs)
            for iprm in range(len(params)):
                grad[iprm] += gs[iprm]
    else:
        func_args=[]
        for ismpl in range(len(samples)):
            func_args.append( (grad_linreg_core,ismpl,ergs,frcs) )
        results= p.map(arg_wrapper,func_args)
        p.close()
        p.join()
        for ismpl in range(len(samples)):
            gs= results[ismpl]
            for iprm in range(len(params)):
                grad[iprm] += gs[iprm]

    if penalty in ('ridge','Ridge','RIDGE'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            grad[n] += 2.0*x[n] *pweight

    elif penalty in ('lasso','LASSO'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            grad[n] += pweight *np.sign(x[n])
    print ' ===> time grad_linreg: {0:12.3f} sec'.format(time.time()-t0) \
          +', {0:12.3f} sec'.format(time.time()-_init_time)
    if grad_scale:
        maxgrad= np.max(np.abs(grad))
        maxx= np.max(np.abs(x))
        print ' maxgrad,maxx=',maxgrad,maxx
        for i in range(len(grad)):
            grad[i]= grad[i] /maxgrad *maxx/10
    #print ' grad after: ',grad
    return grad

def grad_linreg_core(ismpl,ergs,frcs):
    gs= np.zeros(len(params))
    smpl= samples[ismpl]
    sw= 1.0
    if lswgt:
        sw= swgt[ismpl]
    ediff= (ergs[ismpl]-ergrefs[ismpl])/smpl.natm *sw
    basis= bases[ismpl]
    for iprm in range(len(params)):
        bs= 0.0
        for ia in range(smpl.natm):
            bs += basis[ia,iprm]
        gs[iprm] += 2.0*ediff*bs
            
    if fmatch:
        basis= bases[len(samples)+ismpl]
        for iprm in range(len(params)):
            dlprm= 0.0
            dbs= np.zeros(3)
            fdiff= np.zeros(3)
            for ia in range(smpl.natm):
                fdiff= frcs[ismpl][ia] -frcrefs[ismpl][ia]
                dbs= basis[ia,iprm]
                dlprm -= 2.0*( fdiff[0]*dbs[0] \
                               +fdiff[1]*dbs[1] \
                               +fdiff[2]*dbs[2] ) \
                               /(smpl.natm*3) *sw
            gs[iprm] += dlprm
    return gs

def gather_basis_linreg(basedir):
    bdata= []
    #.....read basis data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        f=open(basedir+'/'+dir+'/pmd/out.basis.linreg','r')
        data= f.readline().split()
        natm= int(data[0])
        nelem=  int(data[1])
        basis= np.zeros((smpl.natm,len(params)))
        for ia in range(smpl.natm):
            for ip in range(len(params)):
                data= f.readline().split()
                basis[ia,ip]= float(data[3])
        bdata.append(basis)
        f.close()

    if not fmatch:
        return bdata

    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        g=open(basedir+'/'+dir+'/pmd/out.dbasis.linreg','r')
        data= g.readline().split()
        natm2= int(data[0])
        nelem2=int(data[1])
        dbasis= np.zeros((smpl.natm,len(params),3))
        for ia in range(smpl.natm):
            for ip in range(len(params)):
                data= g.readline().split()
                dbasis[ia,ip,0]= float(data[2])
                dbasis[ia,ip,1]= float(data[3])
                dbasis[ia,ip,2]= float(data[4])
        bdata.append(dbasis)
        g.close()
    return bdata

def regularize_bases_linreg(bases):
    """Regularize bases linearly.
    """
    global bmax

    #.....compute max of each basis
    #print ' Max value of each basis:'
    bmax= np.zeros(len(params))
    for ip in range(len(params)):
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            basis= bases[ismpl]
            for ia in range(smpl.natm):
                bmax[ip]= max(bmax[ip],basis[ia,ip])
        #print '   ip,bmax[ip]=',ip,bmax[ip]

    #.....regularize bases
    for ip in range(len(params)):
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            for ia in range(smpl.natm):
                bases[ismpl][ia,ip] /= bmax[ip]
    if fmatch:
        for ip in range(len(params)):
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                basis= bases[len(samples)+ismpl]
                for ia in range(smpl.natm):
                    basis[ia,ip,0] /= bmax[ip]
                    basis[ia,ip,1] /= bmax[ip]
                    basis[ia,ip,2] /= bmax[ip]
    
    

def output_energy_relation(es,erefs,samples,sdirs, \
                               fname='out.erg.pmd-vs-dft'):
    f= open(fname,'w')
    for i in range(len(erefs)):
        smpl= samples[i]
        sdir= sdirs[i]
        f.write(' {0:15.7e} {1:15.7e} {2:s}\n'.format(erefs[i]/smpl.natm \
                                                    ,es[i]/smpl.natm \
                                                    ,sdir))
    f.close()
    
def output_force_relation(fs,frefs,samples,sdirs, \
                              fname='out.frc.pmd-vs-dft'):
    f= open(fname,'w')
    for i in range(len(samples)):
        sdir= sdirs[i]
        for j in range(samples[i].natm):
            for k in range(3):
                f.write(' {0:15.7e} {1:15.7e} {2:s}\n'.format(frefs[i][j,k] \
                                                                  ,fs[i][j,k] \
                                                                  ,sdir))
    f.close()

def output_statistics(ergs,frcs):
    print '>>>>> statistics:'
    #.....statistics of energies
    demax= 0.0
    desum= 0.0
    for i in range(len(samples)):
        smpl= samples[i]
        de= abs(ergs[i]-ergrefs[i])/smpl.natm
        #print ' ismpl,natm,erg,ergref,de=',i,smpl.natm,ergs[i],ergrefs[i],de
        demax= max(demax,de)
        desum += de**2/len(samples)
    rmse= math.sqrt(desum)
    print '  RMSE of energies        = {0:12.3f} eV/atom'.format(rmse)
    print '  Max residual of energies= {0:12.3f} eV/atom'.format(demax)

    #.....statistics of forces
    dfmax= 0.0
    dfsum= 0.0
    n= 0
    for i in range(len(samples)):
        smpl= samples[i]
        for j in range(smpl.natm):
            for k in range(3):
                df= abs(frcs[i][j,k]-frcrefs[i][j,k])
                dfmax= max(dfmax,df)
                dfsum += df**2
                n += 1
    rmse= math.sqrt(dfsum/n)
    print '  RMSE of forces          = {0:12.3f} eV/A'.format(rmse)
    print '  Max residual of forces  = {0:12.3f} eV/A'.format(dfmax)
    


#================================================== GA wrapper
def fitfunc1(val):
    return math.exp(-val/ga_temp)

def fitfunc2(val):
    return math.log(1.0 +val)

def ga_wrapper():
    ga_check_range()
    ga= GA(ga_nindv,ga_nbit,ga_murate,func,vars,vranges \
               ,fitfunc1,args=(maindir,))
    return ga.run(niter)

def ga_check_range():
    wrong=False
    for i in range(len(vars)):
        min= vranges[i,0]
        max= vranges[i,1]
        if abs(max-min) > 2.0**ga_nbit:
            wrong= True
            print ' A range seems to be too wide [{0},{1}]'.format(min,max)
    if wrong:
        print '{0:*>20}: Some ranges are too wide.'.format(' Error')
        print '  Hoping you know what you are doing...'
        exit()

#============================================= steepest descent dynamics
def sd_dynamics(f,x,args=(),fprime=None,maxiter=10):
    u"""
    Steepest descent dynamics is peformed with using function, f,
    variables, x, and arguments, args.
    """

    #...maximum displacement of weight
    maxdisp= 2.0e-3

    if fprime is None:
        print ' [Error] fprime should be specified in sd_dynamics !'
        exit()

    print '>>>>> sd_dynamics'
    print ' maxiter=',maxiter
    print ' args   =',args
    
    val= f(x,args[0])
    print ' initial value= {0:20.7f}'.format(val)
    grad= fprime(x,args[0])
    maxgrad= np.max(grad)
    alpha= maxdisp /maxgrad
    print " maxgrad,alpha=",maxgrad,alpha

    for it in range(maxiter):
        grad= fprime(x,args[0])
        maxgrad= np.max(grad)
        alpha= maxdisp /maxgrad
        print " maxgrad,alpha=",maxgrad,alpha
        x += -alpha *grad
        val=f(x,args[0])
    print ' final value= {0:20.7f}'.format(val)
    return x
    
#============================================= main routine hereafter
if __name__ == '__main__':
    print "{0:=^72}".format(' FITPOT ')
    _init_time= time.time()
    cwd= os.getcwd()
    #.....inputs: parameters in in.fitpot as a dictionary
    inputs= read_input('in.fitpot')
    show_inputs('in.fitpot')
    #.....params: parameters in in.params.?????
    read_params(maindir+'/'+parfile)
    os.system('cp '+maindir+'/'+parfile+' '+maindir+'/'+parfile+'.ini')
    #write_params(maindir+'/'+parfile+'.ini',vars)
    
    #.....get samples from ##### directories
    sample_dirs= get_sample_dirs()
    sample_dirs.sort()
    if nsmpl != len(sample_dirs):
        print '{0:*>20}: num_samples in in.fitpot is wrong.'.format(' Error')
        exit()
    read_pos()

    #.....initial data
    gather_ref_data(maindir)
    #.....read bases data if needed
    if potential in ('linreg') and not fmethod in ('test','TEST'):
        read_bases(maindir)
        if regularize:
            vars= scale_vars(vars,bmax)
    elif potential in ('NN') and not fmethod in ('test','TEST'):
        NN.init(maindir,params,sample_dirs,samples,nprcs,fmatch \
                ,ergrefs,frcrefs,fmethod,parfile,runmode,rcut,pranges \
                ,vranges)

    #.....1st call of func
    func(vars,maindir)
    if potential in ('linreg') and not fmethod in ('test','TEST'):
        ergs,frcs= calc_ef_from_bases(vars,maindir)
    elif potential in ('NN') and not fmethod in ('test','TEST'):
        ergs,frcs= NN.calc_ef_from_bases(vars)
    else:
        ergs,frcs= gather_pmd_data(maindir)

    if fmethod in ('test','TEST') and potential in ('NN'):
        NN.init(maindir,params,sample_dirs,samples,nprcs,fmatch \
                ,ergrefs,frcrefs,fmethod,parfile,runmode \
                ,rcut,pranges,vranges)

    output_energy_relation(ergs,ergrefs,samples,sample_dirs,fname='out.erg.pmd-vs-dft.ini')
    output_force_relation(frcs,frcrefs,samples,sample_dirs,fname='out.frc.pmd-vs-dft.ini')

    if fmethod in ('cg','CG','conjugate-gradient'):
        print '>>>>> conjugate-gradient was selected.'
        if gradient in ('numerical'):
            solution= opt.fmin_cg(func,vars,args=(maindir,)
                                  ,maxiter=niter,disp=True
                                  ,epsilon=eps,gtol=gtol)
        else:
            if potential in ('linreg'):
                solution= opt.fmin_cg(func,vars,args=(maindir,)
                                      ,fprime=grad_linreg
                                      ,maxiter=niter,disp=True
                                      ,gtol=gtol)
            elif potential in ('NN'):
                solution= opt.fmin_cg(func,vars \
                                      ,args=(maindir,) \
                                      ,fprime=NN.grad \
                                      ,maxiter=niter,disp=True \
                                      ,gtol=gtol)
        print ' CG solution:',solution

    elif fmethod in ('qn','quasi-Newtown','QN','bfgs','BFGS'):
        print '>>>>> quasi-Newton was selected.'
        if gradient in ('numerical'):
            solution= opt.fmin_bfgs(func,vars,args=(maindir,)
                                    ,maxiter=niter,disp=True
                                    ,epsilon=eps,gtol=gtol)
        else:
            if potential in ('linreg'):
                solution= opt.fmin_bfgs(func,vars,args=(maindir,)
                                        ,fprime=grad_linreg
                                        ,maxiter=niter,disp=True
                                        ,gtol=gtol)
            elif potential in ('NN'):
                solution= opt.fmin_bfgs(func,vars \
                                        ,args=(maindir,)
                                        ,fprime=NN.grad
                                        ,maxiter=niter,disp=True
                                        ,gtol=gtol)

        print ' QN solution:',solution

    elif fmethod in ('NM','Nelder-Mead','downhill-simplex'):
        print '>>>>> Nelder-Mead was selected.'
        solution= opt.fmin(func,vars,args=(maindir,)
                           ,maxiter=niter,disp=True)
        print ' NM solution:',solution
    elif fmethod in ('ga','GA','genetic-algorithm'):
        print '>>>>> genetic algorithm was selected.'
        solution= ga_wrapper()
        #...calc best one again
        func(solution,maindir)
    elif fmethod in ('sd_dynamics','SD_dynamics','SD'):
        print '>>>>> SD_dynamics was selected.'
        solution= sd_dynamics(func,vars,args=(maindir,)
                    ,fprime=grad_linreg
                    ,maxiter=niter)
    elif fmethod in ('check_grad'):
        print '>>>>> check_grad was selected.'
        if gradient in ('numerical'):
            print ' Done nothing, because gradient==numerical.'
        else:
            if potential == 'linreg':
                grad= grad_linreg(vars,maindir)
            elif potential == 'NN':
                grad= NN.grad(vars,maindir,)
            agrad= opt.approx_fprime(vars,func,eps,maindir)
            print ''
            print '>>>>> check_grad report:'
            print ' diff =',np.sqrt(np.sum((grad-agrad)**2))
            print '                grad           aprrox_grad        error (%)'
            for i in range(len(grad)):
                print ' {0:20.6f} {1:20.6f} {2:12.2f}'.format(grad[i],\
                                                             agrad[i],\
                                                             abs(grad[i]-agrad[i])/abs(grad[i])*100)
        solution= vars
    elif fmethod in ('test','TEST'):
        print '>>>>> TEST was selected.'
        #func(vars,maindir) # func is already evaluated before 
        if gradient != 'numerical':
            if potential in ('linreg'):
                grad_linreg(vars,maindir)
            elif potential in ('NN'):
                NN.grad(vars,maindir)
        solution= vars

    if not fmethod in ('test','TEST','check_grad'):
        write_params(maindir+'/'+parfile+'.fin',solution)

    if potential in ('linreg'):
        ergs,frcs= calc_ef_from_bases(solution,maindir)
    elif potential in ('NN'):
        ergs,frcs= NN.calc_ef_from_bases(solution)
    else:
        ergs,frcs= gather_pmd_data(maindir)
    output_energy_relation(ergs,ergrefs,samples,sample_dirs,fname='out.erg.pmd-vs-dft.fin')
    output_force_relation(frcs,frcrefs,samples,sample_dirs,fname='out.frc.pmd-vs-dft.fin')
    output_statistics(ergs,frcs)

    print '{0:=^72}'.format(' FITPOT finished correctly ')
    print '   Elapsed time = {0:12.2f}'.format(time.time()-_init_time)
