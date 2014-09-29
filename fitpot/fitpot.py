#!/opt/local/bin/python
# -*- coding: utf-8 -*-

u"""Fit parameters of a certain potential to DFT data.

The potential must be specified in pmd input file, in.pmd.
"""

import os,sys
import time
import glob
import numpy as np
import scipy.optimize as opt
# import matplotlib.pyplot as plt
import math

#.....import from local modules
from MD_system import MD_system
from pga import GA

#.....constants
large= 1.0e+30
tiny = 1.0e-8
max_species= 10

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
l1st_call= True
bases= []

#.....input parameters
nsmpl= 1
niter= 1
fmethod= 'test'
maindir= 'learning_set'
parfile= 'in.params.SW_Si'
runmode= 'serial'
potential= 'none'
gradient= 'none'
xtol= 1e-5
gtol= 1e-5
ftol= 1e-5
eps = 1e-8
eatom= np.zeros(max_species)
fmatch= True
penalty= 'no'
pweight= 1.0
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
    f.write(' {0:6d} {1:10.4f}\n'.format(nprms,rcut))
    for i in range(nprms):
        if abs(prange[i,0]) >= large and abs(prange[i,1]) >= large:
            f.write('{0:22.14e} \n'.format(params[i]))
        elif prange[i,0] == prange[i,1]:
            f.write('{0:22.14e} {1:22.14e}\n'.format(params[i],prange[i,0]))
        else:
            f.write('{0:22.14e} {1:22.14e} {2:22.14e}\n'.format(params[i],
                                                                prange[i,0],
                                                                prange[i,1]))
    f.close()

def params_to_vars(x,xr):
    global vars,vranges
    nvars= 0
    for i in range(nprms):
        if xr[i,0] != xr[i,1]:
            nvars += 1
    vars= np.zeros(nvars)
    vranges= np.zeros([nvars,2])
    j=0
    for i in range(nprms):
        if xr[i,0] != xr[i,1]:
            vars[j]= x[i]
            vranges[j,0]= xr[i,0]
            vranges[j,1]= xr[i,1]
            j += 1

def vars_to_params(x):
    global params
    j=0
    for i in range(nprms):
        if prange[i,0] != prange[i,1]:
            params[i]= x[j]
            j += 1

def read_input(fname='in.fitpot'):
    global nsmpl,niter,fmethod,maindir,parfile,runmode,eps \
           ,xtol,gtol,ftol,fmatch,penalty,pweight,gradient,potential
    global eatom,ga_nindv,ga_nbit,ga_temp
    dict={}
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
            #.....GA parameters
            elif data[0] == 'ga_num_individuals':
                ga_nindv= int(data[1])
            elif data[0] == 'ga_num_bit':
                ga_nbit= int(data[1])
            elif data[0] == 'ga_temperature':
                ga_temp= float(data[1])
            elif data[0] == 'ga_murate':
                ga_murate= float(data[1])
    f.close()

def show_input_params(fname='in.fitpot'):
    print '>>>>> input parameters:'
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
    return (ergs,frcs)

def gather_ref_data(basedir):
    global ergrefs,frcrefs
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

#============================================= function evaluation
def func(x,*args):
    u"""evaluate function L=sum_{samples}[E(pmd)-E(ref)]^2.
    
    This will be called from scipy.optimize.fmin_cg().
    The 1st argument x should be 1-D array of variables.
    """

    #.....write parameters to in.params.????? file
    dir= args[0]

    if fmethod in ('test','TEST') or \
            not potential in ('linreg','NN1','NN2'):
        #.....store original file
        os.system('cp '+dir+'/'+parfile+' '+dir+'/'+parfile+'.tmp')
        write_params(dir+'/'+parfile,x)
        #.....run pmd in all sample directories
        os.chdir(dir)
        #print os.getcwd(),dir
        if runmode in ('serial','Serial','SERIAL','sequential','single'):
            os.system('./serial_run_pmd.sh '+parfile)
        elif runmode in ('parallel','Parallel','PARALLEL'):
            os.system('python ./parallel_run_pmd.py '+parfile)
        else:
            print "{0:*>20}: no such run_mode !!!".format(' Error', runmode)
            exit()
        os.chdir(cwd)
        #.....restore original file
        os.system('cp '+dir+'/'+parfile+'.tmp'+' '+dir+'/'+parfile)
        #.....gather pmd results
        ergs,frcs=gather_pmd_data(dir)
    else:
        #.....calc ergs and frcs from bases data and x (variables)
        read_bases(dir)
        ergs,frcs=calc_ef_from_bases(x,dir)

    #.....calc function value of L
    val= eval_L(ergs,frcs,ergrefs,frcrefs,samples)

    print
    print ' L value=',val

    if penalty in ('ridge','Ridge','RIDGE'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            p += math.sqrt(x[n]**2)
        print ' penalty value=',p*pweight
        val += p*pweight
        print ' total L value=',val

    elif penalty in ('lasso','LASSO'):
        p= 0.0
        lx= len(x)
        for n in range(lx):
            p += abs(x[n])
        print ' penalty value=',p*pweight
        val += p*pweight
        print ' total L value=',val
    sys.stdout.flush()
    return val

def eval_L(cergs,cfrcs,rergs,rfrcs,samples):
    val= 0.0
    for i in range(len(samples)):
        val += (cergs[i]-rergs[i])**2
        if not fmatch:
            continue
        for j in range(samples[i].natm):
            for k in range(3):
                val += (cfrcs[i][j,k]-rfrcs[i][j,k])**2  \
                       /(3*samples[i].natm)
    return val

def calc_ef_from_bases(x,dir):
    """Calculate energies and forces of every samples using bases data.
    """

    #.....initialize variables
    es=np.zeros(len(samples))
    fs= []
    for smpl in samples:
        fs.append(np.zeros((smpl.natm,3)))
    
    if potential in ('linreg'):
        #.....calc energies
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            for ia in range(smpl.natm):
                for iprm in range(len(params)):
                    es[ismpl] += x[iprm] *bases[ismpl][ia,iprm]

        #.....calc forces
        if fmatch:
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                for ia in range(smpl.natm):
                    for iprm in range(len(params)):
                        dbs=bases[len(samples)+ismpl][ia,iprm]
                        fs[ismpl][ia,0] -= x[iprm] *dbs[0]
                        fs[ismpl][ia,1] -= x[iprm] *dbs[1]
                        fs[ismpl][ia,2] -= x[iprm] *dbs[2]

    elif potential in ('NN1'):
        print " NN1 is not yet implemented..."

    return (es,fs)

def read_bases(dir):
    global bases,l1st_call

    #.....read bases from pmd directories
    if l1st_call:
        if potential in ('linreg'):
            bases= gather_basis_linreg(dir)
        elif potential in ('NN1'):
            bases= gather_basis_NN1(dir)
        l1st_call= False

#======================================== derivative of linreg potential
def grad_linreg(x,*args):
    global bases
    dir= args[0]

    read_bases(dir)
    #bases= gather_basis_linreg(dir)
    #.....gather pmd results
    ergs,frcs= calc_ef_from_bases(x,dir)
    # ergs,frcs= gather_pmd_data(dir)
    
    grad= np.zeros(len(params))
    for iprm in range(len(params)):
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            ediff= ergs[ismpl] -ergrefs[ismpl]
            bs= 0.0
            for ia in range(smpl.natm):
                bs += bases[ismpl][ia,iprm]
            grad[iprm] += 2.0*ediff*bs
            
    if fmatch:
        for iprm in range(len(params)):
            dlprm= 0.0
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                dbs= np.zeros(3)
                fdiff= np.zeros(3)
                for ia in range(smpl.natm):
                    fdiff= frcs[ismpl][ia] -frcrefs[ismpl][ia]
                    dbs= bases[len(samples)+ismpl][ia,iprm]
                    dlprm -= 2.0*( fdiff[0]*dbs[0] \
                                       +fdiff[1]*dbs[1] \
                                       +fdiff[2]*dbs[2] ) \
                                       /(smpl.natm*3)
            grad[iprm] += dlprm

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
            
    return grad

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
            
#=========================================== derivative of NN1 potential
def grad_NN1(x,*args):
    dir= args[0]
    #.....read num of symmetry functions and nodes
    nsf,nhl1,wgt1,wgt2= read_NN1(dir)
    #.....check
    if len(params) != (nsf+1)*nhl1 +(nhl1+1):
        print ' [Error] len(params) != (nsf+1)*nhl1 +(nhl1+1)'
        print '   len(params)          = ',len(params)
        print '   nsf                  = ',nsf
        print '   nhl1                 = ',nhl1
        print '   (nsf+1)*nhl1+(nhl1+1)= ',(nsf+1)*nhl1 +(nhl1+1)
        exit()
    #.....gather basis data
    gsf,hl1= gather_basis_NN1(dir,nsf,nhl1)
    #.....gather pmd results
    ergs,frcs= gather_pmd_data(dir)

    grad= np.zeros(len(params))
    iprm= 0
    for isf in range(nsf+1):
        for ihl1 in range(1,nhl1+1):
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                ediff= ergs[ismpl] -ergrefs[ismpl]
                tmp= 0.0
                for ia in range(smpl.natm):
                    gsfs= gsf[ismpl]
                    hl1s= hl1[ismpl]
                    tmp += wgt2[ihl1] *hl1s[ia,ihl1] \
                           *(1.0-hl1s[ia,ihl1]) *gsfs[ia,isf]
                #print ' isf,ihl1,ismpl,ediff,tmp=',isf,ihl1,ismpl,ediff,tmp
                grad[iprm] += 2.0*ediff*tmp
            iprm += 1
    for ihl1 in range(nhl1+1):
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            ediff= ergs[ismpl] -ergrefs[ismpl]
            tmp= 0.0
            for ia in range(smpl.natm):
                tmp += hl1[ismpl][ia,ihl1]
            #print 'ihl1,ismpl,ediff,tmp=',ihl1,ismpl,ediff,tmp
            grad[iprm] += 2.0*ediff*tmp
        iprm += 1
    #print 'grad=',grad
    return grad

def gather_basis_NN1(basedir,nsf,nhl1):
    gsf= []
    hl1= []
    #...read basis data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        f= open(basedir+'/'+dir+'/pmd/out.gsf','r')
        g= open(basedir+'/'+dir+'/pmd/out.hl1','r')
        #.....skip 1st line
        dataf= f.readline().split()
        datag= g.readline().split()
        gsfs= np.zeros((smpl.natm,nsf+1))
        hl1s= np.zeros((smpl.natm,nhl1+1))
        for ia in range(smpl.natm):
            for isf in range(nsf+1):
                dataf= f.readline().split()
                gsfs[ia,isf]= float(dataf[2])
        for ia in range(smpl.natm):
            for ihl1 in range(nhl1+1):
                datag= g.readline().split()
                hl1s[ia,ihl1]= float(datag[2])
        gsf.append(gsfs)
        hl1.append(hl1s)
    return gsf,hl1

def read_NN1(basedir):
    f= open(basedir+'/'+'in.const.NN1')
    g= open(basedir+'/'+'in.params.NN1')
    data= f.readline().split()
    nsf=  int(data[0])
    nhl1= int(data[1])
    data= g.readline().split()
    wgt1= np.zeros((nsf+1,nhl1))
    wgt2= np.zeros(nhl1+1)
    for isf in range(nsf+1):
        for ihl1 in range(nhl1):
            data= g.readline().split()
            wgt1[isf,ihl1]= float(data[0])
    for ihl1 in range(nhl1+1):
        data= g.readline().split()
        wgt2[ihl1]= float(data[0])
    f.close()
    g.close()
    return nsf,nhl1,wgt1,wgt2

#========================================================== output data
def grad_NN2(x,*args):
    dir= args[0]
    #.....read num of symmetry functions and nodes
    nsf,nhl1,nhl2,wgt1,wgt2,wgt3= read_NN2(dir)
    #.....check
    if len(params) != (nsf+1)*nhl1 +(nhl1+1)*nhl2 +(nhl2+1):
        print ' [Error] len(params) != (nsf+1)*nhl1 +(nhl1+1)*nhl2 +(nhl2+1)'
        print '   len(params) = ',len(params)
        print '   nsf         = ',nsf
        print '   nhl1        = ',nhl1
        print '   nhl2        = ',nhl2
        print '   nparams     = ',(nsf+1)*nhl1 +(nhl1+1)*nhl2 +(nhl2+1)
        exit()
    #.....gather basis data
    gsf,hl1,hl2= gather_basis_NN2(dir,nsf,nhl1,nhl2)
    #.....gather pmd results
    ergs,frcs= gather_pmd_data(dir)

    grad= np.zeros(len(params))
    iprm= 0
    for isf in range(nsf+1):
        for ihl1 in range(1,nhl1+1):
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                ediff= ergs[ismpl] -ergrefs[ismpl]
                tmp= 0.0
                gsfs= gsf[ismpl]
                hl1s= hl1[ismpl]
                hl2s= hl2[ismpl]
                for ia in range(smpl.natm):
                    hl1i= hl1s[ia,ihl1]
                    tmp2= hl1i*(1.0-hl1i) *gsfs[ia,isf]
                    for ihl2 in range(1,nhl2+1):
                        hl2i= hl2s[ia,ihl2]
                        tmp += wgt3[ihl2] *hl2i *(1.0-hl2i) \
                            *wgt2[ihl1,ihl2] *tmp2
                grad[iprm] += 2.0*ediff*tmp
            iprm += 1
    for ihl1 in range(nhl1+1):
        for ihl2 in range(1,nhl2+1):
            for ismpl in range(len(samples)):
                smpl= samples[ismpl]
                ediff= ergs[ismpl] -ergrefs[ismpl]
                tmp= 0.0
                hl1s= hl1[ismpl]
                hl2s= hl2[ismpl]
                for ia in range(smpl.natm):
                    hl1i= hl1s[ia,ihl1]
                    hl2i= hl2s[ia,ihl2]
                    tmp += wgt3[ihl2] *hl2i *(1.0-hl2i) *hl1i
            #print 'ihl1,ismpl,ediff,tmp=',ihl1,ismpl,ediff,tmp
                grad[iprm] += 2.0*ediff*tmp
            iprm += 1
    for ihl2 in range(nhl2+1):
        for ismpl in range(len(samples)):
            smpl= samples[ismpl]
            ediff= ergs[ismpl] -ergrefs[ismpl]
            tmp= 0.0
            hl2s= hl2[ismpl]
            for ia in range(smpl.natm):
                hl2i= hl2s[ia,ihl2]
                tmp += hl2i
            grad[iprm] += 2.0*ediff*tmp
        iprm += 1
    #print grad
    return grad

def gather_basis_NN2(basedir,nsf,nhl1,nhl2):
    gsf= []
    hl1= []
    hl2= []
    #...read basis data
    for i in range(len(sample_dirs)):
        dir= sample_dirs[i]
        smpl= samples[i]
        f= open(basedir+'/'+dir+'/pmd/out.gsf','r')
        g= open(basedir+'/'+dir+'/pmd/out.hl1','r')
        h= open(basedir+'/'+dir+'/pmd/out.hl2','r')
        #.....skip 1st line
        dataf= f.readline().split()
        datag= g.readline().split()
        datah= h.readline().split()
        gsfs= np.zeros((smpl.natm,nsf+1))
        hl1s= np.zeros((smpl.natm,nhl1+1))
        hl2s= np.zeros((smpl.natm,nhl2+1))
        for ia in range(smpl.natm):
            for isf in range(nsf+1):
                dataf= f.readline().split()
                gsfs[ia,isf]= float(dataf[2])
        for ia in range(smpl.natm):
            for ihl1 in range(nhl1+1):
                datag= g.readline().split()
                hl1s[ia,ihl1]= float(datag[2])
        for ia in range(smpl.natm):
            for ihl2 in range(nhl2+1):
                datah= h.readline().split()
                hl2s[ia,ihl2]= float(datah[2])
        gsf.append(gsfs)
        hl1.append(hl1s)
        hl2.append(hl2s)
    return gsf,hl1,hl2

def read_NN2(basedir):
    f= open(basedir+'/'+'in.const.NN2')
    g= open(basedir+'/'+'in.params.NN2')
    data= f.readline().split()
    nsf=  int(data[0])
    nhl1= int(data[1])
    nhl2= int(data[2])
    data= g.readline().split()
    wgt1= np.zeros((nsf+1,nhl1))
    wgt2= np.zeros((nhl1+1,nhl2+1))
    wgt3= np.zeros(nhl2+1)
    for isf in range(nsf+1):
        for ihl1 in range(nhl1):
            data= g.readline().split()
            wgt1[isf,ihl1]= float(data[0])
    for ihl1 in range(nhl1+1):
        for ihl2 in range(1,nhl2+1):
            data= g.readline().split()
            wgt2[ihl1,ihl2]= float(data[0])
    for ihl2 in range(nhl2+1):
        data= g.readline().split()
        wgt3[ihl2]= float(data[0])
    f.close()
    g.close()
    return nsf,nhl1,nhl2,wgt1,wgt2,wgt3

#========================================================== output data
def output_energy_relation(ergs,fname='out.erg.pmd-vs-dft'):
    f= open(fname,'w')
    for i in range(len(ergrefs)):
        smpl= samples[i]
        f.write(' {0:15.7e} {1:15.7e}\n'.format(ergrefs[i]/smpl.natm \
                                              ,ergs[i]/smpl.natm ))
    f.close()

def output_force_relation(frcs,fname='out.frc.pmd-vs-dft'):
    f= open(fname,'w')
    for i in range(len(samples)):
        for j in range(samples[i].natm):
            for k in range(3):
                f.write(' {0:15.7e} {1:15.7e}\n'.format(frcrefs[i][j,k],frcs[i][j,k]))
    f.close()

# def plot_energy_relation(ergs,fname='graph.eps'):
#     x=ergrefs
#     y=ergs
#     plt.scatter(x,y)
#     xmin,xmax= plt.xlim()
#     ymin,ymax= plt.ylim()
#     xmax= max(xmax,ymax)
#     xmin= min(xmin,ymin)
#     dat= np.arange(xmin,xmax,50.0)
#     plt.plot(dat,dat,'b--')
#     plt.xlabel('DFT energy (eV)')
#     plt.ylabel('pmd energy (eV)')
#     plt.title('DFT and pmd energy relation')
#     plt.xlim(xmin,xmax)
#     plt.ylim(xmin,xmax)
#     plt.savefig(fname)
#     #plt.show()

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
    #print ' final variabls:'
    #print x
    return x
    
#============================================= main routine hereafter
if __name__ == '__main__':
    print "{0:=^72}".format(' FITPOT ')
    t0= time.time()
    cwd= os.getcwd()
    #.....inputs: parameters in in.fitpot as a dictionary
    inputs= read_input('in.fitpot')
    show_input_params('in.fitpot')
    #.....params: parameters in in.params.?????
    read_params(maindir+'/'+parfile)
    write_params(maindir+'/'+parfile+'.ini',vars)
    
    #.....get samples from ##### directories
    sample_dirs= get_sample_dirs()
    sample_dirs.sort()
    if nsmpl != len(sample_dirs):
        print '{0:*>20}: num_samples in in.fitpot is wrong.'.format(' Error')
        exit()
    read_pos()

    #.....initial data
    gather_ref_data(maindir)
    #.....1st call of func
    func(vars,maindir)
    if potential in ('linreg'):
        ergs,frcs= calc_ef_from_bases(vars,maindir)
    else:
        ergs,frcs= gather_pmd_data(maindir)

    output_energy_relation(ergs,fname='out.erg.pmd-vs-dft.ini')
    output_force_relation(frcs,fname='out.frc.pmd-vs-dft.ini')

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
            elif potential in ('NN1'):
                solution= opt.fmin_cg(func,vars,args=(maindir,)
                                      ,fprime=grad_NN1
                                      ,maxiter=niter,disp=True
                                      ,gtol=gtol)
            elif potential in ('NN2'):
                solution= opt.fmin_cg(func,vars,args=(maindir,)
                                      ,fprime=grad_NN2
                                      ,maxiter=niter,disp=True
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
            elif potential in ('NN1'):
                solution= opt.fmin_bfgs(func,vars,args=(maindir,)
                                        ,fprime=grad_NN1
                                        ,maxiter=niter,disp=True
                                        ,gtol=gtol)
            elif potential in ('grad_NN2','NN2'):
                solution= opt.fmin_bfgs(func,vars,args=(maindir,)
                                        ,fprime=grad_NN2
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
            elif potential == 'NN1':
                grad= grad_NN1(vars,maindir)
            elif potential == 'NN2':
                grad= grad_NN2(vars,maindir)
            agrad= opt.approx_fprime(vars,func,eps,maindir)
            print ' diff =',np.sqrt(np.sum((grad-agrad)**2))
            print ' grad       aprrox_grad'
            for i in range(len(grad)):
                print ' {0:12.6f} {1:12.6f}'.format(grad[i],agrad[i])
        solution= vars
    elif fmethod in ('test','TEST'):
        print '>>>>> TEST was selected.'
        func(vars,maindir)
        if gradient != 'numerical':
            if potential in ('linreg'):
                grad_linreg(vars,maindir)
            elif potential in ('NN1'):
                grad_NN1(vars,maindir)
            elif potential in ('NN2'):
                grad_NN2(vars,maindir)
        solution= vars

    write_params(maindir+'/'+parfile+'.fin',solution)

    if potential in ('linreg'):
        ergs,frcs= calc_ef_from_bases(vars,maindir)
    else:
        ergs,frcs= gather_pmd_data(maindir)
    output_energy_relation(ergs,fname='out.erg.pmd-vs-dft.fin')
    output_force_relation(frcs,fname='out.frc.pmd-vs-dft.fin')

    print '{0:=^72}'.format(' FITPOT finished correctly ')
    print '   Elapsed time = {0:12.2f}'.format(time.time()-t0)
