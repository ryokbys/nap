#  Time-stamp: <2015-02-12 17:46:35 Ryo KOBAYASHI>
"""
Routines related to neural network with one hidden layer.
"""

import math,time,os
import numpy as np
import multiprocessing as mp

#.....import from local modules
import fitpot

#.....constants which may be different from those in fitpot.py
_large= 1.0e+30

#.....global variables
_l1st= True
_nsp=  1
_nsf=  0
_nhl1= 0
_ergs= []
_frcs= []
_gsf= []
_dgsf=[]
_hl1= []
_aml= []
_bml= []
_fmatch= False
_basedir='learning_set'
_samples=[]
_sample_dirs=[]
_nprcs= 1
_ergrefs= []
_frcrefs= []
_fmethod= 'test'
_parfile= 'in.params.NN1'
_runmode= 'serial'
_rcut  = 5.0
_params= []
_pranges= []

#============================================================ routines
def init(*args,**kwargs):
    """
    Initialize variables for the calculation of NN1 parameters.
    This should be called at the first place
    before any other NN1-related routines.
    """

    global _nsf,_nhl1,_wgt1,_wgt2,_gsf,_dgsf,_hl1,_aml,_bml
    global _fmatch,_basedir,_samples,_sample_dirs,_nprcs
    global _ergrefs,_frcrefs,_fmethod,_parfile,_runmode,_rcut
    global _params,_pranges,_vranges

    _basedir = args[0]
    _params   = args[1]
    _sample_dirs= args[2]
    _samples = args[3]
    _nprcs   = args[4]
    _fmatch  = args[5]
    _ergrefs = args[6]
    _frcrefs = args[7]
    _fmethod = args[8]
    _parfile = args[9]
    _runmode = args[10]
    _rcut    = args[11]
    _pranges = args[12]
    _vranges = args[13]

    #.....read NN1 parameters
    f= open(_basedir+'/'+'in.const.NN1')
    data= f.readline().split()
    _nsfc= int(data[0])
    _nhl1= int(data[1])
    _nsp = int(data[2])
    n2= 0
    n3= 0
    for line in f.readlines():
        if int(line.split()[0]) == 1:
            n2 += 1
        elif int(line.split()[0]) == 2:
            n3 += 1
    f.close()
    print ' Number of nodes in hidden-layer-1   =',_nhl1
    print ' Number of species                   =',_nsp
    print ' Number of 2-body symmetry functions =',n2
    print ' Number of 3-body symmetry functions =',n3
    
    # check
    ncmb2= _nsp +factorial(_nsp,2)/2
    ncmb3= _nsp *ncmb2
    _nsf= n2*ncmb2 +n3*ncmb3
    print ' Number of symmetry functions        =',_nsf
    if len(_params) != _nsf*_nhl1 +_nhl1:
        print ' [Error] len(params) != (nsf+1)*nhl1 +(nhl1+1)'
        print '   len(params)          = ',len(_params)
        print '   nsf                  = ',_nsf
        print '   nhl1                 = ',_nhl1
        print '   nsf*nhl1 +nhl1       = ',_nsf*_nhl1 +_nhl1
        exit()

    #.....read bases
    #_gsf,_hl1,_aml,_bml= gather_basis(*args)
    _gsf,_dgsf= gather_bases_new(*args)

def factorial(n,m):
    """
    Returns factorial of n by m-times.
    """
    if m <= 0:
        return 1
    return n*factorial(n-1,m-1)

def sigmoid(x):
    if x < -10.0:
        return 0.0
    elif x > 10.0:
        return 1.0
    return 1.0/(1.0 +math.exp(-x))

def vars2wgts(x):
    wgt1= np.zeros((_nsf+1,_nhl1+1))
    wgt2= np.zeros(_nhl1+1)
    ix= 0
    for isf in range(1,_nsf+1):
        for ihl1 in range(1,_nhl1+1):
            wgt1[isf,ihl1]= x[ix]
            ix += 1
    for ihl1 in range(1,_nhl1+1):
        wgt2[ihl1]= x[ix]
        ix += 1
    return (wgt1,wgt2)

def calc_ef_from_pmd(x,*args):
    dir= args[0]
    cwd= os.getcwd()
    # print ' cwd=',cwd
    # print ' dir=',dir
    #.....store original file
    os.system('cp '+dir+'/'+_parfile+' '+dir+'/'+_parfile+'.tmp')
    write_params(dir+'/'+_parfile,x)
    #.....run pmd in all sample directories
    os.chdir(dir)
    #print os.getcwd(),dir
    if _runmode in ('serial','Serial','SERIAL','sequential','single'):
        os.system('./serial_run_pmd.sh '+_parfile)
    elif _runmode in ('parallel','Parallel','PARALLEL'):
        os.system('python ./parallel_run_pmd.py '+_parfile)
    else:
        print "{0:*>20}: no such run_mode !!!".format(' Error', _runmode)
        exit()
    os.chdir(cwd)
    #.....restore original file
    os.system('cp '+dir+'/'+_parfile+' '+dir+'/'+_parfile+'.current')
    os.system('cp '+dir+'/'+_parfile+'.tmp'+' '+dir+'/'+_parfile)
    #.....gather pmd results
    ergs,frcs= gather_pmd_data(dir)
    return (ergs,frcs)

def calc_ef_from_smd(x,*args):
    dir= args[0]
    cwd= os.getcwd()
    # print ' cwd=',cwd
    # print ' dir=',dir
    #.....store original file
    os.system('cp '+dir+'/'+_parfile+' '+dir+'/'+_parfile+'.tmp')
    write_params(dir+'/'+_parfile,x)
    #.....run smd in all sample directories
    os.chdir(dir)
    #print os.getcwd(),dir
    if _runmode in ('serial','Serial','SERIAL','sequential','single'):
        os.system('./serial_run_smd.sh '+_parfile)
    elif _runmode in ('parallel','Parallel','PARALLEL'):
        os.system('python ./parallel_run_smd.py '+_parfile)
    else:
        print "{0:*>20}: no such run_mode !!!".format(' Error', _runmode)
        exit()
    os.chdir(cwd)
    #.....restore original file
    os.system('cp '+dir+'/'+_parfile+' '+dir+'/'+_parfile+'.current')
    os.system('cp '+dir+'/'+_parfile+'.tmp'+' '+dir+'/'+_parfile)
    #.....gather smd results
    ergs,frcs= gather_smd_data(dir)
    return (ergs,frcs)

def write_params(fname,x):
    
    params,pranges= fitpot.vars_to_params(x,_vranges,_params,_pranges)
    f=open(fname,'w')
    f.write(' {0:6d} {1:10.4f}\n'.format(len(params),_rcut))
    for i in range(len(params)):
        if abs(pranges[i,0]) >= _large and abs(pranges[i,1]) >= _large:
            f.write('{0:22.14e} \n'.format(params[i]))
        elif pranges[i,0] == pranges[i,1]:
            f.write('{0:22.14e} {1:22.14e}\n'.format(params[i],pranges[i,0]))
        else:
            f.write('{0:22.14e} {1:22.14e} {2:22.14e}\n'.format(params[i],
                                                                pranges[i,0],
                                                                pranges[i,1]))
    f.close()

def gather_pmd_data(basedir):
    #.....initialize variables
    ergs=np.zeros(len(_samples))
    frcs= []
    for smpl in _samples:
        frcs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(_sample_dirs)):
        dir= _sample_dirs[i]
        smpl= _samples[i]
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

def gather_smd_data(basedir):
    #.....initialize variables
    ergs=np.zeros(len(_samples))
    frcs= []
    for smpl in _samples:
        frcs.append(np.zeros((smpl.natm,3)))
    #.....read data
    for i in range(len(_sample_dirs)):
        dir= _sample_dirs[i]
        smpl= _samples[i]
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
    return (ergs,frcs)


def calc_ef_from_bases(x,*args):
    """
    Calculate energies and forces of every samples using bases data.
    """
    global _hl1,_l1st,_gsf,_dgsf,_ergs,_frcs,_wgt1,_wgt2

    #.....initialize variables
    _wgt1,_wgt2= vars2wgts(x)
    es=np.zeros(len(_samples))
    fs= []
    for smpl in _samples:
        fs.append(np.zeros((smpl.natm,3)))

    # #.....gather bases from file if this is 1st call
    # if _l1st:
    #     _gsf,_dgsf= gather_bases_new()
    #     _l1st = False

    p= mp.Pool(_nprcs)
    _hl1= []
    if _nprcs == 1:
        for ismpl in range(len(_samples)):
            smpl= _samples[ismpl]
            est,fst,hl1s= calc_ef(ismpl,x,*args)
            _hl1.append(hl1s)
            es[ismpl]= est
            for ia in range(smpl.natm):
                fs[ismpl][ia,0] += fst[ia,0]
                fs[ismpl][ia,1] += fst[ia,1]
                fs[ismpl][ia,2] += fst[ia,2]
    else:
        func_args=[]
        for ismpl in range(len(_samples)):
            func_args.append( (calc_ef,ismpl,x) )
        results= p.map(arg_wrapper,func_args)
        p.close()
        p.join()
        for ismpl in range(len(_samples)):
            smpl= _samples[ismpl]
            est,fst,hl1s= results[ismpl]
            _hl1.append(hl1s)
            es[ismpl]= est
            for ia in range(smpl.natm):
                fs[ismpl][ia,0] += fst[ia,0]
                fs[ismpl][ia,1] += fst[ia,1]
                fs[ismpl][ia,2] += fst[ia,2]

    # print ' es:'
    # print es
    _ergs= es
    _frcs= fs
    return (es,fs)

def arg_wrapper(args):
    return args[0](*args[1:])

def calc_ef(ismpl,x,*args):
    global _wgt1,_wgt2

    smpl= _samples[ismpl]
    gsfs= _gsf[ismpl]
    dgsfs= _dgsf[ismpl]
    es= 0.0
    fs= np.zeros((smpl.natm,3))
    #.....energy
    iprm= 0
    #print ' ismpl=',ismpl
    hl1s= np.zeros((_nhl1+1,smpl.natm))
    for ia in range(smpl.natm):
        #.....calc hidden-layer value using sigmoid function
        tmp= 0.0
        for ihl1 in range(1,_nhl1+1):
            for isf in range(1,_nsf+1):
                hl1s[ihl1,ia]+= _wgt1[isf,ihl1] *gsfs[ia,isf]
            hl1s[ihl1,ia]= sigmoid(hl1s[ihl1,ia])
        for ihl1 in range(1,_nhl1+1):
            es += _wgt2[ihl1] *(hl1s[ihl1,ia]-0.5)

    #.....forces
    if _fmatch:
        dhdg= np.zeros(3)
        for ihl1 in range(1,_nhl1+1):
            w2= _wgt2[ihl1]
            for ja in range(smpl.natm):
                dh=hl1s[ihl1,ja]*(1.0-hl1s[ihl1,ja])
                for isf in range(1,_nsf+1):
                    w1= _wgt1[isf,ihl1]
                    for ia in range(smpl.natm):
                        dhdg[:]=  dh*dgsfs[ja,isf,ia,:]
                        fs[ia,:] -= w1*w2*dhdg[:]
    return (es,fs,hl1s)

def grad(x,*args):
    global _wgt1,_wgt2,_gsf,_hl1,_aml,_bml
    
    t0= time.time()
    dir= args[0]
    #.....get energies and forces
    #ergs,frcs= calc_ef_from_bases(x,*args)
    #ergs,frcs= calc_ef_from_pmd(x,*args)
    #ergs,frcs= gather_pmd_data(dir)
    # for ismpl in range(len(_samples)):
    #     print ' ismpl,ergs,_ergrefs:',ismpl,ergs[ismpl],_ergrefs[ismpl]
    
    #_gsf,_hl1,_aml,_bml= gather_basis(*args)
    #_wgt1,_wgt2= vars2wgts(x)

    p= mp.Pool(_nprcs)
    grad= np.zeros(len(x))
    
    if _nprcs == 1:
        for ismpl in range(len(_samples)):
            gs= grad_core_new(ismpl,x)
            for iprm in range(len(x)):
                grad[iprm] += gs[iprm]
    else:
        func_args=[]
        for ismpl in range(len(_samples)):
            func_args.append( (grad_core_new,ismpl,x))
        results= p.map(arg_wrapper,func_args)
        p.close()
        p.join()
        for ismpl in range(len(_samples)):
            gs= results[ismpl]
            for iprm in range(len(x)):
                grad[iprm] += gs[iprm]

    print ' ===> time NN1.grad: {0:12.3f} sec'.format(time.time()-t0)
    #print ' grad=',grad
    return grad

def grad_core(ismpl,ergs,frcs,*args):
    x      = args[0]

    gs= np.zeros(len(x))
    dgs=np.zeros(len(x))
    smpl= _samples[ismpl]
    ediff= (ergs[ismpl] -_ergrefs[ismpl]) /smpl.natm
    gsfs= _gsf[ismpl]
    hl1s= _hl1[ismpl]
    iprm= 0
    #print ' ismpl=',ismpl
    for isf in range(1,_nsf+1):
        for ihl1 in range(1,_nhl1+1):
            tmp= 0.0
            for ia in range(smpl.natm):
                tmp += _wgt2[ihl1] *hl1s[ia,ihl1] \
                         *(1.0 -hl1s[ia,ihl1]) *gsfs[ia,isf]
            gs[iprm] += 2.0*ediff*tmp 
            iprm += 1
    for ihl1 in range(1,_nhl1+1):
        tmp= 0.0
        for ia in range(smpl.natm):
            tmp += (hl1s[ia,ihl1] -0.5)
        gs[iprm] += 2.0*ediff*tmp
        iprm += 1
#     for ia in range(smpl.natm):
#         for ihl1 in range(_nhl1+1):
#             print '   ia,ihl1,hl1s=',ia,ihl1,hl1s[ia,ihl1]

    if _fmatch:
        amls= _aml[ismpl]
        bmls= _bml[ismpl]
        iprm= 0
        #print ' ismpl=',ismpl
        for isf in range(1,_nsf+1):
            for ihl1 in range(1,_nhl1+1):
                tmp= 0.0
                w2= _wgt2[ihl1]
                for ia in range(smpl.natm):
                    fdiff= frcs[ismpl][ia] -_frcrefs[ismpl][ia]
                    am= amls[ia,ihl1,isf]
                    bm= bmls[ia,ihl1,isf]
                    #print '   isf,ihl1,ia,am,bm=',isf,ihl1,ia,am[:],bm[:]
                    tmp -= 2.0*w2*( fdiff[0]*(am[0]+bm[0]) \
                                    +fdiff[1]*(am[1]+bm[1]) \
                                    +fdiff[2]*(am[2]+bm[2]) ) /smpl.natm/3
                dgs[iprm] += tmp
                iprm += 1
        for ihl1 in range(1,_nhl1+1):
            tmp= 0.0
            for ia in range(smpl.natm):
                fdiff= frcs[ismpl][ia] -_frcrefs[ismpl][ia]
                for isf in range(1,_nsf+1):
                    am= amls[ia,ihl1,isf]
                    w1= _wgt1[isf,ihl1]
                    tmp -= 2.0*w1*( fdiff[0]*am[0] \
                                    +fdiff[1]*am[1] \
                                    +fdiff[2]*am[2] ) /smpl.natm/3
            dgs[iprm] += tmp
            iprm += 1
    # print ' gs,dgs,gs+dgs:'
    # for iprm in range(len(x)):
    #     print ' {0:15.7f} {1:15.7f} {2:15.7f}'.format(gs[iprm],dgs[iprm],gs[iprm]+dgs[iprm])
    gs= gs +dgs
    return gs
                    
def grad_core_new(ismpl,*args):
    x      = args[0]

    gs= np.zeros(len(x))
    dgs=np.zeros(len(x))
    smpl= _samples[ismpl]
    ediff= (_ergs[ismpl] -_ergrefs[ismpl]) /smpl.natm
    gsfs= _gsf[ismpl]
    dgsfs= _dgsf[ismpl]
    hl1s= _hl1[ismpl]
    iprm= _nsf*_nhl1 +_nhl1
    for ihl1 in range(_nhl1,0,-1):
        tmp= 0.0
        for ia in range(smpl.natm):
            tmp += (hl1s[ihl1,ia] -0.5)
        iprm -= 1
        gs[iprm] += 2.0*ediff*tmp
    for isf in range(_nsf,0,-1):
        for ihl1 in range(_nhl1,0,-1):
            tmp= 0.0
            for ia in range(smpl.natm):
                h= hl1s[ihl1,ia]
                tmp += _wgt2[ihl1] *h*(1.0-h) *gsfs[ia,isf]
            iprm -= 1
            gs[iprm] += 2.0*ediff*tmp 

    if _fmatch:
        am= np.zeros((_nsf+1,_nhl1+1,smpl.natm,3))
        bm= np.zeros((_nsf+1,_nhl1+1,smpl.natm,3))
        cm= np.zeros(3)
        iprm= _nsf*_nhl1 +_nhl1
        for ihl1 in range(_nhl1,0,-1):
            tmp= 0.0
            for ja in range(smpl.natm):
                h= hl1s[ihl1,ja]
                dh= h*(1.0-h)
                ddh= h*(1.0-h)*(1.0-2.0*h)
                for isf in range(1,_nsf+1):
                    w1= _wgt1[isf,ihl1]
                    ddhg= ddh*gsfs[ja,isf]
                    for ia in range(smpl.natm):
                        fdiff= (_frcs[ismpl][ia] -_frcrefs[ismpl][ia]) \
                               *2/smpl.natm/3
                        am[isf,ihl1,ia,:]+= dh*dgsfs[ja,isf,ia,:]
                        bm[isf,ihl1,ia,:]+= ddhg*dgsfs[ja,isf,ia,:]
                        tmp -= w1*( fdiff[0]*dh*dgsfs[ja,isf,ia,0] \
                                    +fdiff[1]*dh*dgsfs[ja,isf,ia,1] \
                                    +fdiff[2]*dh*dgsfs[ja,isf,ia,2] )
            iprm -= 1
            dgs[iprm] += tmp
        for isf in range(_nsf,0,-1):
            for ihl1 in range(_nhl1,0,-1):
                tmp= 0.0
                w2= _wgt2[ihl1]
                for ia in range(smpl.natm):
                    fdiff= (_frcs[ismpl][ia] -_frcrefs[ismpl][ia]) \
                        *2 /smpl.natm/3
                    cm[:]= am[isf,ihl1,ia,:] +bm[isf,ihl1,ia,:]
                    tmp -= w2*( fdiff[0]*(cm[0]) \
                                    +fdiff[1]*(cm[1]) \
                                    +fdiff[2]*(cm[2]) )
                iprm -= 1
                dgs[iprm] += tmp
    gs[:]= gs[:] +dgs[:]
    return gs
                    
def gather_basis(*args):
    gsf= []
    hl1= []
    aml= []
    bml= []
    #...read basis data
    for i in range(len(_sample_dirs)):
        dir= _sample_dirs[i]
        smpl= _samples[i]
        f1= open(_basedir+'/'+dir+'/pmd/out.NN1.gsf','r')
        f2= open(_basedir+'/'+dir+'/pmd/out.NN1.hl1','r')
        f3= open(_basedir+'/'+dir+'/pmd/out.NN1.aml','r')
        f4= open(_basedir+'/'+dir+'/pmd/out.NN1.bml','r')
        #.....skip 1st line
        data1= f1.readline().split()
        data2= f2.readline().split()
        data3= f3.readline().split()
        data4= f4.readline().split()
        gsfs= np.zeros((smpl.natm,_nsf+1))
        hl1s= np.zeros((smpl.natm,_nhl1+1))
        amls= np.zeros((smpl.natm,_nhl1+1,_nsf+1,3))
        bmls= np.zeros((smpl.natm,_nhl1+1,_nsf+1,3))
        for ia in range(smpl.natm):
            for isf in range(1,_nsf+1):
                data1= f1.readline().split()
                gsfs[ia,isf]= float(data1[2])
        for ia in range(smpl.natm):
            for ihl1 in range(1,_nhl1+1):
                data2= f2.readline().split()
                hl1s[ia,ihl1]= float(data2[2])
        for ia in range(smpl.natm):
            for ihl1 in range(1,_nhl1+1):
                for isf in range(1,_nsf+1):
                    data3= f3.readline().split()
                    amls[ia,ihl1,isf,0]= float(data3[3])
                    amls[ia,ihl1,isf,1]= float(data3[4])
                    amls[ia,ihl1,isf,2]= float(data3[5])
                    data4= f4.readline().split()
                    bmls[ia,ihl1,isf,0]= float(data4[3])
                    bmls[ia,ihl1,isf,1]= float(data4[4])
                    bmls[ia,ihl1,isf,2]= float(data4[5])
        gsf.append(gsfs)
        hl1.append(hl1s)
        aml.append(amls)
        bml.append(bmls)
        f1.close()
        f2.close()
        f3.close()
        f4.close()

    # print ' hl1:'
    # for ismpl in range(len(_samples)):
    #     smpl= _samples[ismpl]
    #     for ia in range(smpl.natm):
    #         for ihl1 in range(_nhl1+1):
    #             print ' ismpl,ia,ihl1,hl1=',ismpl,ia,ihl1,hl1[ismpl][ia,ihl1]

    
    return gsf,hl1,aml,bml

def gather_bases_new(*args):
    gsf= []
    dgsf=[]
    #...read basis data
    for i in range(len(_sample_dirs)):
        dir= _sample_dirs[i]
        smpl= _samples[i]
        f1= open(_basedir+'/'+dir+'/smd/out.NN1.gsf','r')
        f2= open(_basedir+'/'+dir+'/smd/out.NN1.dgsf','r')
        #.....skip 1st line
        data1= f1.readline().split()
        #data2= f2.readline().split()
        gsfs= np.zeros((smpl.natm,_nsf+1))
        dgsfs=np.zeros((smpl.natm,_nsf+1,smpl.natm,3))
        for ia in range(smpl.natm):
            for isf in range(1,_nsf+1):
                data1= f1.readline().split()
                gsfs[ia,isf]= float(data1[2])
        for ia in range(smpl.natm):
            for isf in range(1,_nsf+1):
                for ja in range(smpl.natm):
                    data2= f2.readline().split()
                    dgsfs[ia,isf,ja,0]= float(data2[3])
                    dgsfs[ia,isf,ja,1]= float(data2[4])
                    dgsfs[ia,isf,ja,2]= float(data2[5])
        gsf.append(gsfs)
        dgsf.append(dgsfs)
        f1.close()
        f2.close()

    # print ' hl1:'
    # for ismpl in range(len(_samples)):
    #     smpl= _samples[ismpl]
    #     for ia in range(smpl.natm):
    #         for ihl1 in range(_nhl1+1):
    #             print ' ismpl,ia,ihl1,hl1=',ismpl,ia,ihl1,hl1[ismpl][ia,ihl1]

    
    return gsf,dgsf
