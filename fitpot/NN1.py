"""
Routines related to neural network with one hidden layer.
"""

import math,time
import numpy as np
import multiprocessing as mp

#.....global variables
_nsf=  0
_nhl1= 0
_gsf= []
_hl1= []
_aml= []
_fmatch= False
_basedir='learning_set'
_samples=[]
_sample_dirs=[]
_nprcs= 1
_ergrefs= []
_frcrefs= []

#============================================================ routines
def init(*args,**kwargs):
    """
    Initialize variables for the calculation of NN1 parameters.
    This should be called at the first place
    before any other NN1-related routines.
    """
    global _nsf,_nhl1,_wgt1,_wgt2,_gsf,_hl1,_aml
    global _fmatch,_basedir,_samples,_sample_dirs,_nprcs
    global _ergrefs,_frcrefs

    _basedir = args[0]
    params   = args[1]
    _sample_dirs= args[2]
    _samples = args[3]
    _nprcs   = args[4]
    _fmatch  = args[5]
    _ergrefs = args[6]
    _frcrefs = args[7]

    #.....read NN1 parameters
    f= open(_basedir+'/'+'in.const.NN1')
    data= f.readline().split()
    _nsf = int(data[0])
    _nhl1= int(data[1])
    f.close()

    # check
    if len(params) != (_nsf+1)*_nhl1 +(_nhl1+1):
        print ' [Error] len(params) != (nsf+1)*nhl1 +(nhl1+1)'
        print '   len(params)          = ',len(params)
        print '   nsf                  = ',_nsf
        print '   nhl1                 = ',_nhl1
        print '   (nsf+1)*nhl1+(nhl1+1)= ',(_nsf+1)*_nhl1 +(_nhl1+1)
        exit()

    #.....read bases
    _gsf,_hl1,_aml= gather_basis(*args)

def sigmoid(x):
    return 1.0/(1.0 +math.exp(-x))

def vars2wgts(x):
    wgt1= np.zeros((_nsf+1,_nhl1+1))
    wgt2= np.zeros(_nhl1+1)
    ix= 0
    for isf in range(_nsf+1):
        for ihl1 in range(1,_nhl1+1):
            wgt1[isf,ihl1]= x[ix]
            ix += 1
    for ihl1 in range(_nhl1+1):
        wgt2[ihl1]= x[ix]
        ix += 1
    return (wgt1,wgt2)

def calc_ef_from_bases(x,*args):
    """
    Calculate energies and forces of every samples using bases data.
    """

    #.....initialize variables
    es=np.zeros(len(_samples))
    fs= []
    for smpl in _samples:
        fs.append(np.zeros((smpl.natm,3)))

    p= mp.Pool(_nprcs)
    
    if _nprcs == 1:
        for ismpl in range(len(_samples)):
            smpl= _samples[ismpl]
            est,fst= calc_ef(ismpl,x,*args)
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
        for ismpl in range(len(_samples)):
            smpl= _samples[ismpl]
            est,fst= results[ismpl]
            es[ismpl]= est
            for ia in range(smpl.natm):
                fs[ismpl][ia,0] += fst[ia,0]
                fs[ismpl][ia,1] += fst[ia,1]
                fs[ismpl][ia,2] += fst[ia,2]

    return (es,fs)


def calc_ef(ismpl,x,*args):
    global _wgt1,_wgt2

    smpl= _samples[ismpl]
    gsfs= _gsf[ismpl]
    es= 0.0
    fs= np.zeros((smpl.natm,3))
    _wgt1,_wgt2= vars2wgts(x)
    #.....energy
    iprm= 0
    #print ' ismpl=',ismpl
    for ia in range(smpl.natm):
        #.....calc hidden-layer value using sigmoid function
        #  This is necessary and one cannot use _hl1 read from out.NN1.hl1
        #  because wgt1 value in hl1 could change... 
        tmp= 0.0
        hl1s= np.zeros(_nhl1+1)
        hl1s[0]= 1.0
        for ihl1 in range(1,_nhl1+1):
            for isf in range(_nsf+1):
                hl1s[ihl1]+= _wgt1[isf,ihl1] *gsfs[ia,isf]
            hl1s[ihl1]= sigmoid(hl1s[ihl1])
        for ihl1 in range(_nhl1+1):
            es += _wgt2[ihl1] *hl1s[ihl1]

    #.....forces
    if _fmatch:
        amls= _aml[ismpl]
        for ia in range(smpl.natm):
            for isf in range(1,_nsf+1):
                for ihl1 in range(1,_nhl1+1):
                    fs[ia,0] -= _wgt1[isf,ihl1]*_wgt2[ihl1]*amls[ia,ihl1,isf,0]
                    fs[ia,1] -= _wgt1[isf,ihl1]*_wgt2[ihl1]*amls[ia,ihl1,isf,1]
                    fs[ia,2] -= _wgt1[isf,ihl1]*_wgt2[ihl1]*amls[ia,ihl1,isf,2]
    return (es,fs)

def grad(x,*args):
    global _wgt1,_wgt2
    
    t0= time.time()
    dir= args[0]
    #.....get energies and forces
    ergs,frcs= calc_ef_from_bases(x,*args)

    p= mp.Pool(_nprcs)
    grad= np.zeros(len(x))
    
    if _nprcs == 1:
        for ismpl in range(len(_samples)):
            gs= grad_core(ismpl,ergs,frcs,x)
            for iprm in range(len(x)):
                grad[iprm] += gs[iprm]
    else:
        func_args=[]
        for ismpl in range(len(_samples)):
            func_args.append( (grad_core,ismpl,ergs,frcs,x))
        results= p.map(arg_wrapper,func_args)
        for ismpl in range(len(_samples)):
            gs= results[ismpl]
            for iprm in range(len(x)):
                grad[iprm] += gs[iprm]

    print ' ===> time NN1.grad: {0:12.3f} sec'.format(time.time()-t0)
    return grad

def grad_core(ismpl,ergs,frcs,*args):
    x      = args[0]

    wgt1,wgt2= vars2wgts(x)
    
    gs= np.zeros(len(x))
    smpl= _samples[ismpl]
    ediff= (ergs[ismpl] -_ergrefs[ismpl]) /smpl.natm
    gsfs= _gsf[ismpl]
    hl1s= _hl1[ismpl]
    iprm= 0
    for isf in range(_nsf+1):
        for ihl1 in range(1,_nhl1+1):
            tmp= 0.0
            for ia in range(smpl.natm):
                tmp += wgt2[ihl1] *hl1s[ia,ihl1] \
                         *(1.0 -hl1s[ia,ihl1]) *gsfs[ia,isf]
            gs[iprm] += 2.0*ediff*tmp 
            iprm += 1
    for ihl1 in range(_nhl1+1):
        tmp= 0.0
        for ia in range(smpl.natm):
            tmp += hl1s[ia,ihl1]
        gs[iprm] += 2.0*ediff*tmp
        iprm += 1
    
    if _fmatch:
        amls= _aml[ismpl]
        iprm= 0
        for isf in range(_nsf+1):
            for ihl1 in range(1,_nhl1+1):
                tmp= 0.0
                w2= wgt2[ihl1]
                for ia in range(smpl.natm):
                    fdiff= frcs[ismpl][ia] -_frcrefs[ismpl][ia]
                    ms= amls[ia,ihl1,isf]
                    tmp -= 2.0*( fdiff[0]*w2*ms[0] \
                                +fdiff[1]*w2*ms[1] \
                                +fdiff[2]*w2*ms[2] ) /smpl.natm/3
                gs[iprm] += tmp
                iprm += 1
        for ihl1 in range(_nhl1+1):
            tmp= 0.0
            for ia in range(smpl.natm):
                fdiff= frcs[ismpl][ia] -_frcrefs[ismpl][ia]
                for isf in range(_nsf+1):
                    ms= amls[ia,ihl1,isf]
                    tmp -= 2.0*( fdiff[0]*wgt1[isf,ihl1]*ms[0] \
                                +fdiff[1]*wgt1[isf,ihl1]*ms[1] \
                                +fdiff[2]*wgt1[isf,ihl1]*ms[2] )/smpl.natm/3
            gs[iprm] += tmp
            iprm += 1
    return gs
                    
def gather_basis(*args):
    gsf= []
    hl1= []
    aml= []
    #...read basis data
    for i in range(len(_sample_dirs)):
        dir= _sample_dirs[i]
        smpl= _samples[i]
        f1= open(_basedir+'/'+dir+'/pmd/out.NN1.gsf','r')
        f2= open(_basedir+'/'+dir+'/pmd/out.NN1.hl1','r')
        f3= open(_basedir+'/'+dir+'/pmd/out.NN1.amsl','r')
        #.....skip 1st line
        data1= f1.readline().split()
        data2= f2.readline().split()
        data3= f3.readline().split()
        gsfs= np.zeros((smpl.natm,_nsf+1))
        hl1s= np.zeros((smpl.natm,_nhl1+1))
        amls= np.zeros((smpl.natm,_nhl1+1,_nsf+1,3))
        for ia in range(smpl.natm):
            for isf in range(_nsf+1):
                data1= f1.readline().split()
                gsfs[ia,isf]= float(data1[2])
        for ia in range(smpl.natm):
            for ihl1 in range(_nhl1+1):
                data2= f2.readline().split()
                hl1s[ia,ihl1]= float(data2[2])
        for ia in range(smpl.natm):
            for ihl1 in range(_nhl1+1):
                for isf in range(_nsf+1):
                    data3= f3.readline().split()
                    amls[ia,ihl1,isf,0]= float(data3[3])
                    amls[ia,ihl1,isf,1]= float(data3[4])
                    amls[ia,ihl1,isf,2]= float(data3[5])
        gsf.append(gsfs)
        hl1.append(hl1s)
        aml.append(amls)
        f1.close()
        f2.close()
        f3.close()
    return gsf,hl1,aml

