#!/bin/env python
"""
NN_io includes some functions related to input/output of NN potential.
"""


def read_const(fname):
    """
    Read in.const.NN file.
    """
    f= open(fname,'r')
    buf= f.readline().split()
    nl= int(buf[0])
    nsp=int(buf[1])
    nhl= []
    for il in range(nl+1):
        nhl.append(int(buf[2+il]))
    combs= []
    consts= []
    itypes= []
    for ihl0 in range(nhl[0]):
        buf= f.readline().split()
        itype= int(buf[0])
        itypes.append(itype)
        if itype <= 100:  # 2-body
            ia= int(buf[1])
            ja= int(buf[2])
            combs.append((ia,ja))
            consts.append(buf[3:])
        else:    # 3-body
            ia= int(buf[1])
            ja= int(buf[2])
            ka= int(buf[3])
            combs.append((ia,ja,ka))
            consts.append(buf[4:])
    f.close()
    print ' reading {0:s} done.'.format(fname)
    return nl,nsp,nhl,itypes,combs,consts

def write_const(fname,nl,nsp,nhl,itypes,combs,consts):
    """
    Write in.const.NN file.
    """
    focnst= open(fname,'w')
    focnst.write(' {0:4d} {1:4d}'.format(nl,nsp))
    for il in range(nl+1):
        focnst.write(' {0:5}'.format(nhl[il]))
    focnst.write('\n')
    for isf in range(nhl[0]):
        #....const
        focnst.write(' {0:3d}  '.format(itypes[isf]))
        for ic in range(len(combs[isf])):
            focnst.write(' {0:2d}'.format(combs[isf][ic]))
        focnst.write('   ')
        for ic in range(len(consts[isf])):
            focnst.write(' {0:s}'.format(consts[isf][ic]))
        focnst.write('\n')
    focnst.close()

def read_params(fname):
    f=open(fname,'r')
    buf= f.readline().split()
    nprm= int(buf[0])
    rcut= float(buf[1])
    rcut3 = float(buf[2])
    prms= []
    for iprm in range(nprm):
        buf= f.readline().split()
        prms.append(buf[0:3])
    f.close()
    print ' reading {0:s} done.'.format(fname)
    return nprm,rcut,rcut3,prms

def write_params(fname,nprm,rcut,rcut3,prms):
    foprms= open(fname,'w')
    foprms.write(' {0:6d} {1:8.4f} {2:7.3f}\n'.format(nprm,rcut,rcut3))
    for p in prms:
        foprms.write(' {0:15.7e}  {1:8.4f}  {2:8.4f}\n'.format(p[0],p[1],p[2]))
    foprms.close()

def read_NN_analysis(fname):
    f= open(fname,'r')
    nnanal= []
    for line in f.readlines():
        buf= line.split()
        nnanal.append((int(buf[1]),float(buf[2])))
    f.close()
    print ' reading {0:s} done.'.format(fname)
    return nnanal

