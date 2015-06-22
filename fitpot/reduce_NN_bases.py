#!/bin/env python
"""
Reduces the number of NN bases using the information from out.NN_analysis.

INPUT:
  * in.const.NN
  * in.params.NN
  * out.NN_analysis

OUTPUT:
  * in.const.NN.new
  * in.params.NN.new
"""

import optparse

incnst="in.const.NN"
inprms="in.params.NN"
innnanal="out.NN_analysis"

outcnst= 'in.const.NN.new'
outprms= 'in.params.NN.new'

def read_const(fname):
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

def read_params(fname):
    f=open(fname,'r')
    buf= f.readline().split()
    nprm= int(buf[0])
    rcut= float(buf[1])
    prms= []
    for iprm in range(nprm):
        buf= f.readline().split()
        prms.append(buf[0:3])
    f.close()
    print ' reading {0:s} done.'.format(fname)
    return nprm,rcut,prms

def read_NN_analysis(fname):
    f= open(fname,'r')
    nnanal= []
    for line in f.readlines():
        buf= line.split()
        nnanal.append((int(buf[1]),float(buf[2])))
    f.close()
    print ' reading {0:s} done.'.format(fname)
    return nnanal

if __name__ == '__main__':

    usage= '%prog [options] in.const.NN in.params.NN out.NN_analysis'

    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-t","--threshold",dest="threshold",type="float",
                      default=1e-8,
                      help="Threshold value less than which are to be elminiated.")
    (options,args)= parser.parse_args()

    threshold= options.threshold
    incnst= args[0]
    inprms= args[1]
    innnanal= args[2]

    nl,nsp,nhl,itypes,combs,consts= read_const(incnst)
    nprm,rcut,prms= read_params(inprms)
    nnanal= read_NN_analysis(innnanal)
    
    if nhl[0] != len(nnanal):
        print '[Error] nhl[0] != len(nnanal)'
        print '  nhl[0],len(nnanal)=',nhl[0],len(nnanal)
        exit()

    focnst= open(outcnst,'w')
    foprms= open(outprms,'w')
    #...count number of bases to be remained
    nsfnew= 0
    for ihl0 in range(nhl[0]):
        if nnanal[ihl0][1] > threshold:
            nsfnew += 1
    print 'num of bases to survive=',nsfnew
    nsfold= nhl[0]
    nhl[0]= nsfnew
    npnew= 0
    for il in range(nl):
        npnew += nhl[il]*nhl[il+1]
    npnew += nhl[nl]

    #....write headers
    focnst.write(' {0:4d} {1:4d}'.format(nl,nsp))
    for il in range(nl+1):
        focnst.write(' {0:5}'.format(nhl[il]))
    focnst.write('\n')
    foprms.write(' {0:6d} {1:8.4f}\n'.format(npnew,rcut))
    
    for isf in range(nsfold):
        if nnanal[isf][1] < threshold:
            continue
        #....const
        focnst.write(' {0:3d}  '.format(itypes[isf]))
        for ic in range(len(combs[isf])):
            focnst.write(' {0:2d}'.format(combs[isf][ic]))
        focnst.write('   ')
        for ic in range(len(consts[isf])):
            focnst.write(' {0:s}'.format(consts[isf][ic]))
        focnst.write('\n')
        #....params
        offset= isf*nhl[1]
        for ip in range(nhl[1]):
            for p in prms[offset+ip]:
                foprms.write('  {0:s}'.format(p))
            foprms.write('\n')
            
    offset= nsfold*nhl[1]
    for prm in prms[offset:]:
        for p in prm:
            foprms.write('  {0:s}'.format(p))
        foprms.write('\n')
    
    focnst.close()
    foprms.close()
