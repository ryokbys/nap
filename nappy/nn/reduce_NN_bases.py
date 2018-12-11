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
from __future__ import print_function

import optparse
import NN_io

incnst="in.const.NN"
inprms="in.params.NN"
innnanal="out.NN_analysis"

outcnst= 'in.const.NN.new'
outprms= 'in.params.NN.new'

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

    nl,nsp,nhl,itypes,combs,consts= NN_io.read_const(incnst)
    nprm,rcut,rcut3,prms= NN_io.read_params(inprms)
    nnanal= NN_io.read_NN_analysis(innnanal)
    
    if nhl[0] != len(nnanal):
        print('[Error] nhl[0] != len(nnanal)')
        print('  nhl[0],len(nnanal)=',nhl[0],len(nnanal))
        exit()

    focnst= open(outcnst,'w')
    foprms= open(outprms,'w')
    #...count number of bases to be remained
    nsfnew= 0
    for ihl0 in range(nhl[0]):
        if nnanal[ihl0][1] > threshold:
            nsfnew += 1
    print('num of bases to survive=',nsfnew)
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
    foprms.write(' {0:6d} {1:8.4f} {2:7.3f}\n'.format(npnew,rcut,rcut3))
    
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
            for i,p in enumerate(prms[offset+ip]):
                fp = float(p)
                if i == 0:
                    foprms.write(' {0:20.10e}'.format(fp))
                else:
                    foprms.write(' {0:9.5f}'.format(fp))
            foprms.write('\n')
            
    offset= nsfold*nhl[1]
    for prm in prms[offset:]:
        for i,p in enumerate(prm):
            fp = float(p)
            if i == 0:
                foprms.write(' {0:20.10e}'.format(fp))
            else:
                foprms.write(' {0:9.5f}'.format(fp))
        foprms.write('\n')
    
    focnst.close()
    foprms.close()
