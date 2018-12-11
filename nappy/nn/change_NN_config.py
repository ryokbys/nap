#!/bin/env python
"""
Change NN configuration, such as the number of layers and
the number of nodes in a layer.
"""
from __future__ import print_function

import sys,optparse,random
import NN_io

if __name__ == '__main__':

    usage='%prog [options] in.const.NN in.params.NN'
    parser= optparse.OptionParser(usage=usage)
    (options,args)= parser.parse_args()

    if len(args) != 2:
        print(usage)
        sys.exit()

    incnst= args[0]
    inprms= args[1]

    nl,nsp,nhl,itypes,combs,consts= NN_io.read_const(incnst)
    nprm,rcut,rcut3,prms= NN_io.read_params(inprms)
    
    #.....print current configuration
    print(" Current NN configuration:")
    print("   Num of layers = ",nl)
    print("   Num of nodes in layers = ",end='')
    for il in range(nl+1):
        print(" {0:d}:{1:d}".format(il,nhl[il]),)
    print("")

    #.....read new configuration
    print("")
    print(" Let's decide new NN configuration.")
    print("")
    print(" Please enter new number of hidden layers:")
    nlnew = int(raw_input(' >>> '))
    print(" Number of hidden layers = ",nlnew)
    nhlnew= []
    nhlnew.append(nhl[0])
    for il in range(1,nlnew+1):
        print("")
        print(" Please enter the number of nodes in layer-{0:d}".format(il))
        ntmp= raw_input(' >>> ')
        nhlnew.append(int(ntmp))
        print(" Number of nodes in layer-{0:d} =".format(il),nhlnew[il])
    print("")
    print(" New NN configuration:")
    print("   Num of layers = ",nlnew)
    print("   Num of nodes in layers = ",end='')
    for il in range(nlnew+1):
        print(" {0:d}:{1:d}".format(il,nhlnew[il]),end='')
    print("")
    
    #....new parameters
    npnew= 0
    for il in range(nlnew):
        npnew += nhlnew[il]*nhlnew[il+1]
    npnew += nhlnew[nlnew]
    print(" New number of parameters =",npnew)
    prmsnew= []
    for iprm in range(npnew):
        prmsnew.append([random.random(),-1.0,1.0])

    outcnst= incnst+'.new'
    outprms= inprms+'.new'
    NN_io.write_const(outcnst,nlnew,nsp,nhlnew,itypes,combs,consts)
    NN_io.write_params(outprms,npnew,rcut,rcut3,prmsnew)

    print("")
    print(" Output files are:")
    print("   * {0:s}".format(outcnst))
    print("   * {0:s}".format(outprms))
    print(" Check these files.  Ciao ;) ")
    
    
