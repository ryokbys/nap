#!/opt/local/bin/python

u"""
Make in.params.linreg to be used in linreg potential.
Since the parameters are determined randomly, 
users should re-parameterize those ones to the system to be considered.
"""

import sys,os,random

paramfname='in.params.linreg'
constfname='in.const.linreg'

pmax= 1.0
pmin=-1.0

usage='Usage: $ python {0} rcut'.format(sys.argv[0])

if len(sys.argv) != 2:
    print ' [Error] num of arguments is wrong !!!'
    print usage
    sys.exit()

if not os.path.isfile(constfname):
    print ' [Error] There must be in.const.linreg in the directory !!!'
    sys.exit()

rcin= float(sys.argv[1])

f=open(constfname,'r')
data= f.readline().split()
nelin= int(data[0])
nexp= int(data[1])
f.close()

nelem= nelin *nexp

f=open(paramfname,'w')
f.write(' {0:10d} {1:10.4f}\n'.format(nelem,rcin))
for i in range(nelem):
    f.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
    f.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
f.close()
