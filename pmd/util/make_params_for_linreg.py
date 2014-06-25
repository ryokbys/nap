#!/opt/local/bin/python

u"""
Make in.params.linreg to be used in linreg potential.
Since the parameters are determined randomly, 
users should re-parameterize those ones to the system to be considered.
"""

import sys,os,random

paramfname='in.params.linreg'

pmax= 1.0
pmin=-1.0

usage='Usage: $ python {0} <nelem> <rc>'

if len(sys.argv) != 3:
    print ' [Error] num of arguments is wrong !!!'
    print usage
    sys.exit()

nelem= int(sys.argv[1])
rcin= float(sys.argv[2])

f=open(paramfname,'w')
f.write(' {0:10d} {1:10.4f}\n'.format(nelem,rcin))
for i in range(nelem):
    f.write(' {0:10.6f}'.format(random.uniform(pmin,pmax)))
    f.write(' {0:10.4f} {1:10.4f}\n'.format(pmin,pmax))
f.close()
