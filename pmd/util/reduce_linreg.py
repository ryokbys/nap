#!/opt/local/bin/python
u"""
Reduce params and constants of linreg potential, and write
new in.params.linreg and in.const.linreg.
"""

import sys,os,random

paramfname='in.params.linreg'
constfname='in.const.linreg'

#.....if the param is smaller than the maximum-value*criterion,
#     the param will be removed.
criterion= 0.01

usage='Usage: $ python {0}'.format(sys.argv[0])

if not len(sys.argv)==1 :
    print ' [Error] num of arguments is wrong !!!'
    print usage
    sys.exit()

if not os.path.isfile(constfname):
    print ' [Error] There must be in.const.linreg in the directory !!!'
    sys.exit()

if not os.path.isfile(paramfname):
    print ' [Error] There must be in.params.linreg in the directory !!!'
    sys.exit()

#.....backup old in-files
os.system('cp in.const.linreg in.const.linreg.old')
os.system('cp in.params.linreg in.params.linreg.old')

f=open(constfname,'r')
data= f.readline().split()
nelin0= int(data[0])
nexp0= int(data[1])
consts= []
for line in f.readlines():
    consts.append(line.split())
f.close()

f=open(paramfname,'r')
data= f.readline().split()
nelem0= int(data[0])
rcin0= float(data[1])
params= []
for line in f.readlines():
    params.append(line.split())
f.close()

#.....search maximum value
pmax= 0.0
for param in params:
    p= abs(float(param[0]))
    pmax= max(pmax,p)

#.....reduce params whose abs values are less than criterion*pmax
to_be_reduced= [ False for i in range(len(params))]
nelem= 0
for i in range(len(params)):
    p= abs(float(params[i][0]))
    if p < pmax*criterion:
        to_be_reduced[i] = True
    else:
        nelem += 1

#.....write reduce const and params
g1= open(constfname,'w')
g1.write(' {0:10d} {1:10d}\n'.format(nelem,nexp0))
g2= open(paramfname,'w')
g2.write(' {0:10d} {1:10.4f}\n'.format(nelem,rcin0))
for i in range(len(params)):
    if to_be_reduced[i]:
        continue
    #...in.consts.linreg
    g1.write(' {0} {1}'.format(consts[i][0],consts[i][1]))
    if int(consts[i][0]) == 1: # Gaussian
        g1.write(' {0} {1}\n'.format(consts[i][2],
                                     consts[i][3],
                                     consts[i][4]))
    elif int(consts[i][0]) == 2: # cosine
        g1.write(' {0} {1}\n'.format(consts[i][2],
                                     consts[i][3]))
    elif int(consts[i][0]) == 3: # angular
        g1.write(' {0}\n'.format(consts[i][2]))
    elif int(consts[i][0]) >= 4: # polynomial
        g1.write(' {0}\n'.format(consts[i][2]))
    #...in.params.linreg
    g2.write(' {0} {1} {2}\n'.format(params[i][0],
                                     params[i][1],
                                     params[i][2]))
g1.close()
g2.close()

