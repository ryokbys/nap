#!/bin/env python
"""
Smear 1D data of equidistance intervals with Gaussians.
"""

import os,sys,math
import optparse

usage= '%prog [options] datafile'

parser= optparse.OptionParser(usage=usage)
parser.add_option("-s","--sigma",dest="sigma",type="float",
                  default=2.0,
                  help="width of Gaussian (sigma) in the unit of dx.")
parser.add_option("-x",dest="xp",type="int",
                  default=1,
                  help="position of x data (start from 1).")
parser.add_option("-y",dest="yp",type="int",
                  default=2,
                  help="position of y data (start from 1).")
(options,args)= parser.parse_args()

xp= options.xp -1
yp= options.yp -1
sigma= options.sigma
infname= args[0]

infile= open(infname,'r')
nline= 0
for line in infile.readlines():
    if line[0] == "#":
        continue
    nline += 1
print ' num of data points in the file=',nline
infile.seek(0)
data= [ [0.0,0.0] for i in range(nline) ]
il= 0
for line in infile.readlines():
    if line[0] == "#":
        continue
    sline= line.split()
    data[il][0]= float(sline[xp])
    data[il][1]= float(sline[yp])
    il += 1
infile.close()

dx= data[1][0]-data[0][0]
sgm= sigma*dx
isgm= int(sigma)+1
pref= 1.0/(math.sqrt(2.0*math.pi)*sgm)

gdat= [ [0.0,0.0] for i in range(nline) ]

for ix in range(nline):
    gdat[ix][0]= data[ix][0]
    for jx in range(-nline+1,nline-1):
        kx= ix+jx
        if kx < 0:
            kx = -kx
        elif kx >= nline:
            kx = nline -(kx-(nline-1))
        gdat[ix][1] += data[kx][1]*pref*math.exp(-(dx*(jx))**2/sgm**2)*dx

outfile= open(infname+'.smeared','w')
for il in range(nline):
    outfile.write(' {0:15.2f} {1:15.7f}\n'.format(gdat[il][0],gdat[il][1]))
outfile.close()

print ' write '+infname+'.smeared'
print ' GAUSSIAN_SMEAR done '

# dtotal= 0.0
# gtotal= 0.0
# for ix in range(nline):
#     if ix == 0 or ix == nline-1:
#         dtotal += data[ix][1]*dx/2
#         gtotal += data[ix][1]*dx/2
#     else:
#         dtotal += data[ix][1]*dx
#         gtotal += data[ix][1]*dx
# print ' dtotal=',dtotal
# print ' gtotal=',gtotal
