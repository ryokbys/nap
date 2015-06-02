#!/bin/env python

import sys,os
from AtomSystem import AtomSystem

infname= sys.argv[1]
print 'input file: ',infname
aSys= AtomSystem()
aSys.read_akr(infname)
natm= len(aSys.atoms)

outfname= infname+".voro"
f= open(outfname,'w')
alc= aSys.alc
ax= aSys.a1[0]
ay= aSys.a2[1]
az= aSys.a3[2]
for ia in range(natm):
    pos= aSys.atoms[ia].pos
    f.write(' {0:4d}'.format(ia+1))
    f.write(' {0:10.3f}'.format(alc*ax*pos[0]))
    f.write(' {0:10.3f}'.format(alc*ay*pos[1]))
    f.write(' {0:10.3f}'.format(alc*az*pos[2]))
    f.write('\n')
f.close()

print 'Use voro++ as follows:'
cmd= 'voro++ -p -c " %i %q %v %s %A"' \
    +' 0.0 {0:5.1f}'.format(alc*ax) \
    +' 0.0 {0:5.1f}'.format(alc*ay) \
    +' 0.0 {0:5.1f}'.format(alc*az) \
    +' {0:s}'.format(outfname)

print '$ '+cmd
os.system(cmd)

