#!/usr/bin/python
# Time-stamp: <2015-01-28 14:21:32 Ryo KOBAYASHI>
#
# Reduce pmd/{erg,frc}0000 by using number of atoms in ./pos file.
#

inerg='pmd/erg0000'
infrc='pmd/frc0000'
inpos='./pos'

outerg='erg.pmd'
outfrc='frc.pmd'

fpos=open(inpos,'r')
alc= fpos.readline().split()[0]
for l in range(6):
    tmp= fpos.readline()
natm=int(fpos.readline().split()[0])
fpos.close()

#.....erg
ferg=open(inerg,'r')
erg= float(ferg.readline().split()[0])
ferg.close()

#.....frc
ffrc=open(infrc,'r')
natmf=int(ffrc.readline().split()[0])

fouterg=open(outerg,'w')
fouterg.write(" {0:22.14e}\n".format(erg/natmf*natm))
fouterg.close()

foutfrc=open(outfrc,'w')

frc= []
foutfrc.write(" {0:10d}\n".format(natm))
for i in range(natm):
    foutfrc.write(ffrc.readline())
ffrc.close()
foutfrc.close()
