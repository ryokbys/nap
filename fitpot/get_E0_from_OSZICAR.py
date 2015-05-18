#!/opt/local/bin/python

fo= open('erg.ref','w')

f=open('OSZICAR','r')
e0= 1.0e+30
for line in f.readlines():
    dat= line.split()
    if 'E0=' in dat:
        idx= dat.index('E0=')
        #print dat[idx+1]
        e0= float(dat[idx+1])
fo.write(' {0:20.7}\n'.format(e0))
f.close()
fo.close()
