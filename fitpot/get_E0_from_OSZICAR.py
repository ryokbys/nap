#!/opt/local/bin/python

fo= open('erg.ref','w')

f=open('OSZICAR','r')
for line in f.readlines():
    dat= line.split()
    if 'E0=' in dat:
        idx= dat.index('E0=')
        #print dat[idx+1]
        fo.write(' {0:20.7}\n'.format(float(dat[idx+1])))
        exit()
f.close()
fo.close()
