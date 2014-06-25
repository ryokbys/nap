
import sys
import numpy as np

def read_old_pmd(fname):
    f= open(fname,'r')
    natm= int(f.readline().split()[0])
    a1= [float(x) for x in f.readline().split()]
    a2= [float(x) for x in f.readline().split()]
    a3= [float(x) for x in f.readline().split()]
    tmp= f.readline().split()
    tmp= f.readline().split()
    tmp= f.readline().split()
    tag= np.zeros(natm)
    pos= np.zeros((natm,3))
    for i in range(natm):
        data= [float(x) for x in f.readline().split()]
        tag[i]= data[0]
        pos[i]= data[1:4]
    f.close()
    return natm,a1,a2,a3,tag,pos

def write_new_pmd(fname,natm,alc,a1,a2,a3,tag,pos):
    f= open(fname,'w')
    f.write(' {0:15.7f}\n'.format(alc))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(a1[0],a1[1],a1[2]))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(a2[0],a2[1],a2[2]))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(a3[0],a3[1],a3[2]))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(0.0, 0.0, 0.0))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(0.0, 0.0, 0.0))
    f.write(' {0:15.7f} {1:15.7f} {2:15.7f}\n'.format(0.0, 0.0, 0.0))
    f.write(' {0:10d}\n'.format(natm))
    for i in range(natm):
        f.write(' {0:22.14e} {1:22.14e} {2:22.14e} {3:22.14e}'.format(tag[i],pos[i,0],pos[i,1],pos[i,2]))
        for j in range(12):
            f.write(' 0.0')
        f.write('\n')
    f.close()

if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print ' [Error] num of args is wrong !!!'
        print '  Usage: $ {0} <old-pmd> <new-pmd>'.format(sys.argv[0])
        sys.exit()

    oldpmd= sys.argv[1]
    newpmd= sys.argv[2]

    natm,a1,a2,a3,tag,pos= read_old_pmd(oldpmd)
    write_new_pmd(newpmd,natm,1.0,a1,a2,a3,tag,pos)
