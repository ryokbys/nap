#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
This class provides read and write functions of VASP POSCAR file.
"""

import numpy as np
import copy

class POSCAR(object):
    """
    POSCAR class enables POSCAR read/write functions.
    """
    
    def __init__(self,fname="POSCAR"):
        self.fname = fname
        self.h = np.zeros((3,3),dtype=float)
        self.num_atoms= []
        self.pos= []
        self.flags= []
        self.species= []

#     def to_dict(self):
#         """
#         Returns a dictionary-type variable of the object.
#         """
#         dict= {}
#         dict.update({'comment':self.c1})
#         dict.update({'afac':self.afac})
#         dict.update({})

    def read(self,fname = 'POSCAR'):
        self.fname = fname
        f= open(fname,'r')
        #.....1st line: comment
        self.c1= f.readline()
        #.....2nd line: multiplying factor
        self.afac= float(f.readline().split()[0])
        #.....3rd-5th lines: lattice vectors
        data= f.readline().split()
        self.h[0]= [ float(x) for x in data ]
        data= f.readline().split()
        self.h[1]= [ float(x) for x in data ]
        data= f.readline().split()
        self.h[2]= [ float(x) for x in data ]
        #.....6th line: num of atoms
        data= f.readline().split()
        if not data[0].isdigit(): # if it is not digit, read next line
            self.species = copy.copy(data)
            data = f.readline().split()
        self.num_atoms= np.array([ int(n) for n in data ])
        #.....7th line: comment
        self.c7= f.readline()
        if self.c7[0] in ('s','S'):
            self.c8= f.readline()
        #.....following lines: atom positions
        for ni in self.num_atoms:
            for j in range(ni):
                data= f.readline().split()
                self.pos.append(np.array([ float(x) for x in data[0:3] ]))
                if len(data) > 3:
                    if len(data) == 6:
                        self.flags.append([ x for x in data[3:6] ])
                    elif len(data) == 4:
                        self.flags.append([ data[3],data[3],data[3] ])
                else:
                    self.flags.append([ 'T', 'T', 'T' ])
        f.close()

    def write(self,fname='POSCAR'):
        f= open(fname,'w')
        f.write(self.c1)
        f.write(' {0:12.7f}\n'.format(self.afac))
        f.write(' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(self.h[0,0],
                                                          self.h[0,1],
                                                          self.h[0,2]))
        f.write(' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(self.h[1,0],
                                                          self.h[1,1],
                                                          self.h[1,2]))
        f.write(' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(self.h[2,0],
                                                          self.h[2,1],
                                                          self.h[2,2]))
        for n in self.num_atoms:
            f.write(' {0:3d}'.format(n))
        f.write('\n')
        
        f.write(self.c7)
        if hasattr(self, 'c8'):
            f.write(self.c8)

        for i in range(len(self.pos)):
            f.write(' {0:12.7f} {1:12.7f} {2:12.7f}'.format(self.pos[i][0],
                                                            self.pos[i][1],
                                                            self.pos[i][2]))
            f.write(' {0} {1} {2}\n'.format(self.flags[i][0],
                                            self.flags[i][1],
                                            self.flags[i][2]))
        f.close()


if __name__ == '__main__':
    poscar= POSCAR()
    poscar.read()
    poscar.write(fname='TEST-POSCAR')

