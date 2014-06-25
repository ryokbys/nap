import numpy as np

class MD_system:
    u"""inclues information: number of atoms, cell vectors,
     atom positions, total energy of the system, and forces on atoms.
    """

    def __init__(self):
        self.natm= 0
        self.a1= np.zeros(3)
        self.a2= np.zeros(3)
        self.a3= np.zeros(3)
        self.initialize_arrays()
        
    def initialize_arrays(self):
        self.tag= np.zeros(self.natm)
        self.pos= np.zeros([self.natm,3])
        self.erg= 0.0
        self.frc= np.zeros([self.natm,3])

    def set_id(self,id):
        self.id= id

    def set_natm(self,natm):
        self.natm= natm
        self.initialize_arrays()

    def num_of_species(self,species):
        num= 0
        for i in range(self.natm):
            if int(self.tag[i]) == species:
                num += 1
        return num

    def read_pmd(self,fname='pmd0000'):
        f=open(fname,'r')
        # 1st: lattice constant
        self.alc= float(f.readline().split()[0])
        # 2nd-4th: cell vectors
        self.a1= [float(x)*self.alc for x in f.readline().split()]
        self.a2= [float(x)*self.alc for x in f.readline().split()]
        self.a3= [float(x)*self.alc for x in f.readline().split()]
        # 5th-7th: velocity of cell vectors
        tmp= f.readline().split()
        tmp= f.readline().split()
        tmp= f.readline().split()
        # 8st: num of atoms
        self.natm= int(f.readline().split()[0])
        self.initialize_arrays()
        # 9th-: atom positions
        for i in range(self.natm):
            data= [float(x) for x in f.readline().split()]
            self.tag[i]= data[0]
            self.pos[i]= data[1:4] # note: data[1:4]=[data[1],data[2],data[3]]
        f.close()

    def write_pmd(self,fname='pmd0000'):
        f=open(fname,'w')
        # lattice constant
        f.write(" {0:15.7f}\n".format(self.alc))
        # cell vectors
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a1[0],\
                                                          self.a1[1],\
                                                          self.a1[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a2[0],\
                                                          self.a2[1],\
                                                          self.a2[2]))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(self.a3[0],\
                                                          self.a3[1],\
                                                          self.a3[2]))
        # velocities of cell vectors
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(0.0, 0.0, 0.0))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(0.0, 0.0, 0.0))
        f.write(" {0:15.7f} {1:15.7f} {2:15.7f}\n".format(0.0, 0.0, 0.0))
        # num of atoms
        f.write(" {0:10d}\n".format(self.natm))
        # atom positions
        for i in range(self.natm):
            f.write(" {} {:11.7f} {:11.7f} {:11.7f}".format(self.tag[i], \
                                                            self.pos[i,0],\
                                                            self.pos[i,1],\
                                                            self.pos[i,2])
                    +" {:.1f} {:.1f} {:.1f}".format(0.0, 0.0, 0.0)
                    +" {:.1f} {:.1f}".format(0.0, 0.0)
                    +" {:.1f} {:.1f} {:.1f}".format(0.0, 0.0, 0.0)
                    +" {:.1f} {:.1f} {:.1f}".format(0.0, 0.0, 0.0)
                    +" {:.1f} {:.1f} {:.1f}".format(0.0, 0.0, 0.0)
                    +"\n")
        f.close()
