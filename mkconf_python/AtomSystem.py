import numpy as np

class AtomSystem(object):
    u"""
    AtomSystem has cell information and atoms.
    """

    def __init__(self):
        self.a1= np.zeros(3)
        self.a2= np.zeros(3)
        self.a3= np.zeros(3)
        self.atoms= []

    def set_lattice(self,alc,a1,a2,a3):
        self.alc= alc
        self.a1= a1
        self.a2= a2
        self.a3= a3

    def add_atom(self,atom):
        self.atoms.append(atom)

    def reset_ids(self):
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)

    def read_pmd(self,fname='pmd0000'):
        f=open(fname,'r')
        # 1st: lattice constant
        self.alc= float(f.readline().split()[0])
        # 2nd-4th: cell vectors
        self.a1= np.array([float(x)*self.alc for x in f.readline().split()])
        self.a2= np.array([float(x)*self.alc for x in f.readline().split()])
        self.a3= np.array([float(x)*self.alc for x in f.readline().split()])
        # 5th-7th: velocity of cell vectors
        tmp= f.readline().split()
        tmp= f.readline().split()
        tmp= f.readline().split()
        # 8st: num of atoms
        natm= int(f.readline().split()[0])
        # 9th-: atom positions
        self.atoms= []
        for i in range(natm):
            data= [float(x) for x in f.readline().split()]
            ai= Atom()
            ai.decode_tag(data[0])
            ai.set_pos(data[1],data[2],data[3])
            self.atoms.append(ai)
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
        f.write(" {0:10d}\n".format(len(self.atoms)))
        # atom positions
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)
            f.write(" {:22.14e} {:11.7f} {:11.7f} {:11.7f}".format(ai.tag(), \
                                                            ai.pos[0],\
                                                            ai.pos[1],\
                                                            ai.pos[2])
                    +"  {:.1f}  {:.1f}  {:.1f}".format(0.0, 0.0, 0.0)
                    +"  {:.1f}  {:.1f}".format(0.0, 0.0)
                    +"  {:.1f}  {:.1f}  {:.1f}".format(0.0, 0.0, 0.0)
                    +"  {:.1f}  {:.1f}  {:.1f}".format(0.0, 0.0, 0.0)
                    +"  {:.1f}  {:.1f}  {:.1f}".format(0.0, 0.0, 0.0)
                    +"\n")
        f.close()

    def write_akr(self,fname='akr0000'):
        f=open(fname,'w')
        # lattice constant
        f.write(" {0:12.4f}\n".format(self.alc))
        # cell vectors
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a1[0],\
                                                          self.a1[1],\
                                                          self.a1[2]))
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a2[0],\
                                                          self.a2[1],\
                                                          self.a2[2]))
        f.write(" {0:12.4f} {1:12.4f} {2:12.4f}\n".format(self.a3[0],\
                                                          self.a3[1],\
                                                          self.a3[2]))
        # num of atoms
        f.write(" {0:10d} {1:4d} {2:4d} {3:4d}\n".format(len(self.atoms),3,0,0))
        # atom positions
        for i in range(len(self.atoms)):
            ai= self.atoms[i]
            ai.set_id(i+1)
            f.write(" {:4d} {:10.5f} {:10.5f} {:10.5f}".format(ai.sid, \
                                                            ai.pos[0],\
                                                            ai.pos[1],\
                                                            ai.pos[2])
                    +"  {:.1f}  {:.1f}  {:.1f}".format(0.0, 0.0, 0.0)
                    +"\n")
        f.close()

