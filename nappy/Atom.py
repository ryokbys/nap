import numpy as np

class Atom(object):
    u"""
    Atom class defines attributions that an atom has to have.
    """

    def __init__(self):
        self.pos= np.zeros((3,))
        self.id= 0
        self.ifmv= 1
        self.sid= 1

    def set_pos(self,x,y,z):
        self.pos= np.array([x,y,z])

    def set_sid(self,sid):
        self.sid= sid

    def set_id(self,id):
        self.id= id

    def tag(self):
        tag= self.sid +self.ifmv*0.1 +self.id*1e-14
        return tag

    def decode_tag(self,tag):
        self.sid= int(tag)
        self.ifmv= int((tag-self.sid)*10)
        self.id= int(((tag-self.sid)*10 -self.ifmv)*1e+14)

