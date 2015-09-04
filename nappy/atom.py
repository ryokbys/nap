#!/usr/bin/env python
"""
Define *atom* class that provides atom information such as id, species (sid),
positions, velocities, and some flags used in `pmd`.
"""

import numpy as np

class atom(object):
    """
    Define attributions of an atom.
    """

    def __init__(self):
        self.pos= np.zeros((3,))
        self.vel= np.zeros((3,))
        self.id= 0
        self.ifmv= 1
        self.sid= 1

    def set_pos(self,x,y,z):
        self.pos= np.array([x,y,z])

    def set_vel(self,x,y,z):
        self.vel= np.array([x,y,z])

    def set_sid(self,sid):
        self.sid= int(sid)

    def set_id(self,id):
        self.id= id

    def tag(self):
        tag= self.sid +self.ifmv*0.1 +self.id*1e-14
        return tag

    def decode_tag(self,tag):
        self.sid= int(tag)
        self.ifmv= int((tag-self.sid)*10)
        self.id= int(((tag-self.sid)*10 -self.ifmv)*1e+14)

