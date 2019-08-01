#!/usr/bin/env python
"""
Define *Atom* class that provides atom information such as id, species (sid),
positions, velocities, and some flags used in `pmd`.

Usage:
    atom.py [options] SYMBOL [SYMBOL...]

Options:
    -h, --help  Show this help message and exit.

"""
from __future__ import print_function

import os,sys
import numpy as np
from docopt import docopt

# sys.path.append(os.path.dirname(__file__))
import elements


def get_symbol_from_number(number):
    if number > len(elements.elements):
        raise ValueError('No element of number {0:d}'.format(number))
    for k,v in elements.elements.items():
        if v['number'] == number:
            return k
    raise ValueError('No element found for number {0:d}'.format(number))


def get_number_from_symbol(symbol):
    try:
        element = elements.elements[symbol]
    except:
        raise ValueError('No symbol {0:s}'.format(symbol))
    return element['number']



class Atom(object):
    """
    Define attributions of an atom.
    """

    def __init__(self,pos=[0.,0.,0.],sid=1,symbol=None):
        self.pos= np.array(pos)
        self.vel= np.zeros((3,))
        self.frc= np.zeros((3,))
        self.strs= np.zeros((6,))
        self.epot = 0.0
        self.ekin = 0.0
        self.id= 0
        self.ifmv= 1
        self.sid= sid
        if symbol:
            self.symbol = symbol
        else:
            self.symbol = get_symbol_from_number(self.sid)
        self.aux= {}

    def set_pos(self,x,y,z):
        self.pos= np.array([x,y,z],dtype=float)

    def set_vel(self,x,y,z):
        self.vel= np.array([x,y,z],dtype=float)


    def set_frc(self,x,y,z):
        """
        Set force on this atom.
        The force should be not scaled and eV/A unit in Cartessian coordinate.
        """
        self.frc= np.array([x,y,z],dtype=float)


    def set_strs(self,xx,yy,zz,yz,xz,xy):
        self.strs= np.array([xx,yy,zz,yz,xz,xy],dtype=float)

    def get_pressure(self):
        p = 0.0
        p = (self.strs[0] +self.strs[1] +self.strs[2])/3
        return p
        
    def set_epot(self,epot):
        self.epot = epot

    def set_ekin(self,ekin):
        self.ekin = ekin

    def set_sid(self,sid):
        self.sid= int(sid)

    def set_symbol(self,symbol):
        self.symbol = symbol

    def set_id(self,id):
        self.id= id

    def set_ifmv(self,ifmv):
        self.ifmv = ifmv
    
    def set_aux(self,key,value):
        self.aux[key] = value
        
    def tag(self):
        tag= self.sid +self.ifmv*0.1 +self.id*1e-14
        return tag

    def decode_tag(self,tag):
        self.sid= int(tag)
        self.ifmv= int((tag-self.sid)*10)
        self.id= int(((tag-self.sid)*10 -self.ifmv)*1e+14)

def write_info_of_element(*symbols):
    for s in symbols:
        try:
            element = elements.elements[s]
        except:
            raise ValueError('No symbol {0:s}'.format(s))
        print('Element: {0:s}'.format(s))
        print('  name: {0:s}'.format(element['name']))
        print('  atomic number: {0:d}'.format(element['number']))
        print('  mass: {0:8.3f}'.format(element['mass']))
        print('  abundance: {0:s}'.format(element['abundance']))
        print('')
    return None

if __name__ == '__main__':
    
    args= docopt(__doc__)
    symbols = args['SYMBOL']
    write_info_of_element(*symbols)


