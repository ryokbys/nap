#!/usr/bin/env python
"""
Define *Atom* class that provides atom information such as id, species (sid),
positions, velocities, and some flags used in `pmd`.

Usage:
    atom.py [options] SYMBOL [SYMBOL...]

Options:
    -h, --help  Show this help message and exit.

"""

import numpy as np
from docopt import docopt

ELEMENTS = {
    'H'  :{'number':  1, 'name':'Hydrogen',      'weight':1.008 },
    'He' :{'number':  2, 'name':'Hellium',       'weight':4.0026},
    'Li' :{'number':  3, 'name':'Lithium',       'weight':6.94  },
    'Be' :{'number':  4, 'name':'Beryllium',     'weight':9.0122},
    'B'  :{'number':  5, 'name':'Boron',         'weight':10.81 },
    'C'  :{'number':  6, 'name':'Carbon',        'weight':12.011},
    'N'  :{'number':  7, 'name':'Nitrogen',      'weight':14.007},
    'O'  :{'number':  8, 'name':'Oxygen',        'weight':15.999},
    'F'  :{'number':  9, 'name':'Fluorine',      'weight':18.998},
    'Ne' :{'number': 10, 'name':'Neon',          'weight':20.180},
    'Na' :{'number': 11, 'name':'Sodium',        'weight':22.990},
    'Mg' :{'number': 12, 'name':'Magnesium',     'weight':24.305},
    'Al' :{'number': 13, 'name':'Aluminum',      'weight':26.982},
    'Si' :{'number': 14, 'name':'Silicon',       'weight':28.085},
    'P'  :{'number': 15, 'name':'Phosphorus',    'weight':30.974},
    'S'  :{'number': 16, 'name':'Sulfur',        'weight':32.06 },
    'Cl' :{'number': 17, 'name':'Chlorine',      'weight':35.45 },
    'Ar' :{'number': 18, 'name':'Argon',         'weight':39.948},
    'K'  :{'number': 19, 'name':'Potassium',     'weight':39.098},
    'Ca' :{'number': 20, 'name':'Calcium',       'weight':40.078},
    'Sc' :{'number': 21, 'name':'Scandium',      'weight':44.956},
    'Ti' :{'number': 22, 'name':'Titanium',      'weight':47.867},
    'V'  :{'number': 23, 'name':'Vanadium',      'weight':50.942},
    'Cr' :{'number': 24, 'name':'Chromium',      'weight':51.996},
    'Mn' :{'number': 25, 'name':'Manganese',     'weight':54.938},
    'Fe' :{'number': 26, 'name':'Iron',          'weight':55.845},
    'Co' :{'number': 27, 'name':'Cobalt',        'weight':58.933},
    'Ni' :{'number': 28, 'name':'Nickel',        'weight':58.693},
    'Cu' :{'number': 29, 'name':'Copper',        'weight':63.546},
    'Zn' :{'number': 30, 'name':'Zinc',          'weight':65.38 },
    'Ga' :{'number': 31, 'name':'Gallium',       'weight':69.723},
    'Ge' :{'number': 32, 'name':'Germanium',     'weight':72.631},
    'As' :{'number': 33, 'name':'Arsenic',       'weight':74.922},
    'Se' :{'number': 34, 'name':'Selenium',      'weight':78.972},
    'Br' :{'number': 35, 'name':'Bromine',       'weight':79.904},
    'Kr' :{'number': 36, 'name':'Krypton',       'weight':84.798},
    'Rb' :{'number': 37, 'name':'-------------', 'weight':85.468},
    'Sr' :{'number': 38, 'name':'-------------', 'weight':87.62 },
    'Y'  :{'number': 39, 'name':'-------------', 'weight':88.906},
    'Zr' :{'number': 40, 'name':'-------------', 'weight':91.224},
    'Nb' :{'number': 41, 'name':'-------------', 'weight':92.906},
    'Mo' :{'number': 42, 'name':'-------------', 'weight':95.95 },
    'Tc' :{'number': 43, 'name':'-------------', 'weight':98.907},
    'Ru' :{'number': 44, 'name':'-------------', 'weight':101.07},
    'Rh' :{'number': 45, 'name':'-------------', 'weight':102.906},
    'Pd' :{'number': 46, 'name':'-------------', 'weight':106.420},
    'Ag' :{'number': 47, 'name':'-------------', 'weight':107.868},
    'Cd' :{'number': 48, 'name':'-------------', 'weight':112.411},
    'In' :{'number': 49, 'name':'-------------', 'weight':114.818},
    'Sn' :{'number': 50, 'name':'-------------', 'weight':118.711},
    'Sb' :{'number': 51, 'name':'-------------', 'weight':121.760},
    'Te' :{'number': 52, 'name':'-------------', 'weight':127.600},
    'I'  :{'number': 53, 'name':'-------------', 'weight':126.904},
    'Xe' :{'number': 54, 'name':'-------------', 'weight':131.294},
    'Cs' :{'number': 55, 'name':'-------------', 'weight':132.905},
    'Ba' :{'number': 56, 'name':'-------------', 'weight':137.328},
    'La' :{'number': 57, 'name':'-------------', 'weight':138.905},
    'Ce' :{'number': 58, 'name':'-------------', 'weight':140.116},
    'Pr' :{'number': 59, 'name':'-------------', 'weight':140.908},
    'Nd' :{'number': 60, 'name':'-------------', 'weight':144.242},
    'Pm' :{'number': 61, 'name':'-------------', 'weight':144.913},
    'Sm' :{'number': 62, 'name':'-------------', 'weight':150.360},
    'Eu' :{'number': 63, 'name':'-------------', 'weight':151.964},
    'Gd' :{'number': 64, 'name':'-------------', 'weight':157.250},
    'Tb' :{'number': 65, 'name':'-------------', 'weight':158.925},
    'Dy' :{'number': 66, 'name':'-------------', 'weight':162.500},
    'Ho' :{'number': 67, 'name':'-------------', 'weight':164.930},
    'Er' :{'number': 68, 'name':'-------------', 'weight':167.259},
    'Tm' :{'number': 69, 'name':'-------------', 'weight':168.934},
    'Yb' :{'number': 70, 'name':'-------------', 'weight':173.055},
    'Lu' :{'number': 71, 'name':'-------------', 'weight':174.967},
    'Hf' :{'number': 72, 'name':'-------------', 'weight':178.490},
    'Ta' :{'number': 73, 'name':'-------------', 'weight':180.948},
    'W'  :{'number': 74, 'name':'-------------', 'weight':183.840},
    'Re' :{'number': 75, 'name':'-------------', 'weight':186.207},
    'Os' :{'number': 76, 'name':'-------------', 'weight':190.230},
    'Ir' :{'number': 77, 'name':'-------------', 'weight':192.217},
    'Pt' :{'number': 78, 'name':'-------------', 'weight':195.085},
    'Au' :{'number': 79, 'name':'-------------', 'weight':196.967},
    'Hg' :{'number': 80, 'name':'-------------', 'weight':200.592},
    'Tl' :{'number': 81, 'name':'-------------', 'weight':204.383},
    'Pb' :{'number': 82, 'name':'-------------', 'weight':207.200},
    'Bi' :{'number': 83, 'name':'-------------', 'weight':208.980},
    'Po' :{'number': 84, 'name':'-------------', 'weight':208.982},
    'At' :{'number': 85, 'name':'-------------', 'weight':209.987},
    'Rn' :{'number': 86, 'name':'-------------', 'weight':222.018},
    'Fr' :{'number': 87, 'name':'-------------', 'weight':223.200},
    'Ra' :{'number': 88, 'name':'-------------', 'weight':226.025},
    'Ac' :{'number': 89, 'name':'-------------', 'weight':227.028},
    'Th' :{'number': 90, 'name':'-------------', 'weight':232.038},
    'Pa' :{'number': 91, 'name':'-------------', 'weight':231.036},
    'U'  :{'number': 92, 'name':'-------------', 'weight':238.029},
    'Np' :{'number': 93, 'name':'-------------', 'weight':237.048},
    'Pu' :{'number': 94, 'name':'-------------', 'weight':244.064},
    'Am' :{'number': 95, 'name':'-------------', 'weight':243.061},
    'Cm' :{'number': 96, 'name':'-------------', 'weight':247.000},
    'Bk' :{'number': 97, 'name':'-------------', 'weight':247.070},
    'Cf' :{'number': 98, 'name':'-------------', 'weight':251.080},
    'Es' :{'number': 99, 'name':'-------------', 'weight':254.000},
    'Fm' :{'number':100, 'name':'-------------', 'weight':257.095},
    'Md' :{'number':101, 'name':'-------------', 'weight':258.100},
    'No' :{'number':102, 'name':'-------------', 'weight':259.101},
    'Lr' :{'number':103, 'name':'-------------', 'weight':262.000},
    'Rf' :{'number':104, 'name':'-------------', 'weight':261.000},
    'Db' :{'number':105, 'name':'-------------', 'weight':262.000},
    'Sg' :{'number':106, 'name':'-------------', 'weight':266.000},
    'Bh' :{'number':107, 'name':'-------------', 'weight':264.000},
    'Hs' :{'number':108, 'name':'-------------', 'weight':269.000},
    'Mt' :{'number':109, 'name':'-------------', 'weight':268.000},
    'Ds' :{'number':110, 'name':'-------------', 'weight':269.000},
    'Rg' :{'number':111, 'name':'-------------', 'weight':272.000},
    'Cn' :{'number':112, 'name':'-------------', 'weight':277.000},
    'Uut':{'number':113, 'name':'-------------', 'weight': -1.000},
    'Fl' :{'number':114, 'name':'-------------', 'weight':289.000},
    'Uup':{'number':115, 'name':'-------------', 'weight': -1.000},
    'Lv' :{'number':116, 'name':'-------------', 'weight':298.000},
    'Uus':{'number':117, 'name':'-------------', 'weight': -1.000},
    'Uuo':{'number':118, 'name':'-------------', 'weight': -1.000},
}

def get_symbol_from_number(number):
    if number > len(ELEMENTS):
        raise ValueError('No element of number {0:d}'.format(number))
    for k,v in ELEMENTS.items():
        if v['number'] == number:
            return k
    raise ValueError('No element found for number {0:d}'.format(number))


def get_number_from_symbol(symbol):
    try:
        element = ELEMENTS[symbol]
    except:
        raise ValueError('No symbol {0:s}'.format(symbol))
    return element['number']


def write_info_of_element(*symbols):
    for s in symbols:
        try:
            element = ELEMENTS[s]
        except:
            raise ValueError('No symbol {0:s}'.format(s))
        print 'Element: {0:s}'.format(s)
        print '  name: {0:s}'.format(element['name'])
        print '  atomic number: {0:d}'.format(element['number'])
        print '  weight: {0:8.3f}'.format(element['weight'])
        print ''


class Atom(object):
    """
    Define attributions of an atom.
    """

    def __init__(self):
        self.pos= np.zeros((3,))
        self.vel= np.zeros((3,))
        self.id= 0
        self.ifmv= 1
        self.sid= 1
        self.symbol = get_symbol_from_number(self.sid)

    def set_pos(self,x,y,z):
        self.pos= np.array([x,y,z],dtype=float)

    def set_vel(self,x,y,z):
        self.vel= np.array([x,y,z],dtype=float)

    def set_sid(self,sid):
        self.sid= int(sid)

    def set_symbol(self,symbol):
        self.symbol = symbol

    def set_id(self,id):
        self.id= id

    def tag(self):
        tag= self.sid +self.ifmv*0.1 +self.id*1e-14
        return tag

    def decode_tag(self,tag):
        self.sid= int(tag)
        self.ifmv= int((tag-self.sid)*10)
        self.id= int(((tag-self.sid)*10 -self.ifmv)*1e+14)

if __name__ == '__main__':
    
    args= docopt(__doc__)
    symbols = args['SYMBOL']
    write_info_of_element(*symbols)


