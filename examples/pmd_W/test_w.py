#!/usr/bin/env python
"""
Test pmd run about Ito3 potential for W-He.

Usage:
  test_w.py [options]

Options:
  -h, --help  Show this message and exit.
"""
import os
from docopt import docopt
import unittest

__author__ = "RYO KOBAYASHI"
__version__ = "200723"

def get_init_epot(fname='out.pmd'):
    with open(fname,'r') as f:
        for line in f.readlines():
            if 'Potential energy' in line:
                dat = line.split()
                epot = float(dat[2])
                break
    return epot


class TestMD(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self._dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(self._dir)
        print('cd ../../pmd && make pmd')
        os.chdir('../../pmd')
        os.system('make pmd')
        os.chdir(self._dir)
        cmd = '../../pmd/pmd > out.pmd 2>&1'
        print(cmd)
        os.system(cmd)
        
    def test_pmd(self):
        os.chdir(self._dir)
        eref = get_init_epot(fname='out.pmd.REF')
        epot = get_init_epot(fname='out.pmd')
        print('eref = {0:24.14e}'.format(eref))
        print('epot = {0:24.14e}'.format(epot))
        ediff = abs(eref-epot)
        print('ediff= {0:24.14e}'.format(ediff))
        self.assertLess(ediff,1.0e-10)


if __name__ == "__main__":

    args = docopt(__doc__)

    unittest.main()
    
