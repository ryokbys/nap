#!/usr/bin/env python
"""
Test code for clmgr.py.

Usage:
  test_clmgr.py [options]

Options:
  -h, --help  Show this message and exit.
"""
import os,sys
from docopt import docopt
import unittest

try:
    import nappy.clutil.clmgr as clmgr
except:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import clmgr as clmgr
    

__author__ = "RYO KOBAYASHI"
__version__ = "170123"

class TestClmgr(unittest.TestCase):

    # def setUp(self):
    #     self._obj = hoge

    def test_find_dirs_to_work(self):
        dirs = clmgr.find_dirs_to_work('.')
        # print(dirs)
        self.assertIsNotNone(dirs)

    def test_already_running(self):
        ans = clmgr.already_running()
        print('already_running? ',ans)
        self.assertIn(ans,(True,False))

    # def dearDown(self):
    #     self._obj.finalize()
        

if __name__ == "__main__":

    args = docopt(__doc__)

    unittest.main()
