# -*- coding: utf-8 -*-
from __future__ import print_function

import os,sys
import unittest

import nappy.sheduler.fujitsu as fujitsu

__copyright__ = u""
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "170125"
__authors__ = "Ryo KOBAYASHI"
__email__   = "kobayashi.ryo@nitech.ac.jp"

text_squeue_to_test = """
  ACCEPT QUEUED  STGIN  READY RUNING RUNOUT STGOUT   HOLD  ERROR   TOTAL
       0     86      0      0     15      0      0      0      0     101
s      0     86      0      0     15      0      0      0      0     101

 JOB ID                      : 390031
 JOB NAME                    : y0.140
 STATE                       : RUN
 USER                        : z41110v
 NODE NUM (REQUIRE)          : 2:noncont
 CPU NUM (REQUIRE)           : 64 <DEFAULT>
 ELAPSE TIME (LIMIT)         : 02:00:00 (7200)
 ELAPSE TIME (USE)           : 01:50:05 (6605)
 JOB START DATE              : 2016/07/29 21:11:29<
 ACCEPT DATE                 : 2016/07/29 20:55:56
 COMMENT                     : 
 RESOURCE GROUP              : fx-small
 NODE ID (USE)               : 0xFFE60008 0xFFE60005

 JOB ID                      : 390033
 JOB NAME                    : y0.160
 STATE                       : RUN
 USER                        : z41110v
 NODE NUM (REQUIRE)          : 2:noncont
 CPU NUM (REQUIRE)           : 64 <DEFAULT>
 ELAPSE TIME (LIMIT)         : 02:00:00 (7200)
 ELAPSE TIME (USE)           : 01:42:19 (6139)
 JOB START DATE              : 2016/07/29 21:19:15<
 ACCEPT DATE                 : 2016/07/29 20:55:56
 COMMENT                     : 
 RESOURCE GROUP              : fx-small
 NODE ID (USE)               : 0xFFE90006 0xFFE90002

 JOB ID                      : 390047
 JOB NAME                    : y0.270
 STATE                       : QUE
 USER                        : z41110v
 NODE NUM (REQUIRE)          : 2:noncont
 CPU NUM (REQUIRE)           : 64 <DEFAULT>
 ELAPSE TIME (LIMIT)         : 02:00:00 (7200)
 ELAPSE TIME (USE)           : 00:00:00 (0)
 JOB START DATE              : (2016/07/29 23:20:00)<
 ACCEPT DATE                 : 2016/07/29 20:55:57
 COMMENT                     : 
 RESOURCE GROUP              : fx-small
 NODE ID (USE)               : 0xFFE80006 0xFFE80002

"""

class TestFujitsuFX(unittest.TestCase):
    """
    Test some functions defined in fujitsu.py.
    """

    def test_joblist_command(self):
        command = fujitsu.get_joblist_command()
        print('command: '+command)
        self.assertTrue('pjstat' in command)

    def test_extract_jobdata(self):
        jobdata = fujitsu.extract_jobdata(text_squeue_to_test)
        # print('jobdata:')
        self.assertTrue( len(jobdata) != 0 )
        
        job = jobdata[0]
        self.assertTrue( 'JOB ID' in job )
        self.assertTrue( 'STATE' in job )
        self.assertTrue( 'USER' in job )
        self.assertTrue( 'NODE NUM' in job )
        self.assertTrue( 'CPU NUM' in job )
        self.assertTrue( 'RESOURCE GROUP' in job )

if __name__ == '__main__':        
    unittest.main()
