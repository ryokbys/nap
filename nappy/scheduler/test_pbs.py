# -*- coding: utf-8 -*-
import os,sys
import unittest

import nappy.scheduler.pbs as pbs

__copyright__ = u""
__version__ = "170127"
__authors__ = "Ryo KOBAYASHI"
__email__   = "kobayashi.ryo@nitech.ac.jp"


text_qstat_to_test_short = """
Job id                    Name             User            Time Use S Queue
------------------------- ---------------- --------------- -------- - -----
114618.mike               bulk             kasugai         20:45:05 R batch          
114620.mike               pmd              kobayashi              0 R default  
"""

text_qstat_to_test = """
Job Id: 114618.mike.bw.nitech.ac.jp
    Job_Name = bulk
    Job_Owner = kasugai@mike.bw.nitech.ac.jp
    resources_used.cput = 20:49:35
    resources_used.mem = 16148kb
    resources_used.vmem = 329084kb
    resources_used.walltime = 20:50:27
    job_state = R
    queue = batch
    server = mike.bw.nitech.ac.jp
    Checkpoint = u
    ctime = Thu Jan 26 19:25:59 2017
    Error_Path = mike.bw.nitech.ac.jp:/home/kasugai/LLTO/test2000/num03/bulk.e
	114618
    exec_host = cnode025.bw.nitech.ac.jp/0
    Hold_Types = n
    Join_Path = oe
    Keep_Files = n
    Mail_Points = a
    mtime = Thu Jan 26 19:26:00 2017
    Output_Path = mike.bw.nitech.ac.jp:/home/kasugai/LLTO/test2000/num03/bulk.
	o114618
    Priority = 0
    qtime = Thu Jan 26 19:25:59 2017
    Rerunable = True
    Resource_List.cput = 7200:00:00
    Resource_List.nodect = 1
    Resource_List.nodes = 1:ppn=1:E5v3
    Resource_List.walltime = 720:00:00
    session_id = 16750
    etime = Thu Jan 26 19:25:59 2017
    submit_args = auto_LLTO_MD.sh
    start_time = Thu Jan 26 19:26:00 2017
    Walltime.Remaining = 2516937
    start_count = 1

Job Id: 114621.mike.bw.nitech.ac.jp
    Job_Name = pmd
    Job_Owner = kobayashi@mike.bw.nitech.ac.jp
    job_state = R
    queue = default
    server = mike.bw.nitech.ac.jp
    Checkpoint = u
    ctime = Fri Jan 27 16:17:00 2017
    Error_Path = mike.bw.nitech.ac.jp:/home/kobayashi/test/pmd.e114621
    exec_host = cnode002.bw.nitech.ac.jp/1+cnode002.bw.nitech.ac.jp/0+cnode003
	.bw.nitech.ac.jp/1+cnode003.bw.nitech.ac.jp/0
    Hold_Types = n
    Join_Path = oe
    Keep_Files = n
    Mail_Points = a
    mtime = Fri Jan 27 16:17:01 2017
    Output_Path = mike.bw.nitech.ac.jp:/home/kobayashi/test/out
    Priority = 0
    qtime = Fri Jan 27 16:17:00 2017
    Rerunable = True
    Resource_List.cput = 7200:00:00
    Resource_List.nodect = 2
    Resource_List.nodes = 2:ppn=2
    session_id = 28089
    Variable_List = PBS_O_QUEUE=default,PBS_O_HOME=/home/kobayashi,
	PBS_O_LANG=UTF-8,PBS_O_LOGNAME=kobayashi,
	PBS_O_PATH=/home/kobayashi/bin:/home/kobayashi/local/bin:/home/kobaya
	shi/local/phonopy-1.9.0.1/bin:/usr/local/bin:/home/kobayashi/p4vasp/bi
	n:/opt/bin:/opt/vasp/bin:/usr/local/test/bin:/usr/lib64/qt-3.3/bin:/us
	r/local/torque-2.3.7/bin:/opt/intel/impi/4.0.0.028/intel64/bin:/opt/in
	tel/Compiler/11.1/073/bin/intel64:/opt/intel/Compiler/11.1/073/bin/int
	el64:/usr/local/bin:/bin:/usr/bin:/usr/local/maui/sbin:/usr/local/maui
	/bin:/home/kobayashi/src/ase/tools,
	PBS_O_MAIL=/var/spool/mail/kobayashi,PBS_O_SHELL=/bin/bash,
	PBS_SERVER=mike.bw.nitech.ac.jp,PBS_O_HOST=mike.bw.nitech.ac.jp,
	PBS_O_WORKDIR=/home/kobayashi/test
    etime = Fri Jan 27 16:17:00 2017
    submit_args = batch.sleep.sh
    start_time = Fri Jan 27 16:17:01 2017
    start_count = 1
"""


class TestPBS(unittest.TestCase):

    def test_joblist_command(self):
        command = pbs.get_joblist_command()
        print('command: '+command)
        self.assertTrue('qstat' in command)

    def test_parse_jobdata(self):
        jobdata = pbs.parse_jobdata(text_qstat_to_test)
        # print(jobdata)
        self.assertTrue( len(jobdata) != 0 )

        job = jobdata[0]
        self.assertTrue( 'Job Id' in job )
        self.assertTrue( 'Job_Name' in job )
        self.assertTrue( 'Job_Owner' in job )
        self.assertTrue( 'Resource_List.nodes' in job )

if __name__ == '__main__':

    unittest.main()
