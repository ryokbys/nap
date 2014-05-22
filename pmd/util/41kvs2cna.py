#!/usr/bin/python
#-----------------------------------------------------------------------
# Rename 'fin###' files to 'ini###'
#-----------------------------------------------------------------------
# Usage:
#   ./41kvs2cna.py

import sys,glob,os

if not os.path.exists("cna"):
    print " error: cna does not exit!!"
    sys.exit()

kvs= "kvs" 
cna= "cna"
for fname in glob.glob(kvs+"[0123456789]*"):
    newfname= cna+fname[3:6]
    print fname," ---> ",newfname
    os.system("./cna < %s > %s" % (fname,newfname))
