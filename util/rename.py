#!/usr/bin/python
#-----------------------------------------------------------------------
# Rename 'fin###' files to 'ini###'
#-----------------------------------------------------------------------
# Usage:
#   ./renam.py ini fin

import sys,glob,os

if len(sys.argv)!= 3:
    print "usage: ./rename.py ini fin"
else:
    ini= sys.argv[1]
    fin= sys.argv[2]
    for fname in glob.glob(ini+"*"):
        newfname= fin+fname[3:6]
        print fname," ---> ",newfname
        os.rename(fname,newfname)
