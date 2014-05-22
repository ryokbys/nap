#!/usr/bin/python


import os

for i in range(10):
    print i
    os.system("./pmd > out.pmd")
    os.system("cp out.stress out.stress.%02d" % i)
    os.system("cp pmd000-010 pmd000-000")
