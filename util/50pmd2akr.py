#!/usr/bin/python
#-----------------------------------------------------------------------
# Convert pmd data files to AKIRA data format.
#-----------------------------------------------------------------------
# Usage
# -----
#  At the directory where there is output directory.
#  $ ./50pmd2akr.py ### ?1 ?2 ?3 ?4 ?5 ...
#      ###  - output number (ccount)
#      ?#   - column number to be included in akr file
#

import sys,glob,re,math

usage=' USAGE: ./50pmd2akr.py 10 8 9 10 11 12 13 14 15 16'
#####  main routine  ###################################################
if len(sys.argv) < 3:
    raise StandardError, usage

if int(sys.argv[1]) >= 1000:
    raise StandardError," # must be less than 1000."

if len(sys.argv) > 11:
    raise StandardError, ' Num of auxiliary data should be less than 9'

nauxdat= len(sys.argv) -2

# iadat= -1
# if len(sys.argv)==3:
#     iadat= int(sys.argv[2])
#     if iadat > 9:
#         raise StandardError," argv[2] must be <= 9"

idfile= int(sys.argv[1])

outfile= open("akr%03d" % idfile, 'w')

#=====Total num of particles
ntot= 0
i=0
nadat= 0
for file in glob.glob('./pmd???-%03d' % idfile ):
    infile= open(file,'r')
    data= infile.readline().split() # 1st
    np= int(data[0])
    #=====If it first file, get system size and check num of auxiliary data
    #aa= zeros((3,3))
    if i == 0:
        data= infile.readline().split() # 2nd
        a1= [float(data[0]), float(data[1]), float(data[2])]
        data= infile.readline().split() # 3rd
        a2= [float(data[0]), float(data[1]), float(data[2])]
        data= infile.readline().split() # 4th
        a3= [float(data[0]), float(data[1]), float(data[2])]
        data= infile.readline().split() # 5th
        data= infile.readline().split() # 6th
        data= infile.readline().split() # 7th
        data= infile.readline().split() # 8th line, 1st atom data
        nadat= len(data) -4
    #    i += 1
    infile.close()
    ntot += np
print " Total num of atoms = %d" % ntot
print " Num of auxiliary data = %d" % nadat
print " Designated num of aux. data = %d" % nauxdat
if nauxdat > nadat:
    raise StandardError, ' nauxdat > nadat'
# if iadat > nadat:
#     raise StandardError, " iadat > nadat"
# if len(sys.argv)==3:
#     nadat= iadat
print " aa:"
print " %15.7e %15.7e %15.7e" % (a1[0],a1[1],a1[2])
print " %15.7e %15.7e %15.7e" % (a2[0],a2[1],a2[2])
print " %15.7e %15.7e %15.7e" % (a3[0],a3[1],a3[2])

print >>outfile," %8d  %2d  0  0" % (ntot,nauxdat)
print >>outfile," %15.7e %15.7e %15.7e" % (a1[0],a1[1],a1[2])
print >>outfile," %15.7e %15.7e %15.7e" % (a2[0],a2[1],a2[2])
print >>outfile," %15.7e %15.7e %15.7e" % (a3[0],a3[1],a3[2])


#=====All atoms
pattern="pmd."
for file in glob.glob('./pmd*-%03d' % idfile ):
    print " %s ---> akr%03d" % (file,idfile)
    infile= open(file,'r')
    data= infile.readline().split() # 1st line
    np= int(data[0])
    data= infile.readline().split() # 2nd line
    data= infile.readline().split() # 3rd line
    data= infile.readline().split() # 4th line
    data= infile.readline().split() # 5th line
    data= infile.readline().split() # 6th line
    data= infile.readline().split() # 7th line
    for i in range(np):
        data= infile.readline().split()
        #print " i,len(data)= %d %d" % (i,len(data))
        for j in range(len(data)):
            data[j]= float(data[j])
        outfile.write(" %2d %12.4e %12.4e %12.4e" % 
                      (int(data[0]),data[1],data[2],data[3]) )
        for j in range(nauxdat):
            k= int(sys.argv[2+j])
            val= data[k]
            #.....Checking 'NAN'
            if( val != val ):
                val= 0.0
                print " val != val"
            outfile.write(" %12.4e" % val )
        outfile.write("\n")

outfile.close()
