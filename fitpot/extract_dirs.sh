#!/bin/bash
#
# Extract directories which include OUTCAR file under the working directory,
# and cp it to the directories named in the order, in 5 digits.
#

n=0
for dir in $(find . -name "OUTCAR")
do 
    dirname=$(echo $dir | sed 's/\/OUTCAR//')
    n=$(expr $n + 1)
    distname=$(printf "%05d\n" $n)
    echo "making $distname"
    mkdir "$distname"
    cp -r "$dirname" "./$distname/vasp"
done
