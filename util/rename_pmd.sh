#!/bin/sh


for file in pmd*-"`printf \"%04d\" $1`"
do
  echo "$file ---> `echo $file | sed 's/-..../-0000/'`"
  cp $file "`echo $file | sed 's/-..../-0000/'`"
done

