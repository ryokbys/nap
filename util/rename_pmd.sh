#!/bin/sh


for file in pmd*-"`printf \"%03d\" $1`"
do
  echo "$file ---> `echo $file | sed 's/-.../-000/'`"
  cp $file "`echo $file | sed 's/-.../-000/'`"
done

