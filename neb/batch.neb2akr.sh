#!/bin/bash
#
# USAGE:
#   ./batch.neb2akr.sh
#

for file in neb[0-9][0-9][0-9]
do
  neb2akr < $file > `echo $file | sed 's/neb/akr/'`
done
