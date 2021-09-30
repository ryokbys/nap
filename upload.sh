#!/bin/bash

cwd=`dirname "${0}"`
if [ $# -ne 1 ]; then
  echo "There must be 1 argument." 1>&2
  exit 1
fi

src=$(cd $(dirname $0); pwd)
dst=$1
echo "${src} ==> ${dst}:src/nap/"

rsync -avz --exclude=".git/" --exclude="*.dSYM/" --exclude="doc/" --include="*/" --exclude="*.o" --exclude="*.mod" --exclude="*~" --exclude="makefile" --exclude="pmd" --exclude="fitpot" --exclude="*.pyc" ${src}/ ${dst}:src/nap/
#...Only the makefile at root dir will be uploaded to the remote host
rsync -avz --include="makefile" --exclude="*" ${src}/ ${dst}:src/nap/
