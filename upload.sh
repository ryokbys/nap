#!/bin/bash

cwd=`dirname "${0}"`
if [ $# -ne 1 ]; then
  echo "There must be 1 argument." 1>&2
  exit 1
fi

src=$(cd $(dirname $0); pwd)
dst=$1
echo "${src} ==> ${dst}:src/nap/"

rsync -avz --exclude=".git/" --exclude="*.dSYM/" --exclude="subdir_*/" --exclude="doc/" --include="*/" --exclude="*.o" --exclude="*.mod" --exclude="*~" --exclude="makefile" --exclude="pmd" --exclude="fitpot" --exclude="*.pyc" --exclude="libpmd.a" ${src}/ ${dst}:src/nap/
#...Only selected makefiles are uploaded to the remote host
rsync -avz --include="examples/*/makefile" --include="makefile" --exclude="*" ${src}/ ${dst}:src/nap/
