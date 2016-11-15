#!/bin/bash

export PATH=/usr/bin:$PATH
time parallel --bar -j2 "mkdir -p {}/pmd; cp {}/pos {}/pmd/pmd0000; cp in.params.NN in.const.NN in.pmd {}/pmd/; cd {}/pmd/; ~/src/nap/pmd/pmd > out.pmd; cd ../.." ::: $*
