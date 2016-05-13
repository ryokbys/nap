#!/bin/bash

time parallel --bar -j12 "cp in.params.NN in.const.NN in.smd {}/smd/; cd {}/smd/; ~/src/nap/pmd/smd > out.smd; cd ../.." ::: $*
