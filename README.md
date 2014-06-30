# What's NAP
Nagoya Atomistic-simulation Package (NAP) consists of several programs that perform parallel molecular dynamics simulation (pmd),
hybrid quantum-mechanical/classical-mechanical simulation (qmcl), nudged elastic band method (neb),
and potential parameter fitting (fitpot).
The program, pmd, includes various interatomic potentials for metals and semiconductors,
and uses spatial decomposition technique for the parallelization, and cell-list method for efficient neighbor search.

# Who made this?
* Ryo KOBAYASHI
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology. (2014-01-07)

# Compilation and usage
See the manual web site below (but Japanese only),
http://locs.bw.nitech.ac.jp/~kobayashi/pmd_manual

For the short test, whether or not you can use this program in your environment,

```bash
$ ./configure --prefix=$(pwd)
$ cd pmd/
$ make
$ cd ../sample
$ ../pmd/pmd
```

If it works, you can use this program.
For details, please see the manual (but Japanese only) or ask me via e-mail.

# LICENSE
This software is released under the MIT License, see LICENSE.
