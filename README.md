<p align="center">
  <img width="256" src="./doc/figs/logo_nap_256.npg">
</p>

# What's NAP
**Nagoya Atomistic-simulation Package (NAP)** consists of several programs and utilities:
- parallel molecular dynamics simulation (*pmd*)
- potential parameter fitting (*fitpot*)
- python modules for pre/post-processes (*nappy*)
The program, *pmd* , includes various interatomic potentials for metals and semiconductors,
and uses spatial decomposition technique for the parallelization, and cell-list method for efficient neighbor search.

# Who made this?
* Ryo KOBAYASHI
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology. (Nov, 2017)

# Compilation and usage

For the short test, whether or not you can use this program in your environment,

```bash
$ ./configure --prefix=$(pwd)
$ cd pmd/
$ make
$ cd ../example/test_W/
$ ../../pmd/pmd
```

If it works, you can use this program.
For details, please see the [manual](http://ryokbys.web.nitech.ac.jp/contents/pmd_usage) or ask me via e-mail (kobayashi.ryo[at]nitech.ac.jp).

# Acknowledgements
This program is supported in part by ["Materials research by Information Integration" Initiative (MI2I)](http://www.nims.go.jp/MII-I/) project of the Support Program for Starting Up Innovation Hub from Japan Science and Technology Agency (JST).


# LICENSE
This software is released under the MIT License, see the LICENSE.

