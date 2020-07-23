<p align="center">
  <img width="256" src="doc/figs/logo_nap_256.png">
</p>

# What's nap
**Nagoya Atomistic-simulation Package (nap)** includes the following programs and utilities:
- parallel molecular dynamics simulation (*pmd*)
- potential parameter fitting (*fitpot* for neural-network potential and *fp.py* for other classical potentials)
- python modules for pre/post-processes (*nappy*)

The program, *pmd*, includes various interatomic potentials for metals and semiconductors,
and uses spatial decomposition technique for the parallelization, and linked-cell method for efficient neighbor search.

# Who made this?
* [Ryo KOBAYASHI](http://ryokbys.web.nitech.ac.jp/index.html)
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology.

# Compilation and usage

For the short test, whether or not you can use this program in your environment,

```bash
$ ./configure --prefix=$(pwd)
$ cd pmd/
$ make
$ cd ../examples/pmd_W/
$ ../../pmd/pmd
```

If it works, you can use this program.
For details, please see the [documentation](http://ryokbys.web.nitech.ac.jp/contents/pmd_usage) or ask me via e-mail (kobayashi.ryo[at]nitech.ac.jp).

# Acknowledgements
This program was supported in part by ["Materials research by Information Integration" Initiative (MI2I)](http://www.nims.go.jp/MII-I/) project of the Support Program for Starting Up Innovation Hub from Japan Science and Technology Agency (JST).


# LICENSE
This software is released under the MIT License, see the LICENSE.

