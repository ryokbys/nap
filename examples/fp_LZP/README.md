# nap/example/fp_LZP

This example shows how to use `fp.py` to optimize the parameters in interatomic potential for Li-Zr-P-O system. The potential forms are screened Coulomb, Morse and angular and the target quantities are RDF, ADF and equilibrium volume, see the Ref. [1] for details. The files required to perform `fp.py` are:

- `in.fitpot` -- information about `fp.py` setting
- `in.vars.fitpot` -- initial values of potential parameters and search ranges of the parameters
- `in.params.Coulomb` -- parameter file for Coulomb potential
- `data.ref.xxx` -- target quantities (RDF, ADF and volume)
- `subjob.sh` -- shell-script that describe what to compute target quantities with the current parameters and
- `in.pmd.NpT`, `pmdini`, and something -- files need to perform MD simulations described in `subjob.sh`

To see command-line options,
```bash
$ python /path/to/fp.py -h
```

To perform `fp.py` using four threads,
```bash
$ /path/to/fp.py --nproc=4 | tee out.fp
```

Since the `fp.py` uses random values to generate trial parameters and the random values would be different every time it is executed, the output obtained will not be identical to that in `out.fp.REF`. But it will be OK if `tail out.fp` looks like the following,
```bash
$ tail out.fp
   iid,Lr,Lth,Lvol,Llat,L=        1    1.3376     0.7049     7.7690     0.0000     9.8116
   iid,Lr,Lth,Lvol,Llat,L=        3    0.8912     0.6027     4.7938     0.0000     6.2877
 step,time,best,vars=      0     40.3    6.2877  1.000  1.142  1.255  1.054  0.868  1.208  2.106  1.778  2.988  2.044  1.921  4.330  2.013  1.536  2.486  1.000
   iid,Lr,Lth,Lvol,Llat,L=        8    3.1944     0.7949    25.6610     0.0000    29.6503
   iid,Lr,Lth,Lvol,Llat,L=        7    1.8878     0.7908    19.5151     0.0000    22.1937
   iid,Lr,Lth,Lvol,Llat,L=        6    1.2444     0.7123     7.4910     0.0000     9.4477
   iid,Lr,Lth,Lvol,Llat,L=        5    0.8285     0.6159     5.3013     0.0000     6.7457
   iid,Lr,Lth,Lvol,Llat,L=        9    2.0158     0.7721    15.0153     0.0000    17.8033
 step,time,best,vars=      1    104.3    6.2877  1.000  1.142  1.255  1.054  0.868  1.208  2.106  1.778  2.988  2.044  1.921  4.330  2.013  1.536  2.486  1.000
elapsed time = 104.290595 sec.
```
And when you use `fp.py`, you had better use more processes than 4 like in this case to efficiently run the program.

## References

1. R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, “High-throughput production of force-fields for solid-state electrolyte materials”, APL Materials 8, 081111 (2020). [link](https://aip.scitation.org/doi/10.1063/5.0015373).
