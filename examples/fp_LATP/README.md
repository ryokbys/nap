# nap/example/fp_LATP

This example shows how to use `fp.py` to optimize the parameters in interatomic potential for Li-Al-Ti-P-O system. The potential forms are screened Coulomb, Morse and angular and the target quantities are RDF, ADF and equilibrium volume, see the Ref. [1] for details. The files required to perform `fp.py` are:

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

To perform `fp.py` using four threads with a specific random seed,
```bash
$ /path/to/fp.py --nproc=4 --random-seed=42 | tee out.fp
```

The output of `tail out.fp` command looks like the following:
```bash
❯ tail out.fp
   iid,losses=        3     0.9394     4.1495     0.0412     5.13005
   iid,losses=        2     0.8077     6.8559     0.0558     7.71932
 step,time,best_iid,best_loss,vars=      0     34.1     4   4.7795  1.000  1.156  0.799  0.854  0.761  0.781  1.015  1.823  2.068  1.523  1.517  2.208  1.346  1.962  1.702  5.038
   iid,losses=        5     0.8177     3.4517     0.0374     4.30683
   iid,losses=        7     1.0320     4.4855     0.0032     5.52068
   iid,losses=        6     0.9540     4.9692     0.0580     5.98121
   iid,losses=        8     0.8334     6.8761     0.0642     7.77373
   iid,losses=        9     0.8904     3.2130     0.0246     4.12806
 step,time,best_iid,best_loss,vars=      1     97.1     9   4.1281  1.000  1.161  1.191  0.778  0.751  0.802  1.242  2.560  1.465  1.717  1.603  2.669  1.921  1.992  2.128  3.879
elapsed time = 97.064136 sec.
```
And when you use `fp.py`, you had better use more processes than 4 like in this case to efficiently run the program.

## Differences from `examples/fp_LZP`

This example used **any-target mode** for target quantities not like **DF-matching mode** in `examples/fp_LZP`. Differences from `examples/fp_LZP` include:

- The format of `data.ref.xxx` files
- `match` entry in `in.fitpot`
- `subjob.sh` script (`--out4fp` option is specified in this case.)


## References

1. R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, “High-throughput production of force-fields for solid-state electrolyte materials”, APL Materials 8, 081111 (2020). [link](https://aip.scitation.org/doi/10.1063/5.0015373).
