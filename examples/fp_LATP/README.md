# nap/example/fp_LZP

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

To perform `fp.py` using four threads,
```bash
$ /path/to/fp.py --nproc=4 | tee out.fp
```

Since the `fp.py` uses random values to generate trial parameters and the random values would be different every time it is executed, the output obtained will not be identical to that in `out.fp.REF`. But it will be OK if `tail out.fp` looks like the following,
```bash
$ tail out.fp
   iid,losses=        6     1.2219     3.3689     0.0120     4.60281
   iid,losses=        9     0.6474     7.5324     0.1334     8.31316
 step,time,best,vars=      1    173.1    2.2234  1.000  0.748  0.783  1.019  1.071  0.689  0.883  1.629  1.829  1.615  1.951  2.442  1.606  1.739  1.791  3.519
   iid,losses=       11     1.7628     0.2426     2.4705     4.47583
   iid,losses=       12     2.0236     1.5161     0.7651     4.30480
   iid,losses=       10     1.2067     0.6369     0.3003     2.14386
   iid,losses=       13     0.6055     9.3178     0.1672    10.09049
   iid,losses=       14     1.0248     8.7687     0.0168     9.81031
 step,time,best,vars=      2    346.8    2.1439  1.000  0.785  0.785  1.030  1.090  0.687  0.911  1.627  1.825  1.615  1.957  2.442  1.653  1.770  1.790  3.429
elapsed time = 346.771660 sec.
```
And when you use `fp.py`, you had better use more processes than 4 like in this case to efficiently run the program.

## Differences from `examples/fp_LZP`

This example used **whatever mode** for target quantities not like **DF-matching mode** in `examples/fp_LZP`. Differences from `examples/fp_LZP` include:

- The format of `data.ref.xxx` files
- `match` entry in `in.fitpot`
- `subjob.sh` script (`--out4fp` option is specified in this case.)


## References

1. R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, “High-throughput production of force-fields for solid-state electrolyte materials”, APL Materials 8, 081111 (2020). [link](https://aip.scitation.org/doi/10.1063/5.0015373).
