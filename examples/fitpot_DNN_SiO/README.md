# nap/example/fp_LZP

This example shows how to use `fitpot` to optimize the parameters in neural-network potential for Si-O system, see the Ref. [1,2] for details of NN potential. The files required to perform `fitpot` are:

- `in.fitpot` -- information about `fitpot` setting.
- `in.params.desc` -- information about descriptors used as inputs for the NN potential.
- `in.params.DNN` -- initial values of potential parameters (NN weights) and search ranges of the parameters
- `in.params.ZBL` -- classical potential added to the NN potential (in this case, ZBL potential)
- `dataset/smpl_XXX/` -- reference dataset

To perform `fitpot` using 2 MPI processes,
```bash
$ mpirun -np 2 /path/to/fitpot | tee out.fitpot
```

If the tail of output shows like the following, at least `fitpot` program finished corectly without errors.
```bash
$ tail out.fitpot
 ENERGY:      100          16.58    0.0000479    0.0000000    0.0000576    0.0000000    1.0000000    0.0000000
 FORCE:       100          16.59    0.0020653    0.0000000    0.0049990    0.0000000    0.9999964    0.0000000
 STRESS:      100          16.59   39.5048633    0.0000000   98.5031577    0.0000000  -99.7924332    0.0000000
 Number of func and grad calls =  203  101
 Memory/proc =         0.420 MB
 Time func =           1.423 sec
 Time grad =          14.367 sec
 Time comm =           0.001 sec
 Time      =          16.588 sec  =   0h00m16s
 Job finished at 10:25:05 on 2020-08-20
```
Please see the nap documentation for more details how to use `fitpot` to obtain the NN potential parameters.

## References

1. Behler, Jörg, and Michele Parrinello. 2007. “Generalized Neural-Network Representation of High-Dimensional Potential-Energy Surfaces.” Physical Review Letters 98 (14): 146401–146401
2. Kobayashi, Ryo, Daniele Giofré, Till Junge, Michele Ceriotti, and William Arthur Curtin. 2017. “Neural Network Potential for Al-Mg-Si Alloys.” Physical Review Materials 1 (5): 53604–11.[link](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.053604)
