# nap/examples/pmd_DNN_SiO

This example shows how to perform an MD simulation of SiO2 system using a neural-network (NN) potential. The files required to run the MD simulation in this case are:

- `pmdini` -- information about simulation cell and positions of atoms
- `in.pmd` -- information about MD simulation setting
- `in.params.desc`, `in.params.DNN`, `in.params.ZBL` -- parameter files for the NN potential for Si-O.

To perform MD,
```bash
$ /path/to/pmd | tee out.pmd
```

Users can check the result by looking at the output `out.pmd` and compare it with `out.pmd.REF` which is the reference output. Users may be compare them by looking at `Potential` and confirm that the results are identical.
```bash
$ grep 'Potential' out.pmd*
out.pmd:   Potential energy=        -71.13849 eV =     -7.904 eV/atom
out.pmd:   Potential energy=        -70.97027 eV =     -7.886 eV/atom
out.pmd.REF:   Potential energy=        -71.13849 eV =     -7.904 eV/atom
out.pmd.REF:   Potential energy=        -70.97027 eV =     -7.886 eV/atom
```

And there are also some `dump_###` files that contain snapshot configurations of atoms and users can see the movie of the simulation using some visualization software such as [Ovito](https://www.ovito.org/about/).
