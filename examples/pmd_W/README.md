# nap/examples/pmd_W

This example shows how to perform an MD simulation of bcc-W system using the Ito potential[1]. The files required to run the MD simulation in this case are:

- `pmdini` -- information about simulation cell and positions of atoms
- `in.pmd` -- information about MD simulation setting

To perform MD,
```bash
$ /path/to/pmd | tee out.pmd
```

Users can check the result by looking at the output `out.pmd` and compare it with `out.pmd.REF` which is the reference output. Users may be compare them by looking at `Potential` and confirm that the results are identical.
```bash
$ grep 'Potential' out.pmd*
out.pmd:   Potential energy=       -463.72014 eV =     -8.587 eV/atom
out.pmd:   Potential energy=       -460.77567 eV =     -8.533 eV/atom
out.pmd.REF:   Potential energy=       -463.72014 eV =     -8.587 eV/atom
out.pmd.REF:   Potential energy=       -460.77567 eV =     -8.533 eV/atom
```


## References

1. A.M. Ito, Y. Yoshimoto, S. Saito, A. Takayama, H. Nakamura, Phys. Scripta T159, 014062 (2014).
