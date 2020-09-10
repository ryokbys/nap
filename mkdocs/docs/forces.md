# Force fields

Force fields (FFs) are specified at `force_type` keyword in `in.pmd`
file. Plural FFs can be specified as a space-separated list,

    force_type    NN Morse Coulomb

Each FF reads one or some parameter files in the working directory,
which are usually named like `in.params.Coulomb` or so, specified by FF.
For example, **DNN** requires two files, `in.params.DNN` (file for
weights in NN) and `in.params.desc` (file for descriptor information).

Available FFs are listed below:

- [force_fields/DNN](force_fields/DNN.md)
- [force_fields/Coulomb](force_fields/Coulomb.md)
- [force_fields/Morse](force_fields/Morse.md)
- [force_fields/SW](force_fields/SW.md)
- Lennard-Jones (Ar)
- Linear regression
