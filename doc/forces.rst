.. Manual for force fields implemented in NAP

========================================
Force fields
========================================

Force fields (FFs) are specified at ``force_type`` keyword in ``in.pmd`` file.
Plural FFs can be specified as a space-separated list,
::

   force_type    NN Morse Coulomb

Each FF reads one or some parameter files in the working directory, which are usually named like ``in.params.Coulomb`` or so, specified by FF. For example, **DNN** requires two files, ``in.params.DNN`` (file for weights in NN) and ``in.params.desc`` (file for descriptor information).

Available FFs are listed below:

- :doc:`force_fields/DNN`
- :doc:`force_fields/Coulomb`
- :doc:`force_fields/Morse`
- :doc:`force_fields/SW`
- Lennard-Jones (Ar)
- Linear regression
