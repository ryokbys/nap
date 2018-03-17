.. Manual for force fields implemented in NAP

========================================
Force fields
========================================

Available force fields are listed below:

- :doc:`force_fields/neural_network`
- :doc:`force_fields/Coulomb`
- :doc:`force_fields/Morse`
- Stillinger-Weber (Si)
- Lennard-Jones (Ar)


Force parameter files
=========================

Force fields (FF) usually require a few to a lot of parameters.
So each FF reads one or some parameter files in the working directory, which are usually named like ``in.params.Coulomb`` or so, specified by FF. For example, **NN** requires two files, ``in.const.NN`` and ``in.params.NN``.


