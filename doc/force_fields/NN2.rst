.. _NN2:

Neural-network (NN2) potential
========================================

.. note::

   The NN potential was updated to ver 2.0 on 7 June 2018 by the force-field name `NN2`.
   The document here is for the NN2 potential. If one has to use ver 1.0, see :doc:`neural_network`.

NN2 potential requires following two input files in the working directory.

- ``in.params.desc``: types of symmetry functions, parameters of the symmetry functions, their cutoff radii, and interacting pairs.
- ``in.params.NN``: NN structures, number of layers, number of nodes for each layer, and weight values of the network.

These two files must be consistent such that the number of weights must correspond to the number of symmetry functions, number of layers, and number of nodes in each layer.
The examples of these files can be found in ``pmd/force_params/NN2_??????`` directories.


in.params.desc
----------------------

.. note::

   The format of ``in.params.desc`` file has changed from the previous one at around May 2019, where the species are written directly by their acronym name such as *Si* not by species-ID (digit 1-9). If you want to use the previous-format ``in.params.desc``, you should modify it by replacing species-ID with species name.

The format of ``in.params.desc`` is like the following,
::

   2    20
   1  W   W   10.000    2.000
   1  W   W   10.000    3.000
   1  W   W   10.000    4.000
   1  W   H   10.000    2.000
   1  W   H   10.000    3.000
   1  W   H   10.000    4.000
   2  W   W   -0.900
   2  W   W   -0.800
   2  W   W   -0.700
   ...

- 1st line has two entries, *number of speceis* and *number of descriptors*.
- Following lines have each descriptor information, the 1st entry is the type of descriptor, 2nd and 3rd are species of interaction pair, from the 4th to the end are parameters of the descriptor. The number of parameters depend on the type of descriptor.



in.params.NN2
------------------------

``in.params.NN2`` file should have the following format.
::

    1   18   10
    -3.64106023330479E-001 -1.0000E-01  1.0000E-01
    -2.01340565152879E+000 -1.0000E-01  1.0000E-01

where three digits in the 1st line are *number of layers*, *number of input nodes*, and *number of nodes in the 1st layer*.
There should be 190 (= 18*10 + 10) following lines (in this case), with *NN weight* and two dummy values.

