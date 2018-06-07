.. _neural_network:

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


``in.params.desc``
----------------------

There are some rules in ``in.params.desc`` as follows:

- Three-body terms must come after all the two-body terms.
- Within two-body terms, terms for the specific pair must be written in consecutive lines such as
  ::

     1   1   1   10.000     3.0000
     1   1   1   10.000     4.0000
     1   1   2   10.000     3.0000
     1   1   2   10.000     4.0000
     1   1   2   10.000     5.0000
     1   1   3   10.000     3.0000

  Not like
  ::

     1   1   1   10.000     3.0000
     1   1   2   10.000     3.0000
     1   1   3   10.000     3.0000
     1   1   1   10.000     4.0000
     1   1   1   10.000     5.0000
     1   1   2   10.000     4.0000

  where the pairs 1-2 and 1-3 (see 2nd and 3rd columns) appear before finishing all the inputs of 1-1 pair.


``in.params.NN2``
------------------------

``in.params.NN2`` file should have the following format.
::

    1   18   10
    -3.64106023330479E-001 -1.0000E-01  1.0000E-01
    -2.01340565152879E+000 -1.0000E-01  1.0000E-01

where three digits in the 1st line are *number of layers*, *number of input nodes*, and *number of nodes in the 1st layer*.
There should be 190 (= 18*10 + 10) following lines (in this case), with *NN weight* and two dummy values.

