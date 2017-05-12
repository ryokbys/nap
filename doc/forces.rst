.. Manual for force fields implemented in NAP

========================================
Force fields
========================================


.. _neural_network:

Neural-network (NN)
========================================

NN potential requires following two input files in the working directory.

- ``in.const.NN``: network structure, types of symmetry functions, parameters of the symmetry functions, and interacting pairs.
- ``in.params.NN``: number of weights, cutoff radii for two- and three-body interactions, and weight values of the network.

These two files must be consistent such that the number of weights must correspond to the number of symmetry functions, number of layers, and number of nodes in each layer.
The examples of these files can be found in ``pmd/force_params/NN_??????`` directories.

There are some rules in ``in.const.NN`` as follows:

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

