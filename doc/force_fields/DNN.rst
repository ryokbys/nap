.. _DNN:

Deep neural-network (DNN) potential
========================================

.. note::

   The NN potentials (NN and NN2) were replaced with this ``DNN`` potential in January 2020.

DNN potential requires following two input files in the working directory.

- ``in.params.desc``: types of symmetry functions, parameters of the symmetry functions, their cutoff radii, and interacting pairs.
- ``in.params.DNN``: NN structures, number of layers, number of nodes for each layer, and weight values of the network.

These two files must be consistent such that the number of descriptors must correspond to the number of inputs in the NN.
The examples of these files can be found in ``pmd/force_params/DNN_??????`` directories.


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



in.params.DNN
------------------------

``in.params.DNN`` file should have the following format.
::

   !  sigtype: 2
   ! 
      3   20   10  10  5
    -3.64106023330479E-001 -1.0000E-01  1.0000E-01
    -2.01340565152879E+000 -1.0000E-01  1.0000E-01

- Lines starting with ``!`` are comment lines. If any of special keyword is found just after ``!``, an option will be passed to the program.
- The 1st entry of the 1st line following comments is the number of hidden layers in the NN.
- The 2nd entry is the number of nodes in 0-th layer, which is called input layer.
- The following digits are the number of nodes in hidden layers. The number of these digits must be the same as the number of hidden layers given by the 1st entry.
- Following lines include weight value, lower and upper bounds, which are only used in ``fitpot``.

