.. _Coulomb:

Coulomb potential (Ewald or short/screend)
===========================================

Coulomb potential requires a parameter file ``in.params.Coulomb`` that has the following format.

::

   charges  fixed
     Si   1.0
     O   -2.0
   interactions
     Si  Si
     Si  O
     O   O
   terms  short
   
   sigma  2.5


In this format, black lines are neglected.
There are some keywords:

``charges`` : ``fixed`` or ``variable`` (or ``qeq``)
  Followed by the lines of species charges, e.g., :math:`q_1 = 1.0`  and :math:`q_2 = -2.0`.

``interactions`` : *optional*
  Followed by the pairs of species. If not specified, all the interactions are taken into account.

``terms`` : ``full``, ``long`` or ``short``/``screened_cut``
  Either full Ewald terms, long-range term only, short-range term only, 
  or short-term with smooth cutoff.

``sigma`` : 
  Width of the Gaussian charge distribution, which is related to the accuracy in case of Ewald method.


Variable charge or QEq
----------------------------

Coulomb potential can treat **variable charge** or **QEq** by specifying ``variable`` or ``qeq`` to the ``charges`` keyword as shown below.

::

   charges  variable
     Si  4.7695    8.7893    0.0   0.0  2.4
     O   7.5405    15.8067   0.0  -1.2  0.0
   interactions
     Si  Si
     Si  O
     O   O
   terms  short
   sigma  2.5
   conv_eps  1.0d-6

Here, ``charges variable`` requires some following lines that have

::

   name,  chi,  Jii,  E0,  qlow,  qup

- ``name``: name of the chemical species
- ``chi``: electronegativity of the species (eV)
- ``Jii``: hardness of the species (eV)
- ``E0``: atomic energy (eV)
- ``qlow``: lower limit of the charge of the species
- ``qup``: upper limit of the charge of the species

``conv_eps``
  Convergence criterion of the charge optimization. 
  Usually it should be very small to achieve good energy conservation.

