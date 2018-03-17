.. _Coulomb:

Coulomb potential (Ewald or short/screend)
===========================================

Coulomb potential requires a parameter file ``in.params.Coulomb`` that has the following format.

::

   charges  fixed
     1   2.0
     2   -1.0
   interactions
    1  1
    1  2
    2  2
   terms  short
   
   sigma  2.5


In this format, black lines are neglected.
There are some keywords:

``charges`` : ``fixed`` or ``variable`` (or ``qeq``)
  Followed by the lines of species charges, e.g., :math:`q_1 = 2.0`  and :math:`q_2 = -1.0`.

``interactions`` : *optional*
  Followed by the pairs of species. If not specified, all the interactions are taken into account.

``terms`` : ``full``, ``long`` or ``short``/``screened``
  Either full Ewald terms, long-range term only, or short-range term only.

``sigma`` : 
  Width of the Gaussian charge distribution, which is related to the accuracy in case of Ewald method.


Variable charge or QEq
----------------------------

Coulomb potential can treat **variable charge** or **QEq** by specifying ``variable`` or ``qeq`` to the ``charges`` keyword as shown below.

::

   charges  variable
     1  Si  4.7695    8.7893    0.0   0.0  2.4
     2  O   7.5405    15.8067   0.0  -1.2  0.0
   interactions
    1  1
    1  2
    2  2
   terms  short
   sigma  2.5
   conv_eps  1.0d-6

Here, ``charges variable`` requires some following lines that have

::

   isp, name,  chi,  Jii,  E0,  qlow,  qup

- ``isp``: species-ID
- ``name``: species-name, which is just for human-readability
- ``chi``: electronegativity of the species (eV)
- ``Jii``: hardness of the species (eV)
- ``E0``: atomic energy (eV)
- ``qlow``: lower limit of the charge of the species
- ``qup``: upper limit of the charge of the species

``conv_eps``
  Convergence criterion of the charge optimization. 
  Usually it should be very small to achieve good energy conservation.

