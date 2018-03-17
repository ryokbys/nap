.. _SW:

Stillinger-Weber (SW) potential
================================

The SW potential requires a parameter file ``in.params.SW`` of the following format.

::

   unit  2.1678  2.0951
   1  1  7.049556277  0.6022245584  4.0  0.0  1.0  1.8
   1  1  1  21.0  1.20

- Line starting with ``unit`` defines the units of energy and length.
- Line with 8 entries is for two-body parameters with the format,
  ::

     isp, jsp, Eij, Aij, Bij, pij, qij, cij, rcij

- Line with 5 entries is for three-body parameters with the format,
  ::

     isp, jsp, ksp, sij, tij

