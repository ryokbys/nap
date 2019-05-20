
.. index:: FAQ
.. _faq:

==================================================
FAQs
==================================================

1. :ref:`faq01`
2. :ref:`faq02`

------------------------

.. _faq01:

There are not enough slots available...
--------------------------------------------------

If you get the following error message when you run ``pmd`` with some number of parallel processes,

::

   --------------------------------------------------------------------------
   There are not enough slots available in the system to satisfy the 4 slots
   that were requested by the application:
     pmd
   
   Either request fewer slots for your application, or make more slots available
   for use.
   --------------------------------------------------------------------------

you may have to increase the number of slots available by the ``--oversubscribe`` option as,
::

   $ mpirun --oversubscribe -np 4 pmd

Reference:

* https://stackoverflow.com/questions/35704637/mpirun-not-enough-slots-available


-------

.. _faq02:

A system call failed during shared memory initialization...
============================================================

The following error message could appear at the end of ``pmd`` output, when the number of MPI processes and the number of automatically determined spatial decomposition divisions are different.

::

   --------------------------------------------------------------------------
   A system call failed during shared memory initialization that should
   not have.  It is likely that your MPI job will now either abort or
   experience performance degradation.
   
     Local host:  mbp-rk-eth.ogt.nitech.ac.jp
     System call: unlink(2) /var/folders/fh/cd90qtsj3d1fr3wmjw954xt40000gn/T//ompi.mbp-rk-eth.501/pid.52045/1/vader_segment.mbp-rk-eth.95f30001.2
     Error:       No such file or directory (errno 2)
   --------------------------------------------------------------------------

But the ``pmd`` results are correct. You can remove it by setting an environment variable as,
::

   export OMPI_MCA_btl=self,tcp


Reference:

* https://github.com/open-mpi/ompi/issues/6518



