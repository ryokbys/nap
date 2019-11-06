==============================
nappy | NAP PYthon utilities
==============================

``nappy`` is set of python utilities. 
It is contained in ``nappy`` directory, but currenly it is not exactly a python module.
Users may have to use those utilities by calling them directly from the shell.


Setup
=======

The ``nappy`` should work with python-2.7 and python-3 series. But, as the maintenance of python-2.x series will stop in the near future, it is recommended to use python-3.

To use ``nappy`` in python program, it is required to add a path to ``nappy`` directory
to the environment variable ``PYTHONPATH``.
In case of ``bash``, you can achieve this by adding the following line to ``~/.bash_profile``,

.. code:: bash

   export PYTHONPATH=${PYTHONPATH}:/path/to/nap

You can check whether the path to ``nappy`` is added to ``PYTHONPATH`` by the following command,
::

   $ python -c 'import nappy; print nappy.__file__'


Quick start
===================

Once ``nappy`` is installed, do the following code on *ipython* or *jupyter notebook*,

.. code:: python

   import nappy.napsys as napsys, analyze
   nsys = napsys.NAPSystem(fname='/path/to/nap/example/test_W/pmdini')
   analyze(nsys)

Above code will show the following result.
::

   a1 vector = [    10.254,      0.000,      0.000]
   a2 vector = [     0.000,      9.613,      0.000]
   a3 vector = [     0.000,      0.000,      9.613]
   a =     10.254 A
   b =      9.613 A
   c =      9.613 A
   alpha =   90.00 deg.
   beta  =   90.00 deg.
   gamma =   90.00 deg.
   volume=    947.617 A^3
   number of atoms   =  54




Read and write files
==============================

To read atomic structures from files of several formats,

.. code:: python

   nsys = NAPSystem(fname='POSCAR')
   nsys = NAPSystem(fname='structure.dat', ffmt='xsf', specorder=['Al','O'])


Current available formats are:

* ``pmd``: file format specific for ``pmd`` program
* ``POSCAR``: VASP POSCAR file
* ``dump``: LAMMPS dump file
* ``xsf``
* ``akr``: file format for Akira viewer

And to write the structure to a file,

.. code:: python

  nsys.write('POSCAR')
  nsys.write('pmdini')
  nsys.write_dump()

The ``write`` function will parse file name and choose appropriate file format from the name.
