==============================
nappy | NAP PYthon utilities
==============================

``nappy`` is set of python utilities. 
It is contained in ``nappy`` directory, but currenly it is not exactly a python module.
Users may have to use those utilities by calling them directly from the shell.


Read and write files
==============================

To read atomic structures from files of several formats,

.. code:: python

   atoms = NAPSystem(fname='POSCAR')
   atoms = NAPSystem(fname='structure.dat', ffmt='xsf', specorder=['Al','O'])


Current available formats are:

* ``pmd``: file format specific for ``pmd`` program
* ``POSCAR``: VASP POSCAR file
* ``dump``: LAMMPS dumps file
* ``xsf``
* ``akr``: file format for Akira viewer

And to write the structure to a file,

.. code:: python

  atoms.write('POSCAR')
  atoms.write('pmdini')
  atoms.write_dump()

The ``write`` function will parse file name and choose appropriate file format from the name.
