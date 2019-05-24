.. _pmd-file:

===========================
Atom-configuration file
===========================
Original file format described here is used in ``pmd`` .


File format
====================

::

   1 :!
   2 :!  specorder:  W  H
   3 :!    
   4 :  2.855300000E+000
   5 :  3.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   6 :  0.00000000000000E+000  3.00000000000000E+000  0.00000000000000E+000
   7 :  0.00000000000000E+000  0.00000000000000E+000  3.00000000000000E+000
   8 :  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   9 :  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   10:  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   11:        55
   12:  1.10000000000001E+000  1.000E-007  1.000E-007  1.000E-007  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
   13:  1.10000000000002E+000  1.666E-001  1.666E-001  1.666E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
   14:  ...
   15:  1.10000000000054E+000  9.00E-001  9.00E-007  9.00E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
   16:  2.10000000000055E+000  5.33E-001  5.33E-001  5.33E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05

Here, line numbers are shown for the ease of explanation.

Line 1-3:
  Lines begin with ``!`` are treated as comment lines.
  There are some keywords that are used to specify some additional data to ``pmd`` when they are at a comment line at the beginning.

  - ``specorder:`` specifies the species order used in ``pmd``.

.. note::

   **specorder** must be specified in the current ``pmd`` (since *rev190515*), as the masses and the interatomic potentials are determined using this information.

Line 4:
  Superficial or apparent lattie constant. This value is to be multiplied to the cell vectors below to obtain absolute cell vectors.

Line 5 to 7:
  Lattice vectors. The 2nd line is *a1* vector, 3rd line for *a2*, and 4th line for *a3*.
  Columns 1, 2, and 3 are *x* , *y* , *z* components of each vector.

Line 8 to 10:
  Velocities of lattice vectors that are used in *NpT* -ensemble simulation which involves lattice deformation. 

Line 11:
  Number of atoms in the system or decomposed region.

After line 11:
  One atom information per one line.

1st column after line 11:
  Tag of an atom. The digit in the one's place means species of the atom.
  The digit in the tenth's place is *ifmv* value which controls the direction of motion of the atom.

2-4th column after line 11:
  *x* , *y* , and *z* coordinates of the atom normalized by the lattice vectors. Thus they should be in (0:1]

5-7th column after line 11:
  *x* , *y* , and *z* components of the atom vector that are also normalized.

8-15the column after line 11:
  Atomic kinetic energy, atomic potential energy, and atomic stress tensor components, 1,2,3,4,5,6.

----------------

Sample Fortran code
==============================

.. code-block:: fortran

   open(ionum,file=cfname,status='replace')
   write(ionum,'(es23.14e3)') hunit
   write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
        ,ib=1,3),l=0,1)
   write(ionum,'(i10)') natm
   write(ionum,'(7es23.14e3,8es22.14)') (tag(i),ra(1:3,i)+sorg(1:3) &
        ,va(1:3,i)/dt &
        ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i) &
        ,epi(i), &
        ,strs(1,1,i)*up2gpa,strs(2,2,i)*up2gpa,strs(3,3,i)*up2gpa, &
        ,strs(2,3,i)*up2gpa,strs(1,3,i)*up2gpa,strs(1,2,i)*up2gpa,i=1,natm)
   close(ionum)

Detail explanations of variables are omitted.
Users can write their own code by following this sample Fortran code.


------------

.. _format_conversion:

Format conversion
===================

There is a python utility, ``napsys.py`` that can convert files amoung the following formats,

  - ``pmd``: input format for **pmd** program.
  - ``POSCAR``: input format of **VASP** program.
  - ``dump``: output format of **LAMMPS** dump command.

You can use like following,
::

  $ python /path/to/nap/nappy/napsys.py convert pmdini POSCAR

Here ``pmdini`` file will be converted to ``POSCAR`` file in POSCAR format, where the file format is determined automatically from file names. Users can also specify input/output file format by the options ``--in-format`` and ``--out-format``.

See the help message for more details as,
::

  $ python /path/to/nap/nappy/napsys.py -h


---------------

.. _cell_maker:

Make crystalline structures
==============================

There is also a python utility, ``cell_maker.py``, which makes typical conventional crystalline structures.
You can make a *pmd* format file of diamond structured cubic system with 8 atoms as,
::

  $ python /path/to/nap/nappy/mkcell/cell_maker.py diamond -l 5.473 -o pmdini

The option ``-l`` specifies the lattice constant of the lattice.
Output format is automatically detected from the file name.
You can also make *fcc*, *bcc*, *sc (simple cubic)*, and *hcp* structures as well.


