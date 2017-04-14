.. _pmd-file:

===========================
Atom-configuration file
===========================
Original file format described here is used in ``pmd`` .


File format
====================

::

   1:  2.855300000E+000
   2:  3.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   3:  0.00000000000000E+000  3.00000000000000E+000  0.00000000000000E+000
   4:  0.00000000000000E+000  0.00000000000000E+000  3.00000000000000E+000
   5:  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   6:  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   7:  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
   8:        55
   9:  1.10000000000001E+000  1.000E-007  1.000E-007  1.000E-007  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
  10:  1.10000000000002E+000  1.666E-001  1.666E-001  1.666E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
  11:  ...
  12:  1.10000000000054E+000  9.00E-001  9.00E-007  9.00E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05
  13:  2.10000000000055E+000  5.33E-001  5.33E-001  5.33E-001  3.62E-004  1.60E-004  9.60E-004  2.29E-01 -4.12E+00 -4.70E-03 -3.82E-15 -5.34E-15 -5.34E-15 -4.70E-03 -7.47E-05

Here, line numbers are shown for the ease of explanation.

Line 1:
  Superficial or apparent lattie constant. This value is to be multiplied to the cell vectors below to obtain absolute cell vectors.

Line 2 to 4:
  Lattice vectors. The 2nd line is *a1* vector, 3rd line for *a2*, and 4th line for *a3*.
  Columns 1, 2, and 3 are *x* , *y* , *z* components of each vector.

Line 5 to 7:
  Velocities of lattice vectors that are used in *NpT* -ensemble simulation which involves lattice deformation. 

Line 8:
  Number of atoms in the system or decomposed region.

After line 9:
  One atom information per one line.

1st column after line 9:
  Tag of an atom. The digit in the one's place means species of the atom.
  The digit in the tenth's place is *ifmv* value which controls the direction of motion of the atom.

2-4th column after line 9:
  *x* , *y* , and *z* coordinates of the atom normalized by the lattice vectors. Thus they should be in (0:1]

5-7th column after line 9:
  *x* , *y* , and *z* components of the atom vector that are also normalized.

8-15the column after line 9:
  Atomic kinetic energy, atomic potential energy, and atomic stress tensor components, 1,2,3,4,5,6.

----------------

Sample Fortran code
==============================

.. code-block:: fortran

      open(ionum,file=cfname,status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      write(ionum,'(7es23.14e3,8es22.14)') (tag(i),ra(1:3,i)+sorg(1:3)
     &     ,va(1:3,i)/dt
     &     ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i)
     &     ,epi(i),
     &     ,strs(1,1,i)*up2gpa,strs(2,2,i)*up2gpa,strs(3,3,i)*up2gpa,
     &     ,strs(2,3,i)*up2gpa,strs(1,3,i)*up2gpa,strs(1,2,i)*up2gpa,i=1,natm)
      close(ionum)

Detail explanations of variables are omitted.
Users can write their own code by following this sample Fortran code.


------------

.. _format_conversion:

Format conversion
==================

There is a python utility, ``pmdsys.py``,  which converts this *pmd* format to *POSCAR* or *akr* format.

You can use like following,
::

  $ python /path/to/nap/nappy/pmdsys.py convert pmdini POSCAR

Here ``pmdini`` file will be converted to ``POSCAR`` file with *POSCAR* format.

See the help message for more details as,
::

  $ python /path/to/nap/nappy/pmdsys.py -h


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


