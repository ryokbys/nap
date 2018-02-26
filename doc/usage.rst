=================
Usage
=================

* :ref:`install`
* :ref:`sample_simulation`
* :ref:`input_files`
* :ref:`units`
* :ref:`make_init_config_file`
* :ref:`make_in_pmd`
* :ref:`run_pmd`
* :ref:`massive_parallel`

-----------



.. index:: install
.. _install:

Install
=====================


Download
--------------

Download the source code of current version of the whole packange from the following Github site:

  https://github.com/ryokbys/nap


Compile
------------
You can download ``nap-master.zip`` file from the site.
And you can get ``nap-master`` directory when you unzip the zip file.
For the ease of following explanation, change the directory name to ``nap`` as,
::

  $ unzip nap-master.zip
  $ mv nap-master nap
  $ cd nap/

Then you can compile the ``pmd`` program as following,
::

  $ ./configure --prefix=$(pwd)
  $ cd pmd
  $ make pmd

If you get an error of finding an MPI Fortran compiler when you are doing ``configure``,
you have to find an MPI Fortran compiler by yourself or asking the system administrator and
do ``configure`` again by specifying the compiler position as,
::

  $ ./configure --prefix=$(pwd) FC=<path/to/fortran-compiler>

The option ``--prefx=$(pwd)`` is not necessary for ``pmd``, but when you use ``fitpot`` program,
it is required.

.. note::

  If you get an error related to C preprocessor such as 
  ::
   
    configure: error: C preprocessor "/lib/cpp" fails sanity check

  you may have to specify true C preprocessor path to ``configure`` command by adding an option like ``CPP=/usr/bin/cpp``.

However the ``pmd`` command should be available with this compilation,
this ``pmd`` may not be so efficient in the system you are compiling.
So you might need to add some options relevant for your system.

gfortran in Mac OSX
_____________________
  In the case of ``gfortran`` with ``openmpi`` ,
  ::
  
    $ ./configure --prefix=$(pwd) FCFLAGS="-O4 -g"
  
  The optimization option ``-O4`` is thought to be the best for this system.


Intel Fortran Compiler
______________________
  If you can use ``ifort`` in your system, for example,
  ::
  
    $ ./configure --prefix=$(pwd) FCFLAGS="-xHOST -O3 -ip -ipo -no-prec-div"
  
  The options ``-ip`` and ``-ipo`` have to do with inline expansions and 
  are relevant to the efficiency of ``pmd``.
  

Fujitsu Fortran in FX?
_____________________________

  It is easier to compile on the computation node not on the login node. Since there are some difference about configuring/compiling on those nodes.
  To configure and compile the ``pmd``, first you need to login to a computation node by doing ``pjsub --interact``.
  ::
    
     $ pjsub --interact
     $ ./configure --prefix=$(pwd) FCFLAGS="-O3"
     $ cd pmd
     $ make pmd
     $ exit
    
  
  
  .. note::
  
     In case that the ``configure`` returns errors and exit without completing the configuration
     and the error message is related to cross compilation,
     you may need to add an option like ``--host=sparc64`` to the above command line.

Fujitsu Fortran in CX400
______________________________
  In the case of Fujitsu Fortran compiler ``mpifrt`` in CX400,
  ::
  
    $ ./configure --prefix=$(pwd) FCFLAGS="-Kfast,parallel"
  

Helios in Rokkasho-mura
_______________________________
  It is Linux OS on Intel CPU, and the compilation seems to be basic one.
  But one needs to add specific options as following,
  ::
  
    $ ./configure --prefix=$(pwd) FC=mpiifort FCFLAGS="-xAVX -O3 -ip -ipo -g -CB"
  
  If you don't specify the ``mpiifor`` explicitly, ``ifort`` is set by default and the compilation does not work correctly.

----------

.. _sample_simulation:

Run a sample simulation
=====================================
There are some input files in ``example/test-W/`` directory. ( ``in.pmd`` and ``pmdini`` )
These input files are for the system of BCC tungsten crystalline structure including one helium atom.

::

    $ cd example/test-W
    $ ../../pmd/pmd

When you run the ``pmd`` command like above, *NVE* -MD simulation of 100 steps is performed.
And the total, kinetic, and potential energies are output in ``out.erg`` file.
So you can look at the evoluation of these energies using ``gnuplot`` command as,
::

    $ gnuplot
    gnuplot> plot 'out.erg' us 1:3 w l, 'out.erg' us 1:4 w l, 'out.erg' us 1:5 w l

In this case, since you are performing *NVE* -MD simulation at 300K of diamond-Si,
the total energy conserves conpensating the deviations of kinetic and potential energies.

.. image:: ./figs/graph_energy-steps.png

And also configurations of atoms at each 10 steps out of 100 steps are written in files, e.g.)
``pmd_0``, ``pmd_10``,..., ``pmd_100``.

.. note::
   The file name format ``pmd####`` has been changed to ``pmd_###`` because this format works for
   larger number of files than previous definition.
   And now the number ``###`` means exactly the MD step in the MD simulation run.

------------

.. _input_files:

Input files needed to run pmd
==================================

``pmd`` must be executed in the directory where these files exist.

in.pmd
  Input file that describes simulation setting.

pmdini
  Cell information and initial positions and velocities of atoms.

.. image:: ./figs/pmd.png

After running ``pmd`` , some output files appear in the same directory.

----------

.. index:: units
.. _units:

Units used in pmd
==========================
* Length: Angstrom
* Time: fs
* Energy: eV
* Mass:  1/12 of carbon atom

---------

.. _make_init_config_file:

Make an initial atom-configuration file
========================================
Please refer :ref:`pmd-file` for details of atom-configuration file.

One has to make an initial atom-configuration file, ``pmdini``, to run ``pmd``.
There are already some programs that make initial atom-configuration files
of some systems ( ``mkconf/mkconf_????.F`` ).
So you can make your own initial atom-configuration file by looking at those program codes.

If there is already a program that makes an atom-configuration file of your target system,
you can make an atom-configuration file as,
::

  $ cd mkconf
  $ emacs makefile
  (find which mkconf_* will be made)
  $ make mkconf_Si_disl
  $ ./mkconf_Si_disl

Here, *Si_disl* is chosen.

Then you get an atom-configuration file ``pmdini``.

.. note::

   If you have to make the program which makes an atom-configuration file,
   copy any program like ``mkconf_BCC.F`` , modify it, add an entry into ``makefile`` , 
   and compile.


---------

.. _make_in_pmd:

Make the *in.pmd* file
============================================
Please refer :ref:`in-pmd` for details of ``in.pmd`` file.

For instance, ``in.pmd`` file for the system of 1000 step MD simulation using ``SW_Si`` potential
is as follows,
::

  #
  #  unit of time  = femto sec
  #  unit of length= Angstrom
  #  unit of mass  = unified atomic mass unit
  #
  num_nodes_x       1
  num_nodes_y       1
  num_nodes_z       1
  
  io_format         ascii
  print_level       1
  
  time_interval     2d0
  num_iteration     1000
  num_out_energy    100
  
  flag_out_pmd      1
  num_out_pmd       10
  
  force_type        SW_Si
  cutoff_radius     3.7712d0
  cutoff_buffer     0.2d0
  
  flag_damping      2
  damping_coeff     0.5d0
  converge_eps      1d-4
  converge_num      3
  
  initial_temperature     -2000d0
  final_temperature     -2000d0
  temperature_control     none
  temperature_target      100d0
  temperature_relax_time  1d0
  
  factor_direction 3 2
    1.000d0  1.000d0  1.000d0
    1.000d0  0.000d0  1.000d0
  
  stress_control       none
  stress_relax_time   100d0
  stress_target
    0.00d0   0.00d0   0.00d0
    0.00d0   0.00d0   0.00d0
    0.00d0   0.00d0   0.00d0
  pressure_target     1.00
  
  shear_stress   0.00
  
  mass  1    28.0855
  mass  2     4.0

Here, the lines begins with ``!`` or ``#`` are treated as comment lines
and blanc lines are skipped.

------------

.. _run_pmd:

Run pmd
================

Run pmd on 1-process
--------------------------------
It is really easy to run ``pmd`` on 1-process.
On the directory where ``in.pmd`` and ``pmdini`` exist, just execute ``pmd`` as,
::

  $ /path/to/pmd/pmd

If you want to perform it background, 
::

  $ /path/to/pmd/pmd > out.pmd 2>&1 &


The following files appear when you perform ``pmd`` :

  ``out.erg``
        Total, kinetic, potential energies, and temperature, volume, pressure.
  ``pmd_##``
        Atom-configurations at a certain MD step. ``##`` means the MD step.



Run pmd on parallel-nodes
------------------------------------
Different from the old version of **pmd** which requires divided atom configuration files for parallel nodes, in the current version (since 05May2016), the parallel simulation can be performed almost the same as the serial run.

Just you need to describe how many divisions on each direction in ``in.pmd`` such as ``num_nodes_x``, ``num_nodes_y`` and ``num_nodes_z`` ,and run **pmd** with ``mpirun`` or ``mpiexec`` command to run MPI executable.
::

  $ mpirun -np 8 --machinefile hosts.list /path/to/pmd > out.pmd 2>&1 &

Here, ``pmd`` is executed on 8-nodes listed in ``hosts.list`` and 
standard output is written in ``out.pmd`` .

If any job-scheduling system is available on the system you are using,
describe the above command in your job script to be submitted.


-------------

.. _massive_parallel:

Notes on performing massively parallel simulation
======================================================
When you perform parallel simulation with over one million atoms,
the data of atom-configuration files becomes considerably large
and reading/writing data takes long time compared with intrinsic computation.
So ``pmd`` can read/write binary version of atom-configuration files 
that are way more small amount of data.
If you want to read/write binary files, describe following in ``in.pmd``,
::

  io_format   binary

And also you have to write code of writing binary atom-configuration file in ``mkconf_*.F``.
In ``mkconf_*.F`` files, there is a line of subroutine call ``call write_pmd0_ascii`` , 
you have replace it to ``call write_pmd0_bin`` and recompile the program.

