
.. index:: in.pmd
.. _in-pmd:

==================================================
Input file: *in.pmd*
==================================================


``pmd`` starts with reading setting file ``in.pmd`` and atom-configuration files ``pmdini``.
So simulation settings except the atom configuration must be described in ``in.pmd``.

Example of *in.pmd*
==================================================
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
    min_iteration     10
    num_out_energy    100
  
    flag_out_pmd      1
    num_out_pmd       10
    flag_sort         1
  
    force_type        NN
    cutoff_radius     5.8d0
    cutoff_buffer     0.2d0
  
    flag_damping      2
    damping_coeff     0.5d0
    converge_eps      1d-4
    converge_num      3
  
    initial_temperature    -2000d0
    final_temperature      -2000d0
    temperature_control     none
    temperature_target      1  100d0
    temperature_relax_time  1d0
  
    factor_direction 3 2
      1.000d0  1.000d0  1.000d0
      1.000d0  0.000d0  1.000d0
  
    stress_control      none
    stress_relax_time   100d0
    stress_target
      0.00d0   0.00d0   0.00d0
      0.00d0   0.00d0   0.00d0
      0.00d0   0.00d0   0.00d0
    pressure_target     1.00

    mass  1    28.0855
    mass  2     4.0

    boundary   ppp

Here, lines begin with ``!`` or ``#`` are treated as comment lines.


Input parameters
========================================
Parameters users can set in ``in.pmd`` are following:

* :ref:`num_nodes`
* :ref:`io_format`
* :ref:`print_level`
* :ref:`time_interval`
* :ref:`num_iteration`
* :ref:`min_iteration`
* :ref:`num_out_energy`
* :ref:`flag_out_pmd`
* :ref:`num_out_pmd`
* :ref:`flag_sort`
* :ref:`force_type`
* :ref:`flag_damping`
* :ref:`damping_coeff`
* :ref:`converge_eps`
* :ref:`converge_num`
* :ref:`initial_temperature`
* :ref:`final_temperature`
* :ref:`temperature_control`
* :ref:`temperature_target`
* :ref:`temperature_relax_time`
* :ref:`pressure_target`
* :ref:`stress_control`
* :ref:`stress_target`
* :ref:`stress_relax_time`
* :ref:`zload_type`
* :ref:`final_strain`
* :ref:`shear_stress`
* :ref:`cutoff_radius`
* :ref:`flag_temp_dist`
* :ref:`num_temp_dist`
* :ref:`boundary`


------------------------

.. _num_nodes:

num_nodes_{x,y,z}
------------------------------

* Default: 1

Number of division in x, y, or z direction.
The product of these, :math:`xyz`, should be the same as the number of divided atom-configuration files
and computer nodes specified when executing ``mpirun`` or ``mpiexec`` command.


------------------------

.. _io_format:

io_format
------------------------------

* Default: ascii

You can choose either ``ascii`` or ``binary`` format of atom-configuration files.
When you perform large scale simulation, you should choose ``binary`` for efficient reading/writing
atom-configuration files.


------------------------

.. _print_level:

print_level
------------------------------

* Default: 1

How much information is written out during the run.

1:
  Normal information for MD simulation run.

10:
  Info for fitpot is written out such as ``out.NN.gsf``, ``out.NN.dgsf``, ``erg0000`` and ``frc0000``.

100:
  Debug info is written out.

------------------------

.. _time_interval: 

time_interval
------------------------------

* Default: 1.0

Time interval in the unit of **femto second**.


------------------------

.. _num_iteration:

num_iteration / num_steps
------------------------------

* Default: 0

Number of MD steps.
Simulation time equals ``time_interval`` times ``num_iteration``.


------------------------

.. _min_iteration:

min_iteration / min_steps
------------------------------

* Default: 0

Minimum number of MD steps.
In the case you want the MD simulation at least *min_iteration*, you can set this parameter.

------------------------

.. _num_out_energy:

num_out_energy
------------------------------

* Default: 1000

Number of outputs of energies.


------------------------

.. _flag_out_pmd:

flag_out_pmd
------------------------------

* Default: 1

A flag whether or not to write atomic configurations to files at certain steps.

0:
  Not to write.

1:
  Write *pmd*-format atomic configurations to files ``pmd####`` where ``####`` indicates sequential number of the files.

2:
  Write LAMMPS *dump*-fomrat atomic configurations to files ``dump####``.

------------------------

.. _num_out_pmd:

num_out_pmd
------------------------------

* Default: 10

Number of atom-configuration files to be written.


------------------------

.. _flag_sort:

flag_sort
------------------------------

* Default: 1

A flag whether or not to sort the order of atoms by tag before writing out atomic configurations to *pmd* files. It might cost some time for large scale simulation.

1:
  Do sorting

2:
  Do not sorting


------------------------

.. _force_type:

force_type
------------------------------

* Default: ``LJ-Ar``

Choice of the interatomic potential.
Available potentials are listed below:

* ``LJ_Ar`` : Lennard-Jones potential for Ar system.
* ``SW_Si`` : Stillinger-Weber potential for Si system.
* ``EDIP_Si`` : Environment Dependent Interatomic Potential for Si system.
* ``Ito3_WHe`` : EAM potential for W-He system made by Ito et al. @NIFS.
* ``NN`` : Neural-network potential that requires two input files ``in.const.NN`` and ``in.params.NN``.

------------------------

.. _flag_damping:

flag_damping
------------------------------

* Default: 0

A flag whether or not damp atom velocities.

0:
  No damping.

1:
  Simple damped MD using the following **damping_coeff**.

2:
  FIRE algorithm. This is usually much faster and stabler than the simple damped MD.


------------------------

.. _damping_coeff:

damping_coeff
------------------------------

* Default: 0.9

Damping coefficient.


------------------------


.. _converge_eps:

converge_eps
------------------------------

* Default: 1d-4

Convergence criterion in eV.


------------------------


.. _converge_num:

converge_num
------------------------------

* Default: 1

Convergence is achieved if the convergence criterion is sufficed this times successively.


------------------------


.. _initial_temperature:

initial_temperature
------------------------------

* Default: -1.0

Initial temperature of all atoms.


------------------------


.. _final_temperature:

final_temperature
------------------------------

* Default: -1.0

Final temperature of all atoms.
If it is set, target temperature of all atoms changes linearly from ``initial_temperature`` as simulation proceeds.


------------------------


.. _temperature_control:

temperature_control
------------------------------

* Default: none

Temperature-control method, ``none``, ``Berendsen``, and ``Langevin`` are now available.


------------------------

.. _temperature_target:

temperature_target
------------------------------

* Default: 300.0

Target temperature (K) of atoms specified by *ifmv*.
For example,
::

   temperature_target   1  300.0

indicates setting the temperature of *ifmv=1* to 300 K.


------------------------

.. _temperature_relax_time:

temperature_relax_time
------------------------------

* Default: 100.0

Relaxation time of Berendsen thermostat (fs).

-------------------------

.. _stress_control:

stress_control
------------------------------

* Default: none

Type of barostat. Following methods are available:

* ``Berendsen`` / ``vc-Berendsen``: variable-cell Berendsen method.
* ``vv-Berendsen``: variable-volume Berendsen method which keeps the cell shape but changes volume by scaling all the cell vectors.

See *Berendsen, C., Postma, P. M., Van Gunsteren, W. F., Dinola, A., & Haak, R. (1984). Molecular dynamics with coupling to an external bath, 81(8), 3684–3690.* for the detail.

-------------------------

.. _pressure_target:

pressure_target
------------------------------

Default:  0.00

Target hydrostatic pressure which only works when *stress_control* is ``vv-Berendsen``.

-------------------------

.. _stress_target:

stress_target
------------------------------

Default:
  | 0.00d0   0.00d0   0.00d0
  | 0.00d0   0.00d0   0.00d0
  | 0.00d0   0.00d0   0.00d0

Target stress tensor which only works when *stress_control* is ``Berendsen`` or ``vc-Berendsen``.

-------------------------

.. _stress_relax_time:

stress_relax_time
------------------------------

* Default: 100d0 [fs]

Relaxation time of the *Berendsen* barostat.


------------------------

.. _zload_type:

zload_type
------------------------------

* Default:  no

How to apply z-direction strain:

* ``atoms`` :  atoms whose ifmv value are 2 and relative z-position are over 0.5 are moved to upward, those of relative z-position under 0.5 are moved downward.
* ``box`` : control z-component of simulation box matrix.
* ``no`` : Not to apply z-direction strain.


------------------------

.. _final_strain:

final_strain
------------------------------

* Default: 0.0

Final strain value (\%).
Thus strain rate can be given as ``final_strain`` / ( ``time_interval`` * ``num_iteration`` ).


------------------------

.. _shear_stress:

shear_stress
------------------------------

* Default: 0.0

Shear stress value applied to the system.


------------------------

.. _cutoff_radius:

cutoff_radius
------------------------------

* Default: 5.0

Cutoff radius (Angstrom) of the interatomic potential used.


------------------------

.. _flag_temp_dist:

flag_temp_dist
------------------------------

* Default: .false.

Flag about whether or not writing out temperature distribution to ``out.temp-dist`` file.



------------------------

.. _num_temp_dist:

num_temp_dist
------------------------------

* Default: 1

Number of bins along *x*-direction where the temperature is calculated.
This value must be a multiple of ``num_nodes_x``.


------------------------

.. _boundary:

boundary
------------------------------

* Default: ``ppp``

Boundary conditions for each axis, 123.

* ``p``: periodic boundary condition
* ``f``: free boundary condition








