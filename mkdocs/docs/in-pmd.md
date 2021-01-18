# Input file: in.pmd

*pmd* starts with reading setting file `in.pmd` and atom-configuration
files `pmdini`. So simulation settings except the atom configuration
must be described in `in.pmd`.

## Example of in.pmd


    #
    #  unit of time  = femto sec
    #  unit of length= Angstrom
    #  unit of mass  = unified atomic mass unit
    #
      io_format         ascii
      print_level       1

      time_interval     2d0
      num_iteration     1000
      min_iteration     10
      num_out_energy    100

      flag_out_pos      1
      num_out_pos       10

      force_type        Morse Coulomb
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

      boundary   ppp

Here, lines begin with `!` or `#` are treated as comment lines.


## Input parameters


### num_nodes_x

- Default: `-1`

Number of division in x, y, or z direction. If one of these is
non-positive (`<=0`), these numbers are automatically estimated from the
system size and the number of MPI processes used. If all of these are
positive, specified values are used. The product of these, $xyz$, should
be the same as the number of divided atom-configuration files and
computer nodes specified when executing *mpirun* or *mpiexec* command.

------------------------------------------------------------------------

### io_format

- Default: `ascii`

You can choose either `ascii` or `binary` format of atom-configuration
files. When you perform large scale simulation, you should choose
`binary` for efficient reading/writing atom-configuration files.

------------------------------------------------------------------------

### print_level

- Default: `1`

How much information is written out during the run.

- `1` -- Normal information for MD simulation run.
- `100` -- Debug info is written out.

------------------------------------------------------------------------

### time_interval

-   Default: 1.0

Time interval in the unit of **femto second**. If negative, it activates
*variable time-step mode* and its absolute value is the maximum time
interval $\Delta t_{\max}$ in the mode. Appropriate range of dt_max
would be 2.0 to 5.0 depending on the minimum mass of ion in the system.

------------------------------------------------------------------------

### vardt_length_scale

-   Default: 0.1

The specific length $L^*$ of the *variable time-step mode* where the
time interval is determined as

$$\begin{equation}
   \Delta t = \min \left( \Delta t_\mathrm{max}, \frac{L^*}{v_\mathrm{max}}\right).
\end{equation}$$

------------------------------------------------------------------------

### num_iteration

- Default: 0
- Alternative: `num_steps`

Number of MD steps. Simulation time equals `time_interval` times
`num_iteration`.

------------------------------------------------------------------------

### min_iteration

- Default: 0
- Alternative: `min_steps`

Minimum number of MD steps. In the case you want the MD simulation at
least *min_iteration*, you can set this parameter.

------------------------------------------------------------------------

### num_out_energy

- Default: 1000

Number of outputs of energies.

------------------------------------------------------------------------

### flag_out_pos

- Default: 1
- Alternative: `flag_out_pmd`

A flag whether or not to write atomic configurations to files at certain
steps.

- `0` -- Not to write.
- `1` -- Write *pmd*-format atomic configurations to files `pmd_####` where `####` indicates sequential number of the files.
- `2` -- Write LAMMPS *dump*-fomrat atomic configurations to files `dump_####`.

------------------------------------------------------------------------

### num_out_pos

- Default: 10
- Alternative: `num_out_pmd`

Number of atom-configuration files to be written.

------------------------------------------------------------------------

### flag_sort

- Default: `1`

A flag whether or not to sort the order of atoms by tag before writing
out atomic configurations to *pmd* files. It might cost some time for
large scale simulation.

- `1` -- Do sorting
- `2` -- Do not sorting

------------------------------------------------------------------------

### force_type

- Default: `None`

Choice of the interatomic potential. Available potentials are listed
below:

- `LJ` : Lennard-Jones potential for Ar system.
- `SW_Si` : Stillinger-Weber potential for Si system.
- `EDIP_Si` : Environment Dependent Interatomic Potential for Si system.
- `Ito3_WHe` : EAM potential for W-He system made by Ito et al. at NIFS.
- `Morse` : Morse potential that requires an input files `in.params.Morse`.
- `Coulomb` : Coulomb potential that requires an input files `in.params.Coulomb`.
- `DNN` : Neural-network potential that requires two input files `in.params.desc` and `in.params.NN2`.

------------------------------------------------------------------------

### flag_damping

- Default: `0`

A flag whether or not damp atom velocities.

- `0` -- No damping.
- `1` -- Simple damped MD using the following **damping_coeff**.
- `2` -- FIRE algorithm. This is usually much faster and stabler than the simple damped MD.

------------------------------------------------------------------------

### damping_coeff

- Default: `0.9`

Damping coefficient.

------------------------------------------------------------------------

### converge_eps

- Default: `1d-4`

Convergence criterion in eV. If it is negative value, not to stop because of convergence.

------------------------------------------------------------------------

### converge_num

- Default: `1`

Convergence is achieved if the convergence criterion is sufficed this times successively.

------------------------------------------------------------------------

### initial_temperature

- Default: `-1.0`

Initial temperature of all atoms.

------------------------------------------------------------------------

### final_temperature

- Default: `-1.0`

Final temperature of all atoms. If it is set, target temperature of all
atoms changes linearly from `initial_temperature` as simulation
proceeds.

------------------------------------------------------------------------

### temperature_control

- Default: `none`

Temperature-control method, `none`, `Berendsen`, and `Langevin` are now
available.

------------------------------------------------------------------------

### temperature_target

- Default: `300.0`

Target temperature (K) of atoms specified by *ifmv*. For example, :

    temperature_target   1  300.0

indicates setting the temperature of *ifmv=1* to 300 K.

------------------------------------------------------------------------

### temperature_relax_time

- Default: `100.0`

Relaxation time of Berendsen thermostat (fs).

------------------------------------------------------------------------

### stress_control

- Default: `none`

Type of barostat. Following methods are available:

- `Berendsen` / `vc-Berendsen`: variable-cell Berendsen method.
- `vv-Berendsen`: variable-volume Berendsen method which keeps the
    cell shape but changes volume by scaling all the cell vectors.

See Berendsen's paper[^Berendsen1984] for the detail.

[^Berendsen1984]: Berendsen, C., Postma, P. M., Van Gunsteren, W. F., Dinola, A., & Haak, R. (1984). Molecular dynamics with coupling to an external bath,
81(8), 3684--3690.

------------------------------------------------------------------------

### pressure_target

Default: `0.0`

Target hydrostatic pressure which only works when *stress_control* is
`vv-Berendsen`.

------------------------------------------------------------------------

### stress_target

Default:

     0.00d0  0.00d0  0.00d0
     0.00d0  0.00d0  0.00d0
     0.00d0  0.00d0  0.00d0

Target stress tensor which only works when *stress_control* is
`Berendsen` or `vc-Berendsen`.

------------------------------------------------------------------------

### stress_relax_time

- Default: `100d0`

Relaxation time (fs) of the *Berendsen* barostat.

------------------------------------------------------------------------

### zload_type

- Default: `no`

How to apply z-direction strain:

- `atoms` -- atoms whose ifmv value are 2 and relative z-position are
    over 0.5 are moved to upward, those of relative z-position under 0.5
    are moved downward.
- `box` -- control z-component of simulation box matrix.
- `no` -- Not to apply z-direction strain.

------------------------------------------------------------------------

### final_strain

- Default: `0.0`

Final strain value (%). Thus strain rate can be given as `final_strain`
/ ( `time_interval` * `num_iteration` ).

------------------------------------------------------------------------

### shear_stress

- Default: `0.0`

Shear stress value applied to the system.

------------------------------------------------------------------------

### cutoff_radius

- Default: `5.0`

Cutoff radius (Angstrom) of the interatomic potential used.

------------------------------------------------------------------------

### flag_temp_dist

- Default: `.false.`

Flag about whether or not writing out temperature distribution to
`out.temp-dist` file.

------------------------------------------------------------------------

### num_temp_dist

- Default: `1`

Number of bins along *x*-direction where the temperature is calculated.
This value must be a multiple of `num_nodes_x`.

------------------------------------------------------------------------

### mass

- Default: masses of given species are set automatically.

If masses of some species need to be set different from those of
elements, masses should be specified as follows. :

    mass   Si  28.0855
    mass   He   4.00

------------------------------------------------------------------------

### boundary

- Default: `ppp`

Boundary conditions for each axis, 123.

- `p`: periodic boundary condition
- `f`: free boundary condition
