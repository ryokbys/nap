# fitpot -- fit parameters of neural-network potential

The validity of MD simulations depends strongly on the accuracy of the
interatomic potential used in the simulation. So, when you think of
doing some simulation of specific system, you have to prepare an
interatomic potential that reproduces the phenomena you are thinking of.

Here we indroduce how to make a neural-network potential and fit
potential parameters with the *fitpot* program included in **nap**
package.

!!! Note
    Currently, the *fitpot* program is used only for neural-network
    potential. For other classical potentials, use *fp.py* instead.


## What does fitpot do?

In the *fitpot*, the following loss function is minimized by optimizing
potential parameters $\{ w \}$.

$$\mathcal{L}(\{w\}) = \frac{1}{\eta N_s}\sum_s^{N_s} \left[ \Delta E^2 +\sum_i^{N^s_\mathrm{a}}\left| \Delta \boldsymbol{F}_i\right|^2 +\left| \Delta \sigma \right|^2\right]$$

in the case of fitting energies and forces. Here, $s$ is the sample
number, $N_s$ the number of samples, $N^s_\mathrm{a}$ the number of
atoms in the sample $s$, $\eta$ the parameter corresponding to how many
properties are fitted (in this case $\eta = 3$ because energy, force and
stress are used.)

To minimize the above loss function, the following gradient-based
methods are available in *fitpot*:

-   Steepest descent (SD)
-   Quasi-Newton method (BFGS)
-   Conjugate gradient (CG)
-   Stochastic gradient descent (SGD)

## Compilation

Since some modules in *pmd* program are required for the compilation of
*fitpot*, compile *pmd* before compiling *fitpot*. :

    $ cd /path/to/nap/
    $ ./configure --prefix=$(pwd)
    $ cd pmd
    $ make pmd
    $ cd ../fitpot/
    $ make fitpot

------------------------------------------------------------------------

## Quick trial with an example

There is an example of *fitpot* with minimal dataset to see how it
works. Go to the directory `examples/fitpot_DNN_SiO/`, read `README.md`,
try running *fitpot*, and look at some output files.

------------------------------------------------------------------------

## Fitting procedure

Hereafter, we assume that the reference data are obtained by using an
*ab-initio* calculation program, VASP.

Potential parameters are fitted as the following procedure:

1. [vasp-data](#prepare-reference-data)
2. [prepare-inputs](#prepare-input-files)
3. [exec-fitpot](#run-fitpot-program)

------------------------------------------------------------------------

### Prepare reference data

Assuming that there are some reference data in `dataset/` directory, and
all the data are stored in the directories whose names start with
`smpl_`.

The following files are required in eacy sample directory (`smpl_*`):

-   `pos`
-   `erg.ref`
-   `frc.ref`
-   `strs.ref`

`pos` is a pmd-format atom-configuration file, `erg.erg` contains a
scalar value of energy of the system, `frc.ref` contains the number of
atoms in the system and forces of all atoms shown as, :

    4
    0.1000   0.0000   0.0000
    0.0000   0.1000   0.0000
    0.0000  -0.1000   0.0000
    -0.1000   0.0000   0.0000

In the case of extracting DFT data from *ab-initio* MD runs with VASP,
positions, energy, forces and stress of each MD step can be obtained
from `vasprun.xml` file as follows, :

    $ python path/to/nap/nappy/vasp/vasprun2fp.py /path/to/dir/that/includes/vasprun.xml/

Then you get directories with names like `#####` including `pos`,
`erg.ref`, `frc.ref` and `strs.ref` files in them.

### Prepare input files

Inputs files needed for *fitpot* are the following:

-   `in.fitpot`
-   `in.params.DNN`
-   `in.params.desc`
-   `in.params.Coulomb` in each `smpl_XXX` directory in some special
     cases

You have to specify the `num_samples` in `in.fitpot` file which is the
number of samples in `dataset/` directory. The number of sample
directories can be counted by the following command,

```bash
$ ls /path/to/dataset | grep smpl_ -c
```

### Run fitpot program

In the directory where `dataset/` directory and `in.fitpot` file exist,
you can run the *fitpot* program as, :

    $ ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

Or if you want it to run in parallel mode, :

    $ mpirun -np 10 ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

There are some output files:

- `out.erg.trn.fin`, `out.erg.tst.fin` -- These files include reference and *pmd* data of energies. To see whether the fitting went well or not, plot these data by using `gnuplot` as,

        $ gnuplot
        gnuplot> plot 'out.erg.trn.fin' us 1:2 w p t 'training set'
        gnuplot> rep 'out.erg.tst.fin' us 1:2 w p t 'test set'

- `out.frc.trn.fin`, `out.frc.tst.fin` --  These files include reference and *pmd* data of forces.

------------------------------------------------------------------------

Input file for *fitpot*
-----------------------

The following code shows an example of the input file `in.fitpot`.

    num_samples       14
    test_ratio        0.1
    num_iteration     100
    num_iter_eval     1

    fitting_method    bfgs
    sample_directory  "../dataset"
    param_file        in.params.DNN
    normalize_input   none

    energy_match       T
    force_match        T
    stress_match       T
    potential          DNN

    ftol              1.0e-5
    xtol              1.0e-4

    penalty           none
    penalty_weight    1d-3

    # Species order:  1) Al, 2) Mg, 3) Si
    specorder    Al  Mg  Si


------------------------------------------------------------------------

### num_samples

Default: *none*

Number of reference samples to be used for training and test.

------------------------------------------------------------------------

### sample_list

Default: `none`

Path to the file that contains a list of samples to be used for training
and test. The format of the list file should be like, :

    smpl_001
    smpl_002
    smpl_003
    ...

or with specifying which samples are training (`1`) or test (`2`) as, :

    smpl_001  1
    smpl_002  2
    smpl_003  1
    ...

If whether training or test is specified in the list,
[test_ratio](#test_ratio) will be neglected.

------------------------------------------------------------------------

### test_ratio

Default: `0.1`

The ratio of test data set $r$ within whole data set $N$. Thus the
number of test data set is $rN$, and the number of training data set is
$(1-r)N$.

------------------------------------------------------------------------

### num_iteration

Default: `1`

Number of iterations of a minimization method.

------------------------------------------------------------------------

### num_iter_eval

Default: `1`

Test data set will be evaluated every *num_iter_eval* iterations.


------------------------------------------------------------------------

### fitting_method

Default: *test*

The method used to fit parameters to the sample data. Available methods
are the following:

- `sd`/`SD` -- Steepest descent algorithm which requires gradient information.
- `cg`/`CG` -- Conjugate gradient algorithm which requires gradient information.
- `bfgs`/`BFGS` -- Quasi-Newton method with BFGS. This requires gradient information.
- `check_grad` -- Comparison of analytical derivative and numerical derivative. Use this to check the implemented analytical gradient.
- `test`/`TEST` -- Just calculate function L and gradient of L w.r.t. fitting parameters.


------------------------------------------------------------------------

### sample_directory

Default: `dataset`

The directory that includes sample data. We call this `dataset` in the
above instruction.

If you want to use `..` to specify the directory relative to the current
working directory, e.g. `../dataset`, you need to enclose with
double-quotation marks like `"../dataset"`.

------------------------------------------------------------------------

### param_file

Default: *in.params.DNN*

The name of the file that has parameter values in it. This is passed to
`pmd` program.

------------------------------------------------------------------------

### ftol

Default: *1.0e-5*

The tolerance of difference of the loss function value.

------------------------------------------------------------------------

### xtol

Default: *1.0e-4*

The tolerance of the change of variables which are optimized. If either
one of [ftol]{.title-ref} or [xtol]{.title-ref} is achieved, the
optimization stops.

------------------------------------------------------------------------

### energy_match, force_match, stress_match

Default: *True* for energy, *False* for force and stress

Whether or not to match forces. ( *True* or *False* ) It is recommended
to match not only energy but also forces, since forces are important for
molecular dynamics.

------------------------------------------------------------------------

### potential or force_field

Default: `DNN`

The potential whose parameters you are going to fit.
Potentials currently available are:

- `DNN` -- Neural-network potential
- `linreg` -- Linear regression potential

------------------------------------------------------------------------

### random_seed

Default: *12345d0*

Initial random seed for the uniform random numbers used in the *fitpot*.
This is used to change the random choice of training and test sets.

------------------------------------------------------------------------

### regularize

Whether or not regularize bases obtained in *linreg* and *DNN*
potentials. ( *True* or *False* )

Default: *False*

------------------------------------------------------------------------

### penalty

Default: *no*

Type of penalty term, *lasso* which is L1-norm penalty or *ridge* which
is L2-norm penalty, or *no* which means no penalty term.


------------------------------------------------------------------------

### penalty_weight

Default: *1.0*

The weight applied to the penalty term. This value also has to be
determined through cross-validation scoring\...

------------------------------------------------------------------------

### sample_error

Default: *0*

The number of samples whose errors are to be given. These errors appear
at the denominators of energy and force in the evaluation function such
that

$$\left( \frac{E^\mathrm{NN}-E^\mathrm{DFT}}{N_\mathrm{a}\varepsilon_\mathrm{e}}\right)^2 +\sum_i^{N_\mathrm{a}} \sum_\alpha^{xyz} \frac{1}{3N_\mathrm{a}}\left( \frac{F^\mathrm{NN}_{i\alpha} -F^\mathrm{DFT}_{i\alpha}}{\varepsilon_\mathrm{f}}\right)^2$$

If the difference between NN energy and DFT energy/force is lower than
this value, this term becomes less than 1.0, which means the
energy/force of the sample is thought to be converged. The initial
values of the errors are 0.001 (eV/atom) and 0.1 (eV/Ang) for energy and
force, respectively.

There must be the same number of following entry lines as the above
value which determine the errors of energy and force of each sample like
the this, :

    sample_error   2
        Al_fcc    0.001  0.2  1.0
        Al_bcc    0.001  0.2  1.0

The each entry has *entry_name*, *error of energy (eV/atom)*, *error of
forces (eV/Ang)* and *error of stresses (GPa)*. The error values are
applied to all the samples that contain *entry_name* in their directory
names.

------------------------------------------------------------------------

### force_denom_type

`relative` or `absolute`

Default: `relative`

Which type of denominator of force term in the loss function is used. If
`absolute` is specified, the *fitpot* uses an *error of forces*
specified in the [sample_error](#sample_error) for the
denominator of force term. If `relative` is specified, the *fitpot* uses
a magnitude of force on the atom in the denominator of force term.

------------------------------------------------------------------------

### specorder

Default: *none*

The order of species common in fitpot. This must be specified before
`atom_energy` entry and must hold for every samples.

------------------------------------------------------------------------

### init_params

Default: `read`

Whether the paramters to be optimized are read from the file or
initialized.

- `read` -- Read parameters from the file.
- `gaussian` --  Parameters are initialized with Gaussian distribution according to *init_params_sgm* and *init_params_mu*.

------------------------------------------------------------------------

### init_params_sgm

Default: `1d0`

Variance of Gaussian distribution of the initial values for parameters.

------------------------------------------------------------------------

### init_params_mu

Default: `0d0`

Mean value of Gaussian distribution of the initial values for
parameters.

------------------------------------------------------------------------

### init_params_rs

Default: `12345.0`

Random seed for the initialization of parameters. This random seed is
only used for this purpose and does not affect random seed for the
choice of training and test sets, which is affected by
`random_seed`{.interpreted-text role="ref"}.

