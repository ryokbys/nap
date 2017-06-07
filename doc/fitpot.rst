.. Manual for potential-parameter fitting program, fitpot

.. index:: fitpot

==================================================
Potential-parameter fitting
==================================================
The validity of MD simulations depends mostly on the accuracy of the interatomic potential used in the simulation.
So, when you think of doing some simulation of specific system, 
you have to prepare an interatomic potential that reproduces the phenomena you are thinking of.

Here we indroduce how to make interatomic potentials and fit potential parameters
by using ``fitpot`` program included in **NAP** package.

.. note::

   ``fitpot`` program now is specific to **Neural-network (NN)** potential and cannot produce other potential parameters.



Compilation
===============
Some shell and python scripts are used, too. 
To make the shell script be available on your system,
you need to run ``configure`` and to change the permission of ``run_pmd.sh`` script as follows,
::

  $ ./configure --prefix=$(pwd)
  $ cd pmd
  $ make pmd
  $ cd ..
  $ make fitpot
  $ chmod 755 fitpot/run_pmd.sh



Fitting procedure
=========================
Hereafter, we assume that the reference data are obtained by using VASP.

Potential parameters are fitted as the following procedure:

  #. :ref:`vasp-data`
  #. :ref:`prepare-pmd`
  #. :ref:`prepare-scripts`
  #. :ref:`prepare-inputs`
  #. :ref:`exec-fitpot`

-----------------

.. _vasp-data:

Prepare reference data
------------------------------
Assuming that there are some reference data sets in ``dataset/`` directory,
and all the data are stored in the directories whose names start with ``smpl_``.

Files needed to perform *fitpot* in eacy sample directory (``smpl_*``) are:

* ``pos``
* ``erg.ref``
* ``frc.ref``
* ``pmd/out.NN.gsf``
* ``pmd/out.NN.dgsf``

``pos`` is pmd format atom-configuration file, ``erg.erg`` contains only one number of energy of the system,
``frc.ref`` contains the number of atoms in the system and forces of all atoms shown as,
::

   4
    0.1000   0.0000   0.0000
    0.0000   0.1000   0.0000
    0.0000  -0.1000   0.0000
   -0.1000   0.0000   0.0000

``out.NN.gsf`` and ``out.NN.dgsf`` are specific to neural-network (NN) potential and written out by *pmd* program.

In case of extracting DFT data from *ab-initio* MD runs with **VASP**, positions, energy, and forces of each MD step 
can be obtained from ``vasprun.xml`` file as follows.
::

  $ python path/to/nap/nappy/vasp/vasprun2fp.py /path/to/dir/that/includes/vasprun.xml/


Then you get directories with names like ``#####`` including ``pos``, ``erg.ref``, and ``frc.ref`` files in it.

.. _prepare-pmd:

Prepare for pmd
---------------
When *fitpot* is performed, results of *pmd* will be used in each sample directory.
The *pmd* can be run by using ``parallel_run_pmd.py`` or ``run_pmd.sh``.
The ``cutoff_raidus`` should be the same or larger than the parameter *rcut* in ``in.params.NN`` file.

Thus, before going to run *fitpot*, users have to run *pmd* in each sample directory.
At first, to run *pmd*, make directory of name ``pmd`` in each sample directory as,
::

  $ for dir in smpl_*; do echo $dir; mkdir -p $dir/pmd; done

And if you are fitting some special potential that needs an auxiliary file like ``in.const.NN``, you need to copy such a file to each ``smpl_*/pmd/`` directory as,
::

  $ for dir in smpl_*; do echo $dir; cp in.const.NN $dir/pmd/; done


.. _prepare-scripts:

Make links of some scripts
-----------------------------------
In order to run ``fitpot`` program there must be the following scripts 
in the ``dataset/`` directory.
You can make links of the scripts as,
::

  $ ln -s /path/to/nap/fitpot/run_pmd.sh dataset/
  $ ln -s /path/to/nap/fitpot/serial_run_pmd.sh dataset/
  $ ln -s /path/to/nap/fitpot/parallel_run_pmd.py dataset/

Here you have to check whether you can run ``pmd`` correctly in every sample directory.
::

  $ ./run_pmd.sh in.params.NN smpl_*

Or if you have a lot of samples, you had better run *pmd* in parallel as,

.. code-block:: bash

  $ cat nodelist.txt
  node1
  node2
  node3
  node4
  $ python ./parallel_run_pmd.py in.params.NN

Then this python script assign nodes listed in ``nodelist.txt`` for the calculation of *pmd*.

.. _prepare-inputs:

Prepare input files
----------------------------------------
Inputs files needed for *fitpot* program are the following:

 * in.fitpot
 * dataset/in.params.NN

where ``NN`` indicates the name of interatomic potential used in *pmd* program.

You have to specify the ``num_samples`` in ``in.fitpot`` file 
which is a number of samples in ``dataset/`` directory.
The number of sample directories can be counted by the following command,

.. code-block:: bash

  $ ls | grep -e '^[0-9]....' | wc -l



.. _exec-fitpot:

Run *fitpot* program
------------------------------------
In the directory where ``dataset/`` directory and ``in.fitpot`` file exist,
you can run *fitpot* program as,
::

  $ ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

Or if you want it to run in parallel mode,
::

  $ mpirun -np 10 ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

There are some output files:

  ``out.erg.fin.[12]``
      These files include reference and *pmd* data of energies.
      To see whether the fitting went well or not, plot these data by using ``gnuplot`` as
      ::
         
         $ gnuplot
         gnuplot> plot 'out.erg.fin.1' us 1:2 w p t 'training set'
         gnuplot> rep 'out.erg.fin.2' us 1:2 w p t 'test set'


  ``out.frc.fin.[12]``
      These files include reference and *pmd* data of forces.


------------------------------

Input file for *fitpot*
================================

The following code shows an example of the input file ``in.fitpot``.
::

   num_samples      14
   test_ratio       0.1
   num_iteration    100
   num_iter_eval    1
   
   fitting_method   bfgs
   main_directory   data_set
   param_file       in.params.NN
   normalize_input  none
   
   force_match      true
   potential        NN

   ftol             1.0e-5
   xtol             1.0e-4
   
   penalty          none
   penalty_weight   1d-3
   
   # 1:Al, 2:Mg, 3:Si
   atom_energy  1  -0.19778
   atom_energy  2  -0.00074
   atom_energy  3  -0.80706




Input parameters for *fitpot*
----------------------------------------
Here are input parameters that users can change in *fitpot* program.

* :ref:`num_samples`
* :ref:`sample_list`
* :ref:`test_ratio`
* :ref:`num_iteration`
* :ref:`num_iter_eval`
* :ref:`fitting_method`
* :ref:`main_directory`
* :ref:`param_file`
* :ref:`ftol`
* :ref:`xtol`
* :ref:`force_match`
* :ref:`potential`
* :ref:`regularize`
* :ref:`penalty_weight`
* :ref:`sample_error`
* :ref:`atom_energy`
* :ref:`init_params`
* :ref:`init_params_sgm`
* :ref:`init_params_mu`
* :ref:`init_params_rs`


.. _num_samples:

num_samples
--------------------
Default: *no default*

Number of reference samples to be used for training and test.


.. _sample_list:

sample_list
--------------------
Default: *(blank)*

Path to the file that contains a list of samples to be used for training and test.
The format of the list file should be like,
::

   smpl_001
   smpl_002
   smpl_003
  ...

or with specifying which samples are training (1) or test (2) as,
::

   smpl_001  1
   smpl_002  2
   smpl_003  1
   ...

If whether training or test is specified in the list, `test_ratio` will be neglected.



.. _test_ratio:

test_ratio
--------------------
Default: *0.1*

The ratio of test data set :math:`r` within whole data set :math:`N`.
Thus the number of test data set is :math:`rN`, and the number of training data set is :math:`(1-r)N`.


.. _num_iteration:

num_iteration
--------------------
Default: *1*

Number of iterations of a minimization method.



.. _num_iter_eval:

num_iter_eval
--------------------
Test data set will be evaluated every *num_iter_eval* iterations.

Default: *1*


.. _fitting_method:

fitting_method
--------------------
Default: *test*

The method used to fit parameters to the sample data.
Available methods are the following:

*cg/CG* :
   Conjugate gradient algorithm which requires gradient information. **CG seems to be slower than BFGS and L-BFGS.**

*bfgs/BFGS* :
   Quasi-Newton method with BFGS. This requires gradient information.

*lbfgs/LBFGS* :
   Limited-memory version of BFGS. This requires gradient information.
   In case of optimizing a lot of parameters, this should be used instead of *BFGS*.

*check_grad* :
   Comparison of analytical derivative and numerical derivative.
   Use this to check the implemented analytical gradient.

*test/TEST* :
   Just calculate function L and gradient of L w.r.t. fitting parameters.



.. _main_directory:

main_directory
--------------------
Default: *dataset*

The directory that includes sample data. We call this ``dataset`` in the above instruction.

If you want to use ``..`` to specify the directory relative to the current working directory, e.g. ``../dataset``, you need to enclose with double-quotation marks like ``"../dataset"``.


.. _param_file:

param_file
--------------------
Default: *in.params.NN*

The name of the file that has parameter values in it. This is passed to ``pmd`` program.


.. _ftol:

ftol
-------
Default: *1.0e-5*

The tolerance of difference of the loss function value.

.. _xtol:

xtol
------
Default: *1.0e-4*

The tolerance of the change of variables which are optimized.
If either one of `ftol` or `xtol` is achieved, the optimization stops.


.. _force_match:

force_match
--------------------
Default: *False*

Whether or not to match forces. ( *True* or *False* )
It is highly recommended to match forces, since forces are important for molecular dynamics.



.. _potential:

potential
--------------------
Default: *none*

The potential whose parameters you are going to fit.
Now folloing potentials are available:

*NN*:
   Neural network potential



.. _regularize:

regularize
--------------------
Whether or not regularize bases obtained in *linreg* and *NN?* potentials. ( *True* or *False* )

Default: *False*


.. _penalty:

penalty
--------------------
Type of penalty term, *lasso* which is L1-norm penalty or *ridge* which is L2-norm penalty,
or *no* which means no penalty term.

Default: *no*



.. _penalty_weight:

penalty_weight
--------------------
The weight applied to the penalty term. This value also has to be determined through 
cross-validation scoring...

Default: *1.0*


.. _sample_error:

sample_error
------------------------------

Default: *0*

The number of samples whose errors are to be given. These errors appear at the denominators of energy and force in the evaluation function such that

.. math::

    \left( \frac{E^\mathrm{NN}-E^\mathrm{DFT}}{N_\mathrm{a}\varepsilon_\mathrm{e}}\right)^2 +\sum_i^{N_\mathrm{a}} \sum_\alpha^{xyz} \frac{1}{3N_\mathrm{a}}\left( \frac{F^\mathrm{NN}_{i\alpha} -F^\mathrm{DFT}_{i\alpha}}{\varepsilon_\mathrm{f}}\right)^2

If the difference between NN energy and DFT energy/force is lower than this value, this term becomes less than 1.0, which means the energy/force of the sample is thought to be converged.
The initial values of the errors are 0.001 (eV/atom) and 0.1 (eV/Ang) for energy and force, respectively.

There must be the same number of following entry lines as the above value which determine the errors of energy and force of each sample like the this,
::

  sample_error   2
      Al_fcc    0.001  0.2
      Al_bcc    0.001  0.2

The each entry has *entry_name*, *error of energy (eV/atom)* and *error of forces (eV/Ang)*.
The error values are applied to all the samples that contain *entry_name* in their directory names.

..
   .. _sample_weight:

   sample_weight
   --------------------
   Default: *False*

   Whether or not to apply weights to samples ( *True* or *False* ).




   .. _sample_weight_erg:

   sample_weight_erg
   --------------------
   Default: *1.0*

   Energy value :math:`E_\text{s}` in eV of the sample weight :math:`\exp (-\Delta E /E_\text{s})`.
   The :math:`\Delta E` is defined as the energy difference (per atom) from the most stable atomic energies.



.. _atom_energy:

atom_energy
--------------------
A DFT atomic energy that will be subtracted from the energies of sample structures.
Since the energy values of sample structures include the energies of atoms that are isolated 
in vacuum or gas phase.
The atomic energies of all atoms in the system should be specified in the following format:
::

  atom_energy   1   -0.808364
  atom_energy   2   -1.20934

where the first argument is species-ID and the second is the atomic energy of the species.

Default: *0.0*


--------------

.. _init_params:

init_params
--------------------
Default: *read*

Whether the paramters to be optimized are read from the file or initialized.

*read*:
   Read parameters from the file.

*gaussian*:
   Parameters are initialized with Gaussian distribution according *init_params_sgm* and *init_params_mu*.


.. _init_params_sgm:

init_params_sgm
--------------------
Default: *1d0*

Variance of Gaussian distribution of the initial values for parameters.


.. _init_params_mu:

init_params_mu
--------------------
Default: *0d0*

Mean value of Gaussian distribution of the initial values for parameters.

.. _init_params_rs:

init_params_rs
--------------------
Default: *12345.0*

Random seed for the initialization of parameters.


