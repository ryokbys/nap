=============
Introduction
=============
*pmd* is an acronym of **parallel molecular dynamics** which
means that molecular dynamics (MD) using spatial decomposition technique
on parallel (ditributed-memory) computers.
And *pmd* is a part of the **Nagoya atomistic-simulation package (NAP)**.

The main features of *pmd* are the following:

* some interatomic potentials for solid state systems that does not include long-range one such as Coulombic potential are available;
* **neural-network (NN) interatomic potential** for metal alloy systems;
* parallel computation using spatial decomposition technique;
* efficient searching of neighbor atoms using linked-list cell method;
* structure relaxation using simple velocity damping or **FIRE** algorithm;
* NVT ensemble simulation using Berendsen thermostat;
* NpT ensemble simulation using Berendsen thermo- and baro-stat;
* non-equilibrium MD (NEMD) for heat flux simulation;
* capable of fitting potential parameters with minimization algorithms such as steepest descent, 
  conjuget gradient, and quasi-Newton (BFGS).

Since this program has been developed for the purpose of personal research tool,
there are not so many functionalities. 
And there are some (open source) MD programs that can do almost the same thing that *pmd* 
can do. 
But there are some features only *pmd* or the parent packange **NAP** can do.
Please feel free to contact me to ask anything
if you want to use *pmd* or **NAP** for your specific purpose.

Bug reports and questions about *pmd* and **NAP** are welcome,
but I am afraid that I might not be able to respond all the reports or questions.


Requirements
====================
*pmd* can be executed in,

* Linux;
* MacOS X;
* Windows (Cygwin).

You also need,

* Fortran compiler;
* MPI library;

installed in the system you want to run ``pmd``.

For some analyses tools,
you may also need *Python 2.7* and some python-utilities such as *numpy* and *scipy*.
