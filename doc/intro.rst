=============
Introduction
=============
*pmd* is an acronym of **parallel molecular dynamics** which
means that molecular dynamics (MD) using spatial decomposition technique
on parallel (ditributed-memory) computers.
And *pmd* is a part of the **Nagoya atomistic-simulation package (nap)**.

The main features of *pmd* are the following:

* several interatomic potentials for solid state systems are available;
* **Deep neural-network (DNN) interatomic potential**;
* **QEq** or **variable charge** Coulombic potential;
* parallel computation using spatial decomposition technique;
* efficient searching of neighbor atoms using linked-list cell method;
* structure relaxation using simple velocity damping or **FIRE** algorithm;
* thermostats: Berendsen and Langevin;
* barostat: Berendsen;
* **variable-timestep** MD for high-energy ion-bombardment simulation;
* **non-equilibrium MD (NEMD)** for heat flux simulation;
* **two-temperature model MD (TTM-MD)** for laser-ablation simulation.

Since this program has been developed for the purpose of personal research tool,
there are not so many functionalities. 
And there are some (open source) MD programs that can do almost the same thing that *pmd* 
can do. 
But there are some features only *pmd* or the parent packange **nap** can do.
Please feel free to contact me to ask anything
if you want to use *pmd* or **nap** for your specific purpose.

Bug reports and questions about *pmd* and **nap** are welcome,
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

For some analysis tools,
you may also need *Python 3.* and some python-utilities such as *numpy*, *scipy* and *pandas*.
