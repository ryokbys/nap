---
title: 'nap: A molecular dynamics package with potential-parameter fitting programs'
tags:
  - Fortran
  - Python
  - materials science
  - molecular dynamics
  - interatomic potential
  - neural-network potential
  - meta-heuristics
authors:
  - name: Ryo KOBAYASHI
    orcid: 0000-0001-8244-5844
    affiliation: 1
affiliations:
 - name: Department of Physical Science and Engineering, Nagoya Institute of Technology
   index: 1
date: 28 August 2020
bibliography: paper.bib
---

# Summary

The **nap** is a package for molecular dynamics (MD) simulation consisting of an MD program (**pmd**) with a spatial-decomposition technique using the MPI library and two parameter-optimization programs: one for classical (CL) potentials (**fp.py**) and another for machine-learning (ML) potentials (**fitpot**).
Since the numbers of parameters to be optimized are much different between CL and ML potentials, optimization approaches for them are also different; meta-heuristic global minimum-search algorithms for the CL potentials where the numbers of parameters are usually much less than one hundred, and gradient-based methods for the ML potentials.
The parameters of CL potentials can be optimized to any target quantity that can be computed using the potentials since meta-heuristic methods do not require the derivatives of the quantity with respect to parameters. On the other hand, ML-potential parameters can be optimized to only energies, forces on atoms and stress components of reference systems, mainly because gradient-based methods require derivatives of other quantities with respect to parameters, and the analytical derivatives and the coding of them are usually painful and sometimes impossible.
Potentials can be used in combination with any other potential, such as pair and angular potentials, short-range and long-range potentials, CL and ML potentials.
With using the **nap** package, users can perform MD simulation of solid-state materials with the choice of different levels of complexity (CL or ML) by creating interatomic potentials optimized to quantum-mechanical calculation data even if no potential is available.

# Statement of need

MD simulation is widely used in many fields such as materials science, chemistry, physics, etc., to study dynamics of atoms or molecules. In order to perform MD simulation of systems including large number of atoms, where quantum-mechanical calculations can not be applied to compute interatomic interactions, empirical interatomic potentials between species are required. And the results of MD simulation are strongly dependent on the property or accuracy of the potentials used in the simulation. Hence, there are a lot of CL potential models have been proposed such as Lennard-Jones (LJ) potential for van der Waals interaction, Coulombic potential for ionic interaction, Morse potential for covalent interactions [@Morse1929-ow], angular-dependent models for angles between covalent bonds, bond-order models for more complex systems, etc.
Recently ML potentials have been also actively studied because they are usually more flexible and can reproduce reference data more precisely than CL potentials.
Even though the potential is flexible or suitable to problems considered, the parameters in the potential model still need to be optimized to well reproduce the properties or phenomena that are in focus.

There are already a lot of MD programs and some of them can use both CL and ML potentials, such as LAMMPS [@Plimpton1995-az] and IMD [@Stadler1997-wr]. And also there are some parameter-optimization programs that can produce parameter sets available in the other MD programs such as aenet [@Artrith2016-mu] for ML potentials and potfit [@Brommer2015-hw] for CL potentials.
However, there is also a demand of combining ML potentials and simpler CL potentials, e.g., ML potential with Coulomb interactions [@Morawietz2013-qq] and ML potential with core-core repulsions [@Wang2019-py], since creating an ML potential that covers very short-range and/or very long-range interactions is very inefficient. Thus it is beneficial if the programs of parameter-optimization for both CL and ML potentials are in one package and highly connected to one MD program and it will be more efficient than using several different programs to optimize parameters of CL and ML potentials.


# Programs and functionalities

- **pmd**: A Fortran program to perform MD simulation. It can perform large-scale MD simulation efficiently by using linked-cell list [@Allen2017-nw] and spatial decomposition with MPI parallelization. Several widely-used interatomic potentials for solid-state materials are implemented, for example, pair potentials such as LJ, Morse, Coulomb and screened Coulomb; angular-dependent potentials such as Stillinger-Weber [@Stillinger1985-gy] and Tersoff [@Tersoff1988-bo]; multi-body potentials such as Finnis-Sinclair [@Finnis1984-hv] and embeded-atom method [@Daw1984-hj]; machine-learning potentials such as linear-regression [@Seko2014-mh] and neural-network (NN) [@Behler2007-hr].
- **fp.py**: A Python program to optimize the parameters of classical potentials using meta-heuristic algorithms such as cuckoo search [@Yang2009-at]. Meta-heuristic algorithms need to perform MD simulations and to evaluate the quantities to be compared with target values of each individual representing a candidate parameter set. These jobs by individuals are performed as child processes of the program by calling a shell script that describes what to compute using given parameters. This program has a functionality of automatic update of the search range of parameters and it allows to optimize parameters efficiently and less dependent on the initial ranges of parameters [@Kobayashi2020-rz].
- **fitpot**: A Fortran program to optimize the parameters of ML potentials, such as linear regression and NN, using gradient-based methods. It can use high-performance computer resources to evaluate energies, forces and stresses of a huge number of sample systems in parallel using the MPI library, which is sometimes crucial since in general ML models require a lot of sampling to avoid over-fitting and make them robust.


# Acknowledgements

This work was supported in part by JSPS KAKENHI Grant Number JP20H05290 (Grant-in-Aid for Scientific Research on Innovative Areas "Interface Ionics") and by the "Materials Research by Information Integration" Initiative (MI2I) project of the "Support Program for Starting Up Innovation Hub" from JST.

# References
