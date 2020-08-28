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
date: 21 July 2020
bibliography: paper.bib
---

# Summary

Molecular dynamics (MD) simulation is widely used to study in many fields such as materials science, chemistry, physics, etc., where dynamics of atoms or molecules is of interest. In order to perform MD simulation of systems including large number of atoms, where *ab-initio* calculations can not be applied to compute interatomic interactions, empirical interatomic potentials between species are required. And the results of MD simulation are heavily dependent on the property or accuracy of the potentials used in the simulation. Hence, there are a lot of potential models have been proposed such as Lennard-Jones (LJ) potential for van der Waals interaction, Coulombic potential for ionic interaction, Morse potential for covalent bond, angular-dependent models for angles between covalent bonds, bond-order models for more complex systems, etc. Recently a new class of potential models, so-called machine-learning (ML) potentials, have been also actively studied because of ML models are usually more flexible and can reproduce reference data more precisely than the classical models. Even though the potential is flexible or suitable to problems considered, the parameters in the potential model should be optimized to well reproduce the properties or phenomena that are in focus.

The *nap* package contains parameter optimization programs of classical potentials and neural-network (NN) potential(REF) (one of ML potentials) that can be used in an MD program *pmd* which is also included in the package. Since the number of parameters to be optimized are much different between classical potentials and NN potential, optimization approaches are different; meta-heuristic global minimum search algorithms for the classical potentials where the number of parameters are usually less than hundred, and gradient-based methods such as quasi-Newton method for the NN potential.
The paremeters of classical potentials can be optimized to any target quantity that can be computed using the potentials since meta-heuristic methods do not require the derivatives of the quantity with respect to parameters. On the other hand, NN-potential parameters can be optimized to only energies, forces on atoms and stress components of reference systems, mainly because gradient-based methods require derivatives of other quantities with respect to parameters, and the analytical derivative and its coding are usually painful and sometimes impossible.

# Programs and functionalities

- **pmd**: A Fortran program to perform MD simulation. It can perform large-scale MD simulation efficiently by using linked-cell list and spatial decomposition with MPI parallelization. Several widely-used interatomic potentials for solid-state materials are implemneted, for example, pair potentials such as Lennard-Jones, Morse, Coulomb and screened Coulomb; angular-dependent potentials such as Stillinger-Weber and Tersoff; multi-body potentials such as Finnis-Sinclair and embeded-atom method; machine-learning potentials such as linear-regression and neural-network. 
- **fp.py**: A Python program to optimize the parameters of classical potentials using meta-heuristic algorithms. Meta-heuristic algorithms require to perform MD simulations and to evaluate the quantities to be compared with target values for each individual representing a candidate parameter set. These works are performed as child processes of the program by calling a shell script that describes what to compute using given parameters. This program has a functionality of automatic update of the search range of parameters and it allows to optimize parameters efficiently and less dependent on the initial ranges of parameters.(REF)
- **fitpot**: A Fortran program to optimize the parameters of NN potential using gradient-based methods. It can use high-performance computer resources to evaluate energies, forces and stresses of huge number of sample systems in parallel using MPI library, which is sometimes crucial since in general ML models require a lot of sampling to avoid over-fitting and make them robust.


# Statement of need

There are already many MD programs and some of them can use both classical and ML potentials, such as LAMMPS.(REF) And there are also some parameter-optimization programs that can produce parameter sets available in other MD programs such as aenet (REF) for NN potentials and potfit (REF) for classical potentials. 

# Acknowledgements

This work was supported in part by JSPS KAKENHI Grant Number JP20H05290 (Grant-in-Aid for Scientific Research on Innovative Areas "Interface Ionics") and by the "Materials Research by Information Integration" Initiative (MI2I) project of the "Support Program for Starting Up Innovation Hub" from JST.

# References
