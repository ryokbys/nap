# What's pmd
pmd is a program package which enables us to perform parallel molecular dynamics simulation.
The package includes various interatomic potentials for metals and semiconductors.
pmd uses spatial decomposition technique for the parallelization, and cell list method for efficient neighbor search.

# Who made this?
* Ryo KOBAYASHI
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology. (2014-01-07)


# Usage
* make MD coordination file
  - compile mkconf
    $ make 10mkconf
  - run 
    $ ./10mkconf

* divide the coordination file to num of nodes
  - compile nconv
    $ make 20nconv
  - run
    $ ./20nconv
  - some inputs are required

* make config file of parallel md program
  - name of config file should be 'pmd.in'
  - compile parallel_md.f by using mpif90
    $ make pmd
  - run with mpirun (e.g., 8 nodes will be used.)
    $ mpirun -n 8 -machinefile hostfile ./a.out
    or
    $ ./30pmdrun.py

* combine pmd###-??? files to kvs???
    $ make 40combine
    $ ./40combine


# Notes:
* Temperature control with simple velocity scaling is added.
  So the format in input file 'pmd.in' is changed a little.
* van der Waals potential for Brenner is not correct,
  because parameters in the vdW potential are fitted with other potential.
* Check choices of force and mass before compilation.
* Check num. of divisions and cutoff length in 'pmd.in' before simulation run.
* Use variable array TAG instead of IS, because of fast parallelization.
  TAG includes species, index of FMV, total id.


# History
* 2012.05.22  Change input format to user readable one,
            and input file name became 'in.pmd'.
* 2009.05.12  Use TAG instead of IS. 
            Reduce num of MPI message passing in BAMOVE and BACOPY.
* 2009.04.28  Add Brenner potential with van der Waals term.
            Also simple velocity scaling is added.
* 2009.03.24  Add smoothing for embedded term, too.
* 2009.03.20  Bugs fixed about copy of IFMV in subroutine BAMOVE.
* 2009.03.20  Add smoothing for 2-body terms in EAM potential.
