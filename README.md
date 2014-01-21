# What's pmd
pmd is a program package which enables us to perform parallel molecular dynamics simulation.
The package includes various interatomic potentials for metals and semiconductors.
pmd uses spatial decomposition technique for the parallelization, and cell list method for efficient neighbor search.

# Who made this?
* Ryo KOBAYASHI
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology. (2014-01-07)


# Compilation
1. Make all the executable files.
  - Make the makefile

         $ ./configure FCFLAGS="-O4 -g"

    FCFLAGS will be added as an option to every compilation.
    So please change the FCFLAGS value according to your system.

  - Make all the executable files

         $ make all


# Usage
1. Make atom coordination file.
  - You have to make 10mkconf executable file in order to create the initial configuration of atoms.
    If there is any mkconf_* file that corresponds to the system you are going to simulate, you can use it and compile the file as following,

         $ make 10mkconf

  - Run the 10mkconf command,

         $ ./10mkconf

    then you get a atom coordination file named 'pmd00000-0000'.

2. Divide the coordination file to num of nodes.
  - Run '20nconv' command and this program leads you interactively to divide the coordination file to parallel configuration.

         $ ./20nconv

3. Make config file of parallel md program.
  - Name of config file should be in.pmd'.
  - Run 'pmd' program with mpirun. For example, running pmd with 8 nodes, as follows

         $ mpirun -n 8 -machinefile hostfile ./a.out

    or

         $ ./30pmdrun.py

4. Combine pmd###-??? files to akr??? for the visualization using Akira.

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
