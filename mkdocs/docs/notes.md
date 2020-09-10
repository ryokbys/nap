Notes
=====

-   Temperature control with simple velocity scaling is added. So the
    format in input file \'in.pmd\' is changed a little.
-   van der Waals potential for Brenner is not correct, because
    parameters in the vdW potential are fitted with other potential.
-   Check num. of divisions and cutoff length in \'in.pmd\' before the
    simulation run.
-   Use variable array TAG instead of IS, because of the efficient
    parallelization. TAG includes species, index of FMV, total id.

History
-------

-   2014.02.11 Input and output files are seperated to folders of name
    of output number.
-   2014.01.?? Units are changed from atomic unit to eV, Angstrom, and
    mass of 1/12.
-   2012.05.22 Change input format to user readable one, and input file
    name became \'in.pmd\'.
-   2009.05.12 Use TAG instead of IS. Reduce num of MPI message passing
    in BAMOVE and BACOPY.
-   2009.04.28 Add Brenner potential with van der Waals term. Also
    simple velocity scaling is added.
-   2009.03.24 Add smoothing for embedded term, too.
-   2009.03.20 Bugs fixed about copy of IFMV in subroutine BAMOVE.
-   2009.03.20 Add smoothing for 2-body terms in EAM potential.
