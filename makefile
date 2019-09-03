
.phony: all clean pmd fitpot FORCE

all: pmd fitpot

clean: FORCE
	cd pmd; make clean
	cd fitpot; make clean

pmd: FORCE
	cd pmd; make pmd

fitpot: FORCE
	cd fitpot; make fitpot

FORCE:
