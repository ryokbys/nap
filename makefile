.phony: all clean pmd fitpot test
SHELL = /bin/bash

all: pmd fitpot test

clean:
	(cd pmd/ && make clean)
	(cd fitpot/ && make clean)

pmd:
	(cd pmd/ && make pmd)

fitpot: pmd
	(cd fitpot/ && make fitpot)

test: pmd fitpot
	@(cd pmd/ && make test)
	@(cd fitpot/ && make test)
	@echo " Run make test-fp after setting up PYTHONPATH for nappy."
	@echo " To setup PYTHONPATH for nappy, see the documentation:"
	@echo "   http://ryokbys.web.nitech.ac.jp/contents/nap_docs/install.html#setup-nappy-required-for-fppy"

test-fp:
	@(cd examples/fp_LZP && make test)

