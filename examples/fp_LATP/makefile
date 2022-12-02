.phony: clean veryclean run test
SHELL = /bin/bash

clean:
	rm -f in.vars.fitpot.[1-9]* *~

veryclean: clean
	rm -rf out.fp out.cs.* subdir_*

run:
	python ../../nappy/fitpot/fp.py --nproc 4 --random-seed 42 2>&1 | tee out.fp

test:
	@echo " Testing fp.py in examples/fp_LATP; this will take a few minutes..."
	@python ../../nappy/fitpot/fp.py --nproc 4 --random-seed 42 2>&1 > out.fp
	@python test.py 2&>1 > /dev/null || (echo " FAIL: examples/fp_LATP" && exit 1)
	@echo " PASS: examples/fp_LATP"

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* data.* pmd_*sh $(resdir)/
