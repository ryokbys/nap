.phony: clean veryclean keep test

clean:
	rm -f dump_* erg.pmd frc.pmd strs.pmd

veryclean: clean
	rm -f out.erg out.pmd pmdfin

test:
	@mpirun -np 1 ../../pmd/pmd 2>&1 > out.pmd
	@grep Potential out.pmd.REF | head -n1 > /tmp/out.REF
	@grep Potential out.pmd | head -n1 > /tmp/out
	@diff -q /tmp/out /tmp/out.REF > /dev/null || (echo " FAIL: examples/pmd_BVS_LLZO" && exit 1)
	@echo " PASS: examples/pmd_BVS_LLZO"

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
