.phony: clean veryclean keep test

clean:
	rm -f dump_* erg.pmd frc.pmd str.pmd

veryclean: clean
	rm -f out.erg out.pmd pmdfin

test:
	@mpirun -np 1 ../../pmd/pmd 2>&1 > out.pmd
	@grep Potential out.pmd.REF | head -n1 > out.REF
	@grep Potential out.pmd | head -n1 > out
	@diff -q out out.REF > /dev/null || (echo " FAIL: examples/pmd_W" && exit 1)
	@rm -f out.REF out
	@echo " PASS: examples/pmd_W"

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
