.phony: clean veryclean keep test

clean:
	rm -f dump_* erg.pmd frc.pmd strs.pmd *~

veryclean: clean
	rm -f out.erg out.pmd.{2,4,8} pmdfin

test:
	@mpirun --oversubscribe -np 2 ../../pmd/pmd 2>&1 > out.pmd.2
	@mpirun --oversubscribe -np 4 ../../pmd/pmd 2>&1 > out.pmd.4
	@mpirun --oversubscribe -np 8 ../../pmd/pmd 2>&1 > out.pmd.8
	@grep Potential out.pmd.REF | tail -n1 > /tmp/out.REF
	@grep Potential out.pmd.2 | tail -n1 > /tmp/out.2
	@grep Potential out.pmd.4 | tail -n1 > /tmp/out.4
	@grep Potential out.pmd.8 | tail -n1 > /tmp/out.8
	@diff -q /tmp/out.2 /tmp/out.REF > /dev/null || (echo " FAIL: parallel 2 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 2 in examples/pmd_parallel"
	@diff -q /tmp/out.4 /tmp/out.REF > /dev/null || (echo " FAIL: parallel 4 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 4 in examples/pmd_parallel"
	@diff -q /tmp/out.8 /tmp/out.REF > /dev/null || (echo " FAIL: parallel 8 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 8 in examples/pmd_parallel"

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
