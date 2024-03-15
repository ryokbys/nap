.phony: clean veryclean keep test

clean:
	rm -f dump_* erg.pmd frc.pmd strs.pmd *~

veryclean: clean
	rm -f out.erg out.pmd.{2,4} pmdfin

test:
	@OMP_NUM_THREADS=1 mpirun -np 2 ../../pmd/pmd 2>&1 > out.pmd.2
	@OMP_NUM_THREADS=1 mpirun -np 4 ../../pmd/pmd 2>&1 > out.pmd.4
	@grep Potential out.pmd.REF | tail -n1 > out.REF
	@grep Potential out.pmd.2 | tail -n1 > out.2
	@grep Potential out.pmd.4 | tail -n1 > out.4
	@diff -q out.2 out.REF > /dev/null || (echo " FAIL: parallel 2 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 2 in examples/pmd_parallel"
	@diff -q out.4 out.REF > /dev/null || (echo " FAIL: parallel 4 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 4 in examples/pmd_parallel"
	@rm -f out.{REF,2,4}

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
