.phony: clean veryclean keep test

clean:
	rm -f dump_* erg.pmd frc.pmd strs.pmd *~

veryclean: clean
	rm -f out.erg out.pmd.{1,2,4} pmdfin

test:
	@OMP_NUM_THREADS=1 mpirun -np 1 ../../pmd/pmd 2>&1 > out.pmd.1
	@OMP_NUM_THREADS=1 mpirun -np 2 ../../pmd/pmd 2>&1 > out.pmd.2
	@OMP_NUM_THREADS=1 mpirun -np 4 ../../pmd/pmd 2>&1 > out.pmd.4
	@grep Potential out.pmd.REF | tail -n1 > out.REF
	@grep Potential out.pmd.1 | tail -n1 > out.1
	@grep Potential out.pmd.2 | tail -n1 > out.2
	@grep Potential out.pmd.4 | tail -n1 > out.4
	@awk 'NR==FNR{ref=$$3;next}{val=$$3;d=val-ref;if(d<0)d=-d;if(d>0.005){exit 1}}' \
	  out.REF out.1 || (echo " FAIL: parallel 1 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 1 in examples/pmd_parallel"
	@awk 'NR==FNR{ref=$$3;next}{val=$$3;d=val-ref;if(d<0)d=-d;if(d>0.005){exit 1}}' \
	  out.REF out.2 || (echo " FAIL: parallel 2 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 2 in examples/pmd_parallel"
	@awk 'NR==FNR{ref=$$3;next}{val=$$3;d=val-ref;if(d<0)d=-d;if(d>0.005){exit 1}}' \
	  out.REF out.4 || (echo " FAIL: parallel 4 in examples/pmd_parallel" && exit 1)
	@echo " PASS: parallel 4 in examples/pmd_parallel"
	@rm -f out.{REF,2,4}

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
