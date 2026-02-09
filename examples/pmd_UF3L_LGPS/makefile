.phony: clean veryclean save plot_erg

NMPI := 1
NOMP := 1

clean:
	rm -f dump_[0-9]* pmd_[0-9]*

veryclean: clean
	rm -f out.* {erg,frc,strs}.pmd traj.extxyz

run:
	OMP_NUM_THREADS=$(NOMP) mpirun -np $(NMPI) pmd | tee out.pmd

test:
	@OMP_NUM_THREADS=$(NOMP) mpirun -np $(NMPI) ../../pmd/pmd 2>&1 > out.pmd
	@grep Potential out.pmd.REF | head -n1 > out.REF
	@grep Potential out.pmd | head -n1 > out
	@diff -q out out.REF > /dev/null || (echo " FAIL: examples/pmd_UF3L_LGPS" && exit 1)
	@echo " PASS: examples/pmd_UF3L_LGPS"
	@rm -f out out.REF

resdir := result_$(shell date "+%y%m%d_%H%M")
save:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin traj.extxyz $(resdir)/
