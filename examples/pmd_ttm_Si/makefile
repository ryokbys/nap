.phony: clean veryclean keep test

clean:
	rm -f dump_* out.Te_* *~

veryclean: clean
	rm -rf out.{pmd,erg,eio_ttm} outs* {erg,frc,strs}.pmd graph_*png pmdfin

test:
	@mpirun -np 1 ../../pmd/pmd 2>&1 > out.pmd
	@grep Potential out.pmd.REF | head -n1 > /tmp/out.REF
	@grep Potential out.pmd | head -n1 > /tmp/out
	@diff -q /tmp/out /tmp/out.REF > /dev/null || (echo " FAIL: examples/pmd_ttm_Si" && exit 1)
	@echo " PASS: examples/pmd_ttm_Si"

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp -r in.* out.* graph_*png outs* pmdini pmdfin dump_* $(resdir)/
