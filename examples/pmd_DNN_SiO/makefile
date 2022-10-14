.phony: clean veryclean keep plot_erg test

clean:
	rm -f dump_*

veryclean: clean
	rm -f out.{erg,pmd} {erg,frc,strs}.pmd pmdfin

test:
	@mpirun -np 1 ../../pmd/pmd 2>&1 > out.pmd
	@grep Potential out.pmd.REF | head -n1 > ./out.REF
	@grep Potential out.pmd | head -n1 > ./out
	@diff -q out out.REF > /dev/null || (echo " FAIL: examples/pmd_DNN_SiO" && exit 1)
	@echo " PASS: examples/pmd_DNN_SiO"
	@rm -f out out.REF

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* data.* dump_* $(resdir)/

plot_erg:
	echo "plot 'out.erg' us 1:3 w l t 'total', '' us 1:4 w l t 'kinetic', '' us 1:5 w l t 'potential'" > plot_erg.gp
	gnuplot -e "load 'plot_erg.gp'; pause -1"
