.phony: clean veryclean run extract keep test

clean:
	rm -f in.params.DNN.* dir_list.txt out.{iter,energy,force,stress} out.*fin
	rm -f *~ fp_* err*

veryclean: clean
	rm -f out.fitpot

NP = 2
run:
	mpirun --oversubscribe -np $(NP) ~/src/nap/fitpot/fitpot | tee out.fitpot

test:
	@/bin/bash ./test_result.sh > /dev/null || (echo " FAIL: examples/fitpot_DNN_SiO" && exit 1)
	@echo " PASS: examples/fitpot_DNN_SiO"

extract: out.fitpot
	grep 'iter,' out.fitpot > out.iter
	grep 'ENERGY:' out.fitpot > out.energy
	grep 'FORCE:' out.fitpot > out.force
	grep 'STRESS:' out.fitpot > out.stress

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* makefile job.*sh $(resdir)/


