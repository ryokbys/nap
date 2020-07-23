.phony: clean veryclean keep

clean:
	rm -f dump_* erg.pmd frc.pmd str.pmd

veryclean: clean
	rm -f out.erg out.pmd pmdfin

resdir := result_$(shell date "+%y%m%d_%H%M")
keep:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin $(resdir)/
