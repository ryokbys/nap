fname = sdrive.o lbfgs.o 

main : $(fname)
	gfortran $(fname) $(extras) -o sdrive

.f.o :; gfortran -c $*.f

