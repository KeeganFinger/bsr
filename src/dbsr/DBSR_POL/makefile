add_executable(dbsr_pol3
mod_pol.f90              dbsr_pol.f90             add_nortb.f90
read_data.f90            r_arg.f90                r_dbsrmat.f90
r_dipmat.f90             solv_mat.f90)

L  =  $(DBS) $(ZCONFJJ) $(ZCOM) $(LAPACK)

O = $(S:.f90=.o)

$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod

clean:
	rm -f *.o *.mod

x:
	rm -f $(N)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o

.f90.o:
	$(FC) $(comp)  $<

# -fcheck=bounds -fbacktrace
