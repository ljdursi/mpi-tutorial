include ../../Makefile.inc

.SUFFIXES:.mod

FFLAGS+=-g
hydro: hydro.o ppm.o solver.o plot.o
	$(MPIF90) -o $@ $^ $(FFLAGS) $(PGPLIBS)

hydro.o: hydro.f90 ppm.mod solver.mod plot.mod
	$(MPIF90) -c hydro.f90 $(FFLAGS)

ppm.mod: ppm.f90 solver.mod
	$(MPIF90) -c $^ $(FFLAGS)

plot.o: plot.f90 solver.mod
	$(MPIF90) -c $^ $(FFLAGS)

%.mod: %.f90
	$(MPIF90) -c $^ $(FFLAGS)

%.o: %.f90
	$(MPIF90) -c $^ $(FFLAGS)

clean:
	rm -f hydro *.o *.dat *.mod *~ *ppm 
