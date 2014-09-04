include ../Makefile.inc

.SUFFIXES:.mod

FFLAGS+=-g
hydro: hydro.o ppm.o solver.o plot.o
	$(FC) -o $@ $^ $(FFLAGS) $(PGPLIBS)

hydro.o: hydro.f90 ppm.mod solver.mod plot.mod
	$(FC) -c hydro.f90 $(FFLAGS)

%.mod: %.f90
	$(FC) -c $^ $(FFLAGS)

clean:
	rm -f hydro *.o *.dat *.mod *~ *ppm 
