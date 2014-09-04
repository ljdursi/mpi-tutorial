include ../Makefile.inc

hydro: hydro.o ppm.o solver.o allocmultid.o plot.o
	$(FC) -o $@ $^ $(FFLAGS) $(CLIBS) $(PGPLIBS)

clean:
	rm -f hydro *.o *.dat *~ *pbm *ppm
