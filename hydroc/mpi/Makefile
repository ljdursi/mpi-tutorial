include ../../Makefile.inc

CC=mpicc
CFLAGS+=-g -Wall -O0

hydro: hydro.o ppm.o solver.o allocmultid.o plot.o
	$(MPICC) -o $@ $^ $(FFLAGS) $(CLIBS) $(PGPLIBS) -lgfortran -lm

clean:
	rm -f hydro *.o *.dat *~ *pbm *ppm
