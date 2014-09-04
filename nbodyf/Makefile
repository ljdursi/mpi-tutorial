PROGS = nbody nbody-allgather

all: $(PROGS)

include ../Makefile.inc

.PREFIX: .f90 .mod

%.o: %.f90
	$(F77) -c $(FFLAGS) $< 

nbody: nbody.o 
	$(F77) $(FFLAGS) nbody.o -o $@  $(FLIBS) $(PGPLIBS) -lm

nbody-allgather.o: nbody-allgather.f90
	$(MPIF77) $(FFLAGS) nbody-allgather.f90 -c -o $@ 

nbody-allgather: nbody-allgather.o 
	$(MPIF77) $(FFLAGS) nbody-allgather.o -o $@  $(FLIBS) $(PGPLIBS) -lm

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(PROGS)
	rm -f *~
