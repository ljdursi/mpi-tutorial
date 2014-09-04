CPROGS = nbody nbody_omp

all: $(CPROGS)

include ../Makefile.inc

OPENMP=-fopenmp 
UTILS=allocmultid.o

%.o:%.c
	$(CC) -g $(CFLAGS) $(OPENMP) -c $< -o $@

nbody: nbody.o allocmultid.c

nbody-allgather.o: nbody-allgather.c 
	$(MPICC) $(CFLAGS) -c $<

nbody-pipeline.o: nbody-pipeline.c 
	$(MPICC) $(CFLAGS) -c $< 

nbody-nonblocking-pipeline.o: nbody-nonblocking-pipeline.c
	$(MPICC) $(CFLAGS) -c $<

nbody-gridparticles.o: nbody-gridparticles.c
	$(MPICC) $(CFLAGS) -c $<

nbody: nbody.o $(UTILS)
	$(F77) $(CFLAGS) nbody.o -g $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

nbody_omp: nbody_omp.o $(UTILS)
	$(F77) $(CFLAGS) nbody_omp.o -g -fopenmp $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

nbody-allgather: nbody-allgather.o $(UTILS) $(MPIUTILS)
	$(MPIF77) $(CFLAGS) nbody-allgather.o $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

nbody-pipeline: nbody-pipeline.o $(UTILS) $(MPIUTILS)
	$(MPIF77) $(CFLAGS) nbody-pipeline.o $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

nbody-nonblocking-pipeline: nbody-nonblocking-pipeline.o $(UTILS) $(MPIUTILS)
	$(MPIF77) $(CFLAGS) nbody-nonblocking-pipeline.o $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

nbody-gridparticles: nbody-gridparticles.o $(UTILS) $(MPIUTILS)
	$(MPIF77) $(CFLAGS) nbody-gridparticles.o $(UTILS) -o $@  $(CLIBS) $(PGPLIBS) -lm

clean:
	rm -f *.o
	rm -f $(CPROGS)
	rm -f myparticles.txt
