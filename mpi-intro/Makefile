all: minmeanmax minmeanmax-allreduce minmeanmax-mpi
include ../Makefile.inc

minmeanmax: minmeanmax.f90
	${F77} ${FFLAGS} -o $@ $<

minmeanmax-mpi: minmeanmax-mpi.f90
	${MPIF77} ${FFLAGS} -g -o $@ $<

minmeanmax-allreduce: minmeanmax-allreduce.f90
	${MPIF77} ${FFLAGS} -o $@ $<


clean:
	rm -rf minmeanmax
	rm -rf minmeanmax-mpi
	rm -rf minmeanmax-allreduce
	rm -rf *.o
	rm -rf *~
