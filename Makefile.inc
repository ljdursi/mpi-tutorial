ROOTDIR=${HOME}/mpi

USEPGPLOT = -DPGPLOT
PGPLIBS   = -L/usr/X11R6/lib64 -L${PGPLOT_DIR}   -lcpgplot -lpgplot -lX11  -lpng

INCLUDES=-I${PGPLOT_DIR}

CC=gcc
F77=gfortran
F90=gfortran
CFLAGS=-std=c99 -O3 -Wall ${USEPGPLOT} ${INCLUDES}
FFLAGS=-O3 -Wall ${INCLUDES} -g -fbounds-check
CFLAGS_FAST=-std=c99 -O3 ${USEPGPLOT} ${INCLUDES}
CLIBS=-lm
LDFLAGS=-lm

OMPCC=gcc
OMPF77=gfortran
OMPCFLAGS=${CFLAGS} -fopenmp
OMPCFLAGS_FAST = ${CFLAGS_FAST} -fopenmp
OMPFFLAGS=-fopenmp
OMPLDFLAGS=${LDFLAGS} -fopenmp

MPICC=mpicc
MPIF90=mpif90
MPIF77=mpif90
MPIUTILS=
