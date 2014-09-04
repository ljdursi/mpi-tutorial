/*  Original F90 code  copyright (C) 2001,2003, Ue-Li Pen
 *  written November 2001 by Ue-Li Pen, pen@cita.utoronto.ca
 *
 * see http://arxiv.org/abs/astro-ph/0305088
 * or http://www.cita.utoronto.ca/~pen/MHD
 *
 * This code is licencensed under the GPL
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "ppm.h"
#include "plot.h"
#include "solver.h"
#include "allocmultid.h"

int main(int argc, char **argv) {
  int   n, nx, ny;
  float t, dt;
  int   iter;
  float ***u;
  int ierr, size, rank;
  MPI_Comm commcart;
  int gridsize[2];
  int gridcoords[2];
  int periods[2] = {1,1};
  int locnx;
  int left,right;
  int start,end;
  char filename[30];

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
  gridsize[0] = size; gridsize[1] = 1;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, gridsize,
                       periods, 1, &commcart);
  ierr = MPI_Comm_rank(commcart, &rank);
  ierr = MPI_Cart_shift(commcart, 0, 1, &left, &right);
  ierr = MPI_Cart_coords(commcart, rank, 2, gridcoords);

  if (argc != 2) {
      printf("Usage: %s n, where n is size of domain.\n", argv[0]);
      exit(-1);
  } else {
      n = atoi(argv[1]);
      if (n < 2 ||  n > 500) {
         n = 100;
         printf("%s: invalid n = %s; using %d.\n", argv[0], argv[1], n);
      }
  }

  /* calculate an nx and ny from n, rank */
  locnx = n/gridsize[0];
  start = rank*locnx;
  end   = rank*locnx + locnx-1;
  if (gridcoords[0] == gridsize[0]-1) end = n-1;
  locnx = end-start+1;
 
  nx = locnx+2*NGUARD; /* guardcells on either side for BCs */
  ny = n+2*NGUARD;
  u = alloc3d_float(ny,nx,NVARS);

  initialconditions(u, nx, ny, n, start);
  sprintf(filename,"%d-ics.ppm",rank);
  outputppm(u,nx,ny,NVARS,filename,IDENS);
  //outputppm_mpiio(u,nx,ny,n,start,NVARS,"ics.ppm",IDENS);
  t=0.;
  for (iter=0; iter < 12*n; iter++) {
      timestep(u,nx,ny,&dt,commcart,left,right);
      t += 2*dt;
      if ((iter % 10) == 0) {
        if (rank == 0) 
            printf("%4d dt = %f, t = %f\n", iter, dt, t);
        plot(u, nx, ny, start);
      }
  }
  sprintf(filename,"%d-dens.ppm",rank);
  outputppm(u,nx,ny,NVARS,filename,IDENS);
  //outputppm_mpiio(u,nx,ny,n,start,NVARS,"dens.ppm",IDENS);
  closeplot();

  free3d_float(u,ny);
  MPI_Finalize();
  return 0;
}
