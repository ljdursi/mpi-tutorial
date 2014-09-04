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
#include "ppm.h"
#include "plot.h"
#include "solver.h"
#include "allocmultid.h"

int main(int argc, char **argv) {
  int   n, nx, ny;
  float t, dt;
  int   iter;
  float ***u;

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

  nx = n+4; /* two cells on either side for BCs */
  ny = n+4;
  u = alloc3d_float(ny,nx,NVARS);

  initialconditions(u, nx, ny);
  outputppm(u,nx,ny,NVARS,"ics.ppm",IDENS);
  t=0.;
  for (iter=0; iter < 6*nx; iter++) {
      timestep(u,nx,ny,&dt);
      t += 2*dt;
      if ((iter % 10) == 1) {
        printf("%4d dt = %f, t = %f\n", iter, dt, t);
        plot(u, nx, ny);
      }
  }
  outputppm(u,nx,ny,NVARS,"dens.ppm",IDENS);
  closeplot();

  free3d_float(u,ny);
  return 0;
}
