#include <stdio.h>
#include <mpi.h>
#include "allocmultid.h"
#include "solver.h"

void color(float val, float lo, float hi, char *rgb)  {
  float scaleval = (val-lo)/(hi-lo) * 255.;
  int red;
  int green;

  red = (int)(scaleval*3);
  if (red > 255) red = 255;
  rgb[0] = (char)red;

  green = (int)scaleval;

  rgb[1] = (char)green;
  rgb[2] = (char)0;
}

void outputppm(float ***data, int nx, int ny, int nvars, const char *filename, const int varnum) {
  FILE *out;
  char ***rgb;
  float lo, hi;
  float globlo, globhi;
  int i, j;

  lo = 1e19;
  hi = -1e19;
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      if (data[j][i][varnum] < lo) lo = data[j][i][varnum];
      if (data[j][i][varnum] > hi) hi = data[j][i][varnum];
    }
  }

  MPI_Allreduce(&lo,&globlo, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&hi,&globhi, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  rgb = alloc3d_char(ny,nx,3);
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      color(data[j][i][varnum], globlo, globhi, &(rgb[ny-1-j][i][0]));
    }
  }

  out = fopen(filename,"wb");
  fprintf(out, "P6\n%d %d\n255\n", nx, ny);
  fwrite(&(rgb[0][0][0]),3,ny*nx,out);
  fclose(out);
 
  free3d_char(rgb,ny);
}

void outputppm_mpiio(float ***data, int nx, int ny, int n, int start, int nvars, char *filename, const int varnum) {
  FILE *out;
  char ***rgb;
  float lo, hi;
  float globlo, globhi;
  int i, j; 
  int globalarray[]={n,3*n};
  int localarray[]={n,3*(nx-2*NGUARD)};
  int starts[]={0,3*start};
  MPI_Datatype fileregion;
  int rank;
  MPI_File fh;
  int offset;
  MPI_Status status;
  
  MPI_Type_create_subarray(2, globalarray, localarray, starts, MPI_ORDER_C, MPI_CHARACTER, &fileregion);
  MPI_Type_commit(&fileregion);

  lo = 1e19;
  hi = -1e19;
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      if (data[j][i][varnum] < lo) lo = data[j][i][varnum];
      if (data[j][i][varnum] > hi) hi = data[j][i][varnum];
    }
  }

  MPI_Allreduce(&lo,&globlo, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&hi,&globhi, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  rgb = alloc3d_char(n, nx-2*NGUARD, 3);
  for (j=NGUARD; j<ny-NGUARD; j++) {
     for (i=NGUARD; i<nx-NGUARD; i++) {
      color(data[j][i][varnum], globlo, globhi, &(rgb[ny-NGUARD-1-j][i-NGUARD][0]));
    }
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
      out = fopen(filename,"wb");
      fprintf(out, "P6\n%d %d\n255\n", n, n);
      fclose(out);
  }

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
  offset = 18;
  
  MPI_File_set_view(fh, offset, MPI_CHARACTER, fileregion, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, &(rgb[0][0][0]), n*(nx-2*NGUARD)*3, MPI_CHARACTER, &status);
  MPI_File_close(&fh);
 
  free3d_char(rgb,ny-2*NGUARD);
}
