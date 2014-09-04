#include <stdio.h>
#include "allocmultid.h"

void color(float val, float lo, float hi, int rgb[3])  {
  float scaleval = (val-lo)/(hi-lo) * 255.;

  rgb[0] = (int)(scaleval*3);    
  if (rgb[0] > 255) rgb[0] = 255;

  rgb[1] = (int)(scaleval);
  rgb[2] = 0;
}

void outputppm(float ***data, int nx, int ny, int nvars, const char *filename, const int varnum) {
  FILE *out;
  int ***rgb;
  float lo, hi;
  int i, j;

  lo = 1e19;
  hi = -1e19;
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      if (data[j][i][varnum] < lo) lo = data[j][i][varnum];
      if (data[j][i][varnum] > hi) hi = data[j][i][varnum];
    }
  }

  rgb = alloc3d_int(nx,ny,3);
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      color(data[j][i][varnum], lo, hi, &(rgb[i][j][0]));
    }
  }

  out = fopen(filename,"wb");
  fprintf(out, "P6\n%d %d\n255\n", nx, ny);
  for (j=ny-1; j>=0; j--) {
    for (i=0; i<nx; i++) {
      fwrite(rgb[i][j], 1, 3,  out);
    }
  }
  fclose(out);
 
  free3d_int(rgb,nx);
}
