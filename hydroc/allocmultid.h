#ifndef __ALLOCMULTID_H
#define __ALLOCMULTID_H 1

int ***alloc3d_int(int n, int m, int l);
void free3d_int(int ***array, int n);

float ***alloc3d_float(int n, int m, int l);
void free3d_float(float ***array, int n);

int **alloc2d_int(int n, int m);
void free2d_int(int **array);

float **alloc2d_float(int n, int m);
void free2d_float(float **array);

#endif
