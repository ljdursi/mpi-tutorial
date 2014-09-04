#include <stdlib.h>

float ***alloc3d_float(int n, int m, int l) {
    float *data = malloc(n*m*l*sizeof(float));
    if (!data) return NULL;

    float ***array = (float ***)malloc(n*sizeof(float **));
    for (int i=0; i<n; i++) {
        array[i] = (float **)malloc(m*sizeof(float *));
        for (int j=0; j<m; j++) {
            array[i][j] = &(data[(i*m+j)*l]);
        }
    } 
    return array;
}

void free3d_float(float ***array, int n) {
    free(&(array[0][0][0]));
    for (int i=0; i<n; i++) {
        free(array[i]);
    } 
    free(array);
}

int ***alloc3d_int(int n, int m, int l) {
    int *data = malloc(n*m*l*sizeof(int));
    if (!data) return NULL;

    int ***array = (int ***)malloc(n*sizeof(int **));
    for (int i=0; i<n; i++) {
        array[i] = (int **)malloc(m*sizeof(int *));
        for (int j=0; j<m; j++) {
            array[i][j] = &(data[(i*m+j)*l]);
        }
    } 
    return array;
}

void free3d_int(int ***array, int n) {
    free(&(array[0][0][0]));
    for (int i=0; i<n; i++) {
        free(array[i]);
    } 
    free(array);
}

float **alloc2d_float(int n, int m) {
    float *data = malloc(n*m*sizeof(float));
    if (!data) return NULL;

    float **array = (float **)malloc(n*sizeof(float *));
    for (int i=0; i<n; i++) {
        array[i] = &(data[i*m]);
    } 
    return array;
}

void free2d_float(float **array) {
    free(&(array[0][0]));
    free(array);
}

int **alloc2d_int(int n, int m) {
    int *data = malloc(n*m*sizeof(int));
    if (!data) return NULL;

    int **array = (int **)malloc(n*sizeof(int *));
    for (int i=0; i<n; i++) {
        array[i] = &(data[i*m]);
    } 
    return array;
}

void free2d_int(int **array) {
    free(&(array[0][0]));
    free(array);
}

char ***alloc3d_char(int n, int m, int l) {
    char *data = malloc(n*m*l*sizeof(char));
    if (!data) return NULL;

    char ***array = (char ***)malloc(n*sizeof(char **));
    for (int i=0; i<n; i++) {
        array[i] = (char **)malloc(m*sizeof(char *));
        for (int j=0; j<m; j++) {
            array[i][j] = &(data[(i*m+j)*l]);
        }
    } 
    return array;
}

void free3d_char(char ***array, int n) {
    free(&(array[0][0][0]));
    for (int i=0; i<n; i++) {
        free(array[i]);
    } 
    free(array);
}

