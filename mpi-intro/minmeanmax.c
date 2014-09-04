#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    const int nx=1500;
    float *dat;
    int i;
    float datamin, datamax, datamean;

    /*
     * generate random data
     */

    dat = (float *)malloc(nx * sizeof(float));
    srand(0);
    for (i=0;i<nx;i++) {
        dat[i] = 2*((float)rand()/RAND_MAX)-1.;
    }

    /*
     * find min/mean/max
     */ 

    datamin = 1e+19;
    datamax =-1e+19;
    datamean = 0;

    
    for (i=0;i<nx;i++) {
        if (dat[i] < datamin) datamin=dat[i];
        if (dat[i] > datamax) datamax=dat[i];
        datamean += dat[i];
    }
    datamean /= nx;
    free(dat);

    printf("Min/mean/max = %f,%f,%f\n", datamin,datamean,datamax);

    return 0;
}
