#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    const int nx=1500;
    float *dat;
    int i;
    float datamin, datamax, datamean;
    float minmeanmax[3];
    float globminmeanmax[3];
    int ierr;
    int rank, size;
    int tag=1;
    int masterproc=0;
    MPI_Status status;


    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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

    minmeanmax[0] = datamin;
    minmeanmax[2] = datamax;
    minmeanmax[1] = datamean;

    if (rank != masterproc) {
       ierr = MPI_Ssend( /* ... ? ... */ );
    } else {
        globminmeanmax[0] = datamin;
        globminmeanmax[2] = datamax;
        globminmeanmax[1] = datamean;
        for (i=1;i<size;i++) {
            ierr = MPI_Recv( /* ... ? ... */);

            globminmeanmax[1] += minmeanmax[1];

            if (minmeanmax[0] < globminmeanmax[0])
                globminmeanmax[0] = minmeanmax[0];

            if (minmeanmax[2] > globminmeanmax[2])
                globminmeanmax[2] = minmeanmax[2];

        }
        globminmeanmax[1] /= size;
        printf("Min/mean/max = %f,%f,%f\n", globminmeanmax[0],
               globminmeanmax[1],globminmeanmax[2]);
    }
 
    ierr = MPI_Finalize();

    return 0;
}
