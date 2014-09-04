#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "cpgplot.h"

int main(int argc, char **argv) {
    /* simulation parameters */
    const int totpoints=10000000;
    int locpoints;
    int starti, endi;
    const float xleft = -12., xright = +12.;
    float locxleft, locxright;
    const float kappa = 1.;

    const int nsteps=1000;
    const int plotsteps=50;
    const int righttag = 1;
    const int lefttag = 2;

    /* data structures */
    float *x;
    float **temperature;
    float *theory;

    /* parameters of the original temperature distribution */
    const float ao=1., sigmao=1.;
    float a, sigma;

    float fixedlefttemp, fixedrighttemp;

    int old, new;
    int step, i;
    int red, grey,white;
    float time;
    float dt, dx;
    float error;

    int ierr;
    int size, rank;
    int leftneighbour, rightneighbour;
    MPI_Status statuses[4];
    MPI_Request request[4];
    float globalerror;
    const int printtask=0;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    locpoints = totpoints/size;

    starti = locpoints*rank;
    endi   = starti+locpoints-1;
    if (rank == size-1) {
        endi = totpoints-1;
        locpoints = endi-starti+1;
    }

    leftneighbour = rank - 1;
    rightneighbour = rank + 1;
    if (leftneighbour < 0) leftneighbour = MPI_PROC_NULL;
    if (rightneighbour > size-1) rightneighbour = MPI_PROC_NULL;
    

    /* set parameters */

    dx = (xright-xleft)/(totpoints-1);
    dt = dx*dx * kappa/10.;

    locxleft = starti*dx + xleft;
    locxleft = starti*dx + xleft;
    locxright = endi*dx + xleft;
    locxright = endi*dx + xleft;
    /*
     * allocate data, including ghost cells: old and new timestep
     * theory doesn't need ghost cells, but we include it for simplicity
     */

    theory = (float *)malloc((locpoints+2)*sizeof(float));
    x      = (float *)malloc((locpoints+2)*sizeof(float));
    temperature = (float **)malloc(2 * sizeof(float *));
    temperature[0] = (float *)malloc((locpoints+2)*sizeof(float));
    temperature[1] = (float *)malloc((locpoints+2)*sizeof(float));
    old = 0;
    new = 1;


    /* setup initial conditions */

    time = 0.;
    for (i=0; i<locpoints+2; i++) {
        x[i] = locxleft + (i-1)*dx;
        temperature[old][i] = ao*exp(-(x[i]*x[i]) / (2.*sigmao*sigmao));
        theory[i]           = ao*exp(-(x[i]*x[i]) / (2.*sigmao*sigmao));
    }
    fixedlefttemp = ao*exp(-(locxleft-dx)*(locxleft-dx) / (2.*sigmao*sigmao));
    fixedrighttemp= ao*exp(-(locxright+dx)*(locxright+dx)/(2.*sigmao*sigmao));

#ifdef USEPGPLOT
    cpgbeg(0, "/xwindow", 1, 1);
    cpgask(0);
    cpgenv(locxleft, locxright, 0., 1.5*ao, 0, 0);
    cpglab("x", "Temperature", "Diffusion Test");
    red = 2;  cpgscr(red,1.,0.,0.);
    grey = 3; cpgscr(grey,.2,.2,.2);
    white=4;cpgscr(white,1.0,1.0,1.0);
 
    cpgsls(1); cpgslw(1); cpgsci(grey);
    cpgline(locpoints+2, x, theory);
    cpgsls(2); cpgslw(3); cpgsci(red);
    cpgline(locpoints+2, x, temperature[old]);
#endif

    /* evolve */

    for (step=0; step < nsteps; step++) {
        /* boundary conditions: keep endpoint temperatures fixed. */

        temperature[old][0] = fixedlefttemp;
        temperature[old][locpoints+1] = fixedrighttemp;

        /* sending rightward boundary condition fill */
        ierr = MPI_Isend(&(temperature[old][locpoints]), 1, MPI_FLOAT, rightneighbour, righttag, MPI_COMM_WORLD, &request[0]);
        ierr = MPI_Irecv(&(temperature[old][0]), 1, MPI_FLOAT, leftneighbour, righttag, MPI_COMM_WORLD, &request[1]);
                   
        /* sending leftward leftward */
        ierr = MPI_Isend(&(temperature[old][1]), 1, MPI_FLOAT, leftneighbour, lefttag, MPI_COMM_WORLD, &request[2]);
        ierr = MPI_Irecv(&(temperature[old][locpoints+1]), 1, MPI_FLOAT, rightneighbour, lefttag, MPI_COMM_WORLD, &request[3]);

        for (i=2; i<locpoints; i++) {
            temperature[new][i] = temperature[old][i] + dt*kappa/(dx*dx) *
                         (temperature[old][i+1] - 2.*temperature[old][i] + 
                          temperature[old][i-1]) ;
        }
 
        ierr = MPI_Waitall(4, request, statuses);
        i = 1;
        temperature[new][i] = temperature[old][i] + dt*kappa/(dx*dx) *
                       (temperature[old][i+1] - 2.*temperature[old][i] + 
                        temperature[old][i-1]) ;
        i = locpoints;
        temperature[new][i] = temperature[old][i] + dt*kappa/(dx*dx) *
                       (temperature[old][i+1] - 2.*temperature[old][i] + 
                        temperature[old][i-1]) ;
    

        time += dt;

#ifdef USEPGPLOT
        if (step % plotsteps == 0) {
	  cpgbbuf();
	  //cpgenv(locxleft, locxright, 0., 1.5*ao, 0, 0);
	  cpgeras();
	  cpgsls(2); cpgslw(12); cpgsci(red);
	  cpgline(locpoints+2, x, temperature[new]);
        }
#endif
        /* update correct solution */

        sigma = sqrt(2.*kappa*time + sigmao*sigmao);
        a = ao*sigmao/sigma;
        for (i=0; i<locpoints+2; i++) {
            theory[i] = a*exp(-(x[i]*x[i]) / (2.*sigma*sigma));
        }

#ifdef USEPGPLOT
        if (step % plotsteps == 0) {
               cpgsls(1); cpgslw(6); cpgsci(white);
               cpgline(locpoints+2, x, theory);
	       cpgebuf();
        }
#endif  
        error = 0.;
        for (i=1;i<locpoints+1;i++) {
            error += (theory[i] - temperature[new][i])*(theory[i] - temperature[new][i]);
        }
        ierr = MPI_Reduce(&error, &globalerror, 1, MPI_FLOAT, MPI_SUM, printtask, MPI_COMM_WORLD);
        if (rank == printtask && (step & plotsteps == 0)) {
            globalerror = sqrt(globalerror); 
            printf("Step = %d, Time = %g, Error = %g\n", step, time, error);
        }
        old = new;
        new = 1 - old;
    }


    /*
     * free data
     */

    free(temperature[1]);
    free(temperature[0]);
    free(temperature);
    free(x);
    free(theory);
    
    ierr = MPI_Finalize();

    return 0;
}
