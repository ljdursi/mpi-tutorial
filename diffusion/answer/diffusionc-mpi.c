#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include <mpi.h>

int main(int argc, char **argv) {
    /* simulation parameters */
    const int totpoints=1000;
    int locpoints;
    const float xleft = -12., xright = +12.;
    float locxleft, locxright;
    const float kappa = 1.;

    const int nsteps=100000;
    const int plotsteps=50;

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
    int rank, size;
    int start,end;
    int left, right;
    int lefttag=1, righttag=2;
    MPI_Status status;

    /* MPI Initialization */
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    locpoints = totpoints/size;
    start = rank*locpoints;
    end   = (rank+1)*locpoints - 1;
    if (rank == size-1) 
        end = totpoints-1;
    locpoints = end-start+1;

    left = rank-1;
    if (left < 0) left = MPI_PROC_NULL;
    right= rank+1;
    if (right >= size) right = MPI_PROC_NULL;

    /* set parameters */

    dx = (xright-xleft)/(totpoints-1);
    dt = dx*dx * kappa/10.;

    locxleft = xleft + start*dx;
    locxright = xleft + end*dx;
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


    /* evolve */

    for (step=0; step < nsteps; step++) {
        /* boundary conditions: keep endpoint temperatures fixed. */

        temperature[old][0] = fixedlefttemp;
        temperature[old][locpoints+1] = fixedrighttemp;

        /* update internal boundaries */

        /* send data rightwards */
        MPI_Sendrecv(&(temperature[old][locpoints]), 1, MPI_FLOAT, right, righttag, 
                     &(temperature[old][0]), 1, MPI_FLOAT, left,  righttag, MPI_COMM_WORLD, &status);
        
        /* send data leftwards */
        MPI_Sendrecv(&(temperature[old][1]), 1, MPI_FLOAT, left, lefttag, 
                     &(temperature[old][locpoints+1]), 1, MPI_FLOAT, right,  lefttag, MPI_COMM_WORLD, &status);

        for (i=1; i<locpoints+1; i++) {
            temperature[new][i] = temperature[old][i] + dt*kappa/(dx*dx) *
                (temperature[old][i+1] - 2.*temperature[old][i] + 
                 temperature[old][i-1]) ;
        }


        time += dt;

        if (step % plotsteps == 0) {
            cpgbbuf();
            cpgeras();
            cpgsls(2); cpgslw(12); cpgsci(red);
            cpgline(locpoints+2, x, temperature[new]);
        }

        /* update correct solution */

        sigma = sqrt(2.*kappa*time + sigmao*sigmao);
        a = ao*sigmao/sigma;
        for (i=0; i<locpoints+2; i++) {
            theory[i] = a*exp(-(x[i]*x[i]) / (2.*sigma*sigma));
        }

        if (step % plotsteps == 0) {
            cpgsls(1); cpgslw(6); cpgsci(white);
            cpgline(locpoints+2, x, theory);
            cpgebuf();
        }

        error = 0.;
        for (i=1;i<locpoints+1;i++) {
            error += (theory[i] - temperature[new][i])*(theory[i] - temperature[new][i]);
        }
        error = sqrt(error);

        if (rank == 0) 
            printf("Step = %d, Time = %g, Error = %g\n", step, time, error);

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
