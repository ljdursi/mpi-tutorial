#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"

int main(int argc, char **argv) {
    /* simulation parameters */
    const int totpoints=1000;
    const float xleft = -12., xright = +12.;
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

    /* set parameters */

    dx = (xright-xleft)/(totpoints-1);
    dt = dx*dx * kappa/10.;

    /*
     * allocate data, including ghost cells: old and new timestep
     * theory doesn't need ghost cells, but we include it for simplicity
     */

    theory = (float *)malloc((totpoints+2)*sizeof(float));
    x      = (float *)malloc((totpoints+2)*sizeof(float));
    temperature = (float **)malloc(2 * sizeof(float *));
    temperature[0] = (float *)malloc((totpoints+2)*sizeof(float));
    temperature[1] = (float *)malloc((totpoints+2)*sizeof(float));
    old = 0;
    new = 1;

    /* setup initial conditions */

    time = 0.;
    for (i=0; i<totpoints+2; i++) {
        x[i] = xleft + (i-1+0.5)*dx;
        temperature[old][i] = ao*exp(-(x[i]*x[i]) / (2.*sigmao*sigmao));
        theory[i]           = ao*exp(-(x[i]*x[i]) / (2.*sigmao*sigmao));
    }
    fixedlefttemp = ao*exp(-(xleft-dx)*(xleft-dx) / (2.*sigmao*sigmao));
    fixedrighttemp= ao*exp(-(xright+dx)*(xright+dx)/(2.*sigmao*sigmao));

#ifdef PGPLOT
    cpgbeg(0, "/xwindow", 1, 1);
    cpgask(0);
    cpgenv(xleft, xright, 0., 1.5*ao, 0, 0);
    cpglab("x", "Temperature", "Diffusion Test");
    red = 2;  cpgscr(red,1.,0.,0.);
    grey = 3; cpgscr(grey,.2,.2,.2);
    white=4;cpgscr(white,1.0,1.0,1.0);

    cpgsls(1); cpgslw(1); cpgsci(grey);
    cpgline(totpoints+2, x, theory);
    cpgsls(2); cpgslw(3); cpgsci(red);
    cpgline(totpoints+2, x, temperature[old]);
#endif

    /* evolve */

    for (step=0; step < nsteps; step++) {
        /* boundary conditions: keep endpoint temperatures fixed. */

        temperature[old][0] = fixedlefttemp;
        temperature[old][totpoints+1] = fixedrighttemp;

        for (i=1; i<totpoints+1; i++) {
            temperature[new][i] = temperature[old][i] + dt*kappa/(dx*dx) *
                (temperature[old][i+1] - 2.*temperature[old][i] + 
                 temperature[old][i-1]) ;
        }


        time += dt;

#ifdef PGPLOT
        if (step % plotsteps == 0) {
            cpgbbuf();
            cpgeras();
            cpgsls(2); cpgslw(12); cpgsci(red);
            cpgline(totpoints+2, x, temperature[new]);
        }
#endif

        /* update correct solution */

        sigma = sqrt(2.*kappa*time + sigmao*sigmao);
        a = ao*sigmao/sigma;
        for (i=0; i<totpoints+2; i++) {
            theory[i] = a*exp(-(x[i]*x[i]) / (2.*sigma*sigma));
        }

#ifdef PGPLOT
        if (step % plotsteps == 0) {
            cpgsls(1); cpgslw(6); cpgsci(white);
            cpgline(totpoints+2, x, theory);
            cpgebuf();
        }
#endif
        error = 0.;
        for (i=1;i<totpoints+1;i++) {
            error += (theory[i] - temperature[new][i])*(theory[i] - temperature[new][i]);
        }
        error = sqrt(error);

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

    return 0;
}
