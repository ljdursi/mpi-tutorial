#include <stdlib.h>
#include <mpi.h>
#include "solver.h"

#ifdef PGPLOT
#include <cpgplot.h>
#endif

#define NCONTOURS 5

int plot(float ***u, const int nx, const int ny, const int start) {

#ifdef PGPLOT
	const int ng=NGUARD;
	double xl, xr, yl, yr, dx, dy;
	float  denscontours[NCONTOURS], prescontours[NCONTOURS];
	float *dens, *pres;
	double maxd=-1e19, mind=+1e-19;
	double maxp=-1e19, minp=+1e-19;
    double globmaxd, globmind, globmaxp, globminp;
	static int called = 0;
	int i, j, count;	
	float tr[6] = {0.,0.,0.,0.,0.,0.};
	
	xl = start*1.; xr = start+nx-2*ng-1;
	yl = 0.; yr = ny-2*ng-1;
    dx = 1.; dy = 1.;

	dens = (float *)malloc((nx-2*ng)*(ny-2*ng)*sizeof(float));
	pres = (float *)malloc((nx-2*ng)*(ny-2*ng)*sizeof(float));
	count = 0;
	for (j=ny-ng-1; j>=ng; j--) {
		for (i=ng; i<nx-ng; i++) {
			dens[count] = u[j][i][IDENS];
            float vx = u[j][i][IMOMX]/dens[count];
            float vy = u[j][i][IMOMY]/dens[count];
			pres[count] = (u[j][i][IENER] - 0.5*(vx*vx+vy*vy)*dens[count]);
			if (dens[count] > maxd) maxd = dens[count];
			if (dens[count] < mind) mind = dens[count];
			if (pres[count] > maxp) maxp = pres[count];
			if (pres[count] < minp) minp = pres[count];
			count++;
		}
	}
	
    MPI_Allreduce(&mind, &globmind,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&minp, &globminp,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&maxd, &globmaxd,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&maxp, &globmaxp,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for (i=0; i<NCONTOURS; i++) {
		denscontours[i] = globmind + (i+1)*(globmaxd-globmind)/(NCONTOURS+1);
		prescontours[i] = globminp + (i+1)*(globmaxp-globminp)/(NCONTOURS+1);
	}

	tr[0] = xl; tr[3] = yl;
	tr[1] = dx; tr[5] = dy;

	if (!called) {
		cpgbeg(0, "/XWINDOW", 1, 1);
		called = 1;
		cpgask(0);
	}
	cpgenv(xl, xr, yl, yr, 1, 1);

	cpgsci(2);
	cpgcont(dens, nx-2*ng, ny-2*ng, 1, nx-2*ng, 1, ny-2*ng, denscontours, NCONTOURS, tr);
	if (minp != maxp) {
		cpgsci(3);
	    cpgcont(pres, nx-2*ng, ny-2*ng, 1, nx-2*ng, 1, ny-2*ng, prescontours, NCONTOURS, tr);
	}
	cpgsci(1);
	
	free(dens); free(pres);
#endif
	return 0;
}

int closeplot(){
#ifdef PGPLOT
	cpgend();
#endif
	return 0;
}
