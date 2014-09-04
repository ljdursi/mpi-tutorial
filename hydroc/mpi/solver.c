#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "solver.h"
#include "allocmultid.h"

#define MAX(x,y) ((y) > (x) ? (y) : (x))

const float eosgamma=5./3.;

void bufferGuardcells(float ***u, const int nx, const int ny, 
                      const char xydim, 
                      MPI_Comm cartcomm, int left, int right) {
    int datasize;
    float *sendbuffer, *recvbuffer;
    int i,j,v,count;
    MPI_Status rstatus;
    const int lefttag = 1;
    const int righttag = 2;

    datasize = (ny-2*NGUARD)*NGUARD*NVARS;
    sendbuffer = (float *)malloc(datasize*sizeof(float));
    recvbuffer = (float *)malloc(datasize*sizeof(float));
    if (!sendbuffer || !recvbuffer) {
        fprintf(stderr,"Error mallocing send/recv buffer(!?)\n");
        return;
    }

    /* rightward */
    /* copy data into a buffer */

    /* MPI_Sendrecv(sendbuffer, datasize,... ?
                    recvbuffer, datasize,... ? */
    /* copy data out of a buffer */

    /* leftward */
    /* copy data into a buffer */

    /* MPI_Sendrecv(sendbuffer, datasize,... ?
                    recvbuffer, datasize,... ? */
    /* copy data out of a buffer */


    free(sendbuffer);
    free(recvbuffer);
}

void vectorGuardcells(float ***u, const int nx, const int ny, const char xydim,
                      MPI_Comm cartcomm, const int left, const int right) {
    const int lefttag=1, righttag=2;
    int ierr;
    MPI_Status rstatus;
    MPI_Datatype xbctype;

    /* ierr = MPI_Type_vector(...? */
    MPI_Type_commit(&xbctype);

    /* ierr = MPI_Sendrecv(...?  cartcomm, &rstatus); */

    /* ierr = MPI_Sendrecv(...?  cartcomm, &rstatus); */

    MPI_Type_free(&xbctype);
}

void periodicBCs(float ***u, const int nx, const int ny, const char xydim) {

    if (xydim == 'y') {
        int n = ny;
        for (int j=0; j<NGUARD; j++) {
            for (int i=0; i<nx; i++) { 
                for (int var=0 ; var<NVARS; var++) {
                    u[j][i][var] = u[n-2*NGUARD+j][i][var];
                    u[n-NGUARD+j][i][var] = u[j+NGUARD][i][var];
                }
            } 
        }
    } else if (xydim == 'x') {
        int n = nx;
        for (int j=0; j<ny; j++) { 
            for (int i=0; i<NGUARD; i++) { 
                for (int var=0 ; var<NVARS; var++) {
                    u[j][i][var] = u[j][n-2*NGUARD+i][var];
                    u[j][n-NGUARD+i][var] = u[j][i+NGUARD][var];
                }
            }
        } 
    }
}


float cfl(float ***u, const int nx, const int ny) {
   float c;
   float globalc;

   c = 0.;
   
   for (int j=0; j<ny; j++) {
       for (int i=0; i<nx; i++)  {
            float vx = fabsf(u[j][i][IMOMX]/u[j][i][IDENS]);
            float vy = fabsf(u[j][i][IMOMY]/u[j][i][IDENS]);
            float v = vx;
            if (vy > vx) v=vy;
            float p = (u[j][i][IENER] - 0.5*(vx*vx+vy*vy)*u[j][i][IDENS])*(eosgamma-1.);
            float cloc = v + sqrtf(fabsf(p*eosgamma/u[j][i][IDENS]));
            c = MAX(c,cloc);
        }
   }

   /* find global max of all cs */
   MPI_Allreduce(&c, &globalc, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   return (1./globalc);
}

void initialconditions(float ***u,const int nx, const int ny, const int n, const int start) {
  int i, j;
  const float backgrounddens=1., projdens=50., projvel=4., p0=1.;
  float r;
  float x,y;
  
 
  for (j=0;j<ny;j++) {
      y = ((j-NGUARD)*1.-n/2.);
      for (i=0;i<nx;i++) {
          x = ((i-NGUARD+start)*1.-n/2.);
          r =sqrtf(x*x+y*y);

          if (r < 0.1*sqrtf(n*n*1.+n*n*1.)) {
            u[j][i][IDENS] = projdens;
            u[j][i][IMOMX] = projvel*projdens;
            u[j][i][IMOMY] = 0.;
            u[j][i][IENER] = 0.5*(projdens*projvel*projvel)+p0/(eosgamma-1.);
          } else {
            u[j][i][IDENS] = backgrounddens;
            u[j][i][IMOMX] = 0.;
            u[j][i][IMOMY] = 0.;
            u[j][i][IENER] = p0/(eosgamma-1.);
          }
      }
  }
}

void hydroflux(float **v, float *rc, float **u, const int n) {
  float c = 0.;
  
  for (int i=0; i<n; i++) {
    float vx = u[i][IMOMX]/u[i][IDENS];
    float p  = u[i][IENER] - 0.5*(u[i][IMOMX]*u[i][IMOMX]+u[i][IMOMY]*u[i][IMOMY])
                                   /u[i][IDENS];
    p *= (eosgamma-1.);
  
    v[i][IDENS] = u[i][IMOMX];
    v[i][IMOMX] = u[i][IMOMX]*vx + p;
    v[i][IMOMY] = u[i][IMOMY]*vx;
    v[i][IENER] = (u[i][IENER]+p)*vx;

    float newc = (fabs(vx)+sqrt(fabs(eosgamma*p/u[i][IDENS])));
    c = MAX((c),newc);
  }

  if (c > 0.) {
      for (int i=0; i<n; i++) {
          for (int var=0; var<NVARS; var++) {
            v[i][var] /= c;
          }
      }
  }

  *rc = c;
}


void tvd1d(float **u, const int n, const float dt) {
    float **v, **u1, **wr, **wl;
    float **dfl, **dfr, **flux;
    float c;
    const float smallp=1.e-3;

    v  = alloc2d_float(n, NVARS);
    
    hydroflux(v, &c, u, n);

    wl = alloc2d_float(n, NVARS);
    wr = alloc2d_float(n, NVARS);
    for (int i=0; i<n; i++) {
        for (int var=0; var<NVARS; var++) {
            wr[i][var] = u[i][var] + v[i][var];
            wl[i][var] = u[i][var] - v[i][var];
        }
    }

    flux = alloc2d_float(n, NVARS);
    for (int i=0; i<n; i++) {
        int irshift = i+1;
        if (irshift >= n) irshift -= n;
        for (int var=0; var<NVARS; var++) {
            flux[i][var] = c*(wr[i][var] - wl[irshift][var])/2.;
        }
    }
    
    u1 = alloc2d_float(n, NVARS);
    for (int i=0; i<n; i++) {
        int ilshift = i-1;
        if (ilshift < 0) ilshift += n;
        for (int var=0; var<NVARS; var++) {
            u1[i][var] = u[i][var] - (flux[i][var] - flux[ilshift][var])*dt/2.;
        }
    }
    
    hydroflux(v, &c, u1, n);
    for (int i=0; i<n; i++) {
        for (int var=0; var<NVARS; var++) {
            wr[i][var] = u1[i][var] + v[i][var];
            wl[i][var] = u1[i][var] - v[i][var];
        }
    }

    dfr = alloc2d_float(n, NVARS);
    for (int i=0; i<n; i++) {
        int irshift = i+1;
        if (irshift >= n) irshift -= n;
        int ilshift = i-1;
        if (ilshift < 0) ilshift += n;
        for (int var=0; var<NVARS; var++) {
            dfr[i][var] = 0.;
            float dfrp = c*(wr[irshift][var] - wr[i][var])/2.; 
            float dfrm = c*(wr[i][var] - wr[ilshift][var])/2.; 
            if (dfrp*dfrm > 0.) 
                dfr[i][var] = (2.*(dfrp*dfrm)/(dfrp+dfrm));
        }
    }
    
    dfl = alloc2d_float(n, NVARS);
    for (int i=0; i<n; i++) {
        int irshift = i+1;
        if (irshift >= n) irshift -= n;
        int irrshift = i+2;
        if (irrshift >= n) irrshift -= n;
        for (int var=0; var<NVARS; var++) {
            dfl[i][var] = 0.;
            float dflp = c*(wl[irshift][var] - wl[irrshift][var])/2.; 
            float dflm = c*(wl[i][var] - wl[irshift][var])/2.; 
            if (dflp*dflm > 0.) 
                dfl[i][var] = (2.*(dflp*dflm)/(dflp+dflm));
        }
    }
    
    for (int i=0; i<n; i++) {
        int irshift = i+1;
        if (irshift >= n) irshift -= n;
        for (int var=0; var<NVARS; var++) {
            flux[i][var] = (c*(wr[i][var] - wl[irshift][var]) 
                             + (dfr[i][var]-dfl[i][var]))/2.;
        }
    }
    

    for (int i=NGUARD; i<n-NGUARD; i++) {
        int ilshift = i-1;
        if (ilshift < 0) ilshift += n;
        for (int var=0; var<NVARS; var++) {
            u[i][var] -= (flux[i][var] - flux[ilshift][var])*dt;
        }
        float rho=u[i][IDENS];
        float vx=u[i][IMOMX]/rho;
        float vy=u[i][IMOMY]/rho;
        if ( (u[i][IENER] - 0.5*(vx*vx+vy*vy)*rho)*(eosgamma-1.) < smallp ) 
            u[i][IENER] = smallp/(eosgamma-1.) +  0.5*(vx*vx+vy*vy)*rho ;
    }
    
    free2d_float(v);
    free2d_float(u1);
    free2d_float(wr);
    free2d_float(wl);
    free2d_float(dfl);
    free2d_float(dfr);
    free2d_float(flux);
    return;
}

void xytranspose(float ***ut, float ***u, const int nx, const int ny) {
    int i,j;
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            ut[i][j][IDENS] = u[j][i][IDENS];
            ut[i][j][IENER] = u[j][i][IENER];
            ut[i][j][IMOMX] = u[j][i][IMOMY];
            ut[i][j][IMOMY] = u[j][i][IMOMX];
        }
    }
}


void xsweep(float ***u, const int nx, const int ny, const float dt){
  for (int j=0; j<ny; j++) {
     tvd1d(u[j],nx,dt);
  }
}

void timestep(float ***u, const int nx, const int ny, float *dt, 
              MPI_Comm cartcomm, int left, int right) {
    float ***ut;

    ut = alloc3d_float(nx, ny, NVARS);
    *dt=0.25*cfl(u,nx,ny);

    /* the x sweep */
    periodicBCs(u,nx,ny,'x');
    //should be bufferGuardcells(u,nx,ny,'x',cartcomm,left,right);
    xsweep(u,nx,ny,*dt);

    /* the y sweeps */
    xytranspose(ut,u,nx,ny);
    periodicBCs(ut,ny,nx,'x');
    xsweep(ut,ny,nx,*dt);
    periodicBCs(ut,ny,nx,'x');
    xsweep(ut,ny,nx,*dt);

    /* 2nd x sweep */
    xytranspose(u,ut,ny,nx);
    periodicBCs(u,nx,ny,'x');
    //should be bufferGuardcells(u,nx,ny,'x',cartcomm,left,right);
    xsweep(u,nx,ny,*dt);

    free3d_float(ut,nx);
}

