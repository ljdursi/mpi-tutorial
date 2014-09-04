#include <math.h>
#include <stdlib.h>
#include "solver.h"
#include "allocmultid.h"

const float eosgamma=5./3.;

#define MAX(x,y) ((y) > (x) ? (y) : (x))

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
   float c = 0.;
   
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

   return (1./c);
}

void initialconditions(float ***u,const int nx, const int ny) {
  int i, j;
  const float backgrounddens=1., projdens=50., projvel=4., p0=1.;
  float r;
  int n = MAX(nx,ny);
  float x,y;
 
  for (j=0;j<ny;j++) {
      y = ((j-NGUARD)*1.-n/2.);
      for (i=0;i<ny;i++) {
          x = ((i-NGUARD)*1.-n/2.);
          r =sqrt(x*x+y*y);

          if (r < 0.1*sqrt(nx*nx*1.+ny*ny*1.)) {
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
  int j;

  for (j=0; j<ny; j++) {
     tvd1d(u[j],nx,dt);
  }
}

void timestep(float ***u, const int nx, const int ny, float *dt) {
    float ***ut;

    ut = alloc3d_float(ny, nx, NVARS);
    *dt=0.25*cfl(u,nx,ny);

    /* the x sweep */
    periodicBCs(u,nx,ny,'x');
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
    xsweep(u,nx,ny,*dt);

    free3d_float(ut,ny);
}

