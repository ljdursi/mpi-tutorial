#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <getopt.h>
#include "allocmultid.h"

#ifdef PGPLOT
#include <cpgplot.h>
#endif

#define NDIM 3
#define EPS 0.1

#define GRAVCONST 1.         /* sets timescale */
#define SIM_RANDOM  1        /* which ICs are we going to use? */

/*--------------------------------------------------------------------------------*/


typedef struct nbody_struct_s {
      double x[NDIM];    /*the particle positions*/
      double v[NDIM];    /*the particle velocities*/
      double f[NDIM];    /*the forces on the particles*/
      double mass;
      double PE;        /* potential energy */
} NBody;


/*--------------------------------------------------------------------------------*/

NBody *allocate_nbody(int n)
{
   NBody *mydata;
   mydata=(NBody *)malloc(n * sizeof(NBody));
   return mydata;
}

/*--------------------------------------------------------------------------------*/

void initialize_particles(NBody *data, int n, int simulation) 
{
   int i;
   double vamp = 0.1; /* KE/U ~ 5/6 vamp^2 */

   switch (simulation) {
      case SIM_RANDOM:
         i = 0;
         while (i < n) {
            double r = 0.;
            for (int j=0; j<NDIM; j++) {
               data[i].x[j] = 2.0*rand()/RAND_MAX - 1;
               r += data[i].x[j] * data[i].x[j];
          
               data[i].v[j]=vamp*(2.0*rand()/RAND_MAX - 1.);
               data[i].f[j]=0;
            } 
            if (r > 1.) continue;

            data[i].mass = 1./n;
            i++;
         }
         break;

   } 
   return;
}


/*--------------------------------------------------------------------------------*/

void calculate_forces_fastest(NBody *data, int n)
{
   for (int i=0;i<n;i++) {
      for (int j=0;j<NDIM;j++) {
         data[i].f[j]=0; 
      }
      data[i].PE = 0.;
   }

   for (int i=0;i<n;i++){
      for (int j=0;j<i;j++) {
         double rsq=EPS*EPS, dx[NDIM],forcex;
         for (int k=0;k<NDIM;k++){
            dx[k]=data[j].x[k]-data[i].x[k];
            rsq+=dx[k]*dx[k];
         }
         double ir =1./sqrt(rsq);
         rsq=ir/rsq;
         for (int k=0;k<NDIM;k++)
         {
            forcex=rsq*dx[k] * data[i].mass * data[j].mass * GRAVCONST;
            data[i].f[k] += forcex;
            data[j].f[k] -= forcex;
         }
         data[i].PE -= GRAVCONST * data[i].mass * data[j].mass * ir;
         data[j].PE -= GRAVCONST * data[i].mass * data[j].mass * ir;
      }
   }
}

/*--------------------------------------------------------------------------------*/

void calculate_energy(NBody *data, int n, double  *energy) {
   *energy = 0.;
  
   for (int i=0; i<n; i++) {
      *energy += 0.5*data[i].PE;
      double ke = 0.; 
      for (int j=0; j<NDIM; j++) 
         ke += (data[i].v[j])*(data[i].v[j]);
      ke *= 0.5 * data[i].mass;
      *energy += ke;
   }

   return;
}

/*--------------------------------------------------------------------------------*/

void calculate_forces(NBody *data, int n)
{
   calculate_forces_fastest(data,n);
   //calculate_forces_fastest_omp(data);
				  
   return;
}

/*--------------------------------------------------------------------------------*/

/* kick-drift-kick */
void nbody_step(NBody *data, int n, double dt)
{
   for (int i=0;i<n;i++)
      for (int j=0;j<NDIM;j++) {
         data[i].x[j] += data[i].v[j]*dt;
      }
   calculate_forces(data, n);
   for (int i=0;i<n;i++)
      for (int j=0;j<NDIM;j++) {
         double a = data[i].f[j] / data[i].mass;
         data[i].v[j] += a*dt;
      }

   return;
}

/*--------------------------------------------------------------------------------*/

void display_particles(NBody *data, double max, int n)
/*assume a top-down view (i.e. ignore the z-component*/
{
#ifdef PGPLOT
   int npix=256;
   float tr[6];
   static int firsttime=1;

   if (firsttime)
   {
      assert(NDIM>=2);
      cpgopen("?");
      cpgenv(1,npix,1,npix,0,1);
      firsttime=0;
   }
   tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
   float **mat=alloc2d_float(npix,npix);
   for (int i=0;i<npix;i++)
      for (int j=0;j<npix;j++)
         mat[i][j]=0;
  
   for (int i=0;i<n;i++)
   {
      int ii= npix*((data[i].x[0]+max)/(2*max));
      int jj= npix*((data[i].x[1]+max)/(2*max));
      if ((ii>=0)&&(ii<npix)&&(jj>=0)&&(jj<npix))
         mat[ii][jj]++;
   }

   cpgimag(mat[0],npix,npix,1,npix,1,npix,0,5.0,tr);
#endif
  
}


/*--------------------------------------------------------------------------------*/

int get_options(int argc, char **argv, int *npts, int *nsteps, int *dooutput, int *outevery, int *simulation, double *dt) {

   int output = *dooutput;

   const struct option long_options[] = {
      {"output"   , no_argument, &output, 1},
      {"nooutput" , no_argument, &output, 0},
      {"npts"     , required_argument, 0, 'n'},
      {"nsteps"   , required_argument, 0, 't'},
      {"simulation", required_argument, 0, 's'},
      {"outevery", required_argument, 0, 'o'},
      {"dt", required_argument, 0, 'd'},
      {"nproc", required_argument, 0, 'P'},
      {"help",      no_argument, 0, 'h'},
      {0, 0, 0, 0}};

   char c;
   int option_index;
   int tempint;
   double tempflt;

   while (1) {
      c = getopt_long(argc, argv, "n:t:s:o:P:h:e", long_options, &option_index);
      if (c == -1) break;

      switch (c) {
         case 0: if (long_options[option_index].flag != 0)
            break;

         case 'n': tempint = atoi(optarg);
            if (tempint < 1 || tempint > 100000) {
               fprintf(stderr,"%s: Cannot use number of points %s;\n  Using %d\n", argv[0], optarg, *npts);
            } else {
               *npts = tempint;
            }
            break;


         case 't': tempint = atoi(optarg);
            if (tempint < 1 ) {
               fprintf(stderr,"%s: Cannot use number of steps  %s;\n  Using %d\n", argv[0], optarg, *nsteps);
            } else {
               *nsteps  = tempint;
            }
            break;

         case 'o': tempint = atoi(optarg);
            if (tempint < 1) {
               fprintf(stderr,"%s: invalid output frequency %s;\n  Using %d\n", argv[0], optarg, *outevery);
            } else {
               *outevery = tempint;
            }
            break;

         case 's': tempint = atoi(optarg);
            if (tempint != SIM_RANDOM) {
               fprintf(stderr,"%s: invalid simulation number %s;\n  Using %d\n", argv[0], optarg, *simulation);
            } else {
               *simulation = SIM_RANDOM;
            }
            break;


         case 'd': tempflt = atof(optarg);
            if (tempflt <= 0. || tempflt > 1.) {
               fprintf(stderr,"%s: invalid timestep parameter %s;\n  Using %g\n", argv[0], optarg, *dt);
            } else {
               *dt = tempflt;
            }
            break;

         default: printf("Invalid option %s\n", optarg);

         case 'h':
            puts("Options: ");
            puts("    --npts=N      (-n N): Set the number of particles.");
            puts("    --nsteps=N    (-t N): Set the number of time steps taken.");
            puts("    --simulation=N(-s N): Set the simulation; 0 == Random ICs.");
            puts("    --outevery=N  (-o N): Set output every N timesteps.");
            puts("    --dt=X  (-d X): Timestep size dt.");
            puts("");
            return +1;

            break;
      }
   }

   *dooutput = output;
   return 0;
}


/*================================================================================*/

int main(int argc, char *argv[])
{
   int npts = 1500;
   int nsteps = 500;
   int outevery = 1;
   int simulation = SIM_RANDOM;
   int status;
   int output = 1;
   double dt = 0.01;
   double time = 0.;
   double tote = 0.;

   status = get_options(argc, argv, &npts, &nsteps, &output, &outevery, &simulation, &dt);
   if (status) return 0;

   NBody *mydata;


   mydata=allocate_nbody(npts);
   printf("mydata has %d particles in a %d-dimensional space.\n",npts,NDIM);
  
   initialize_particles(mydata, npts, simulation);
  
   calculate_forces(mydata, npts);
   calculate_energy(mydata, npts, &tote);

   for (int i=0;i<nsteps;i++)
   {
      nbody_step(mydata,npts,dt);
      calculate_energy(mydata, npts, &tote);
      time += dt;
      if (output) {
         printf("%i\t%g\t%g\t%g\n", i, dt, time, tote);
         if (!(i % outevery)) display_particles(mydata,1.3,npts);
      }
   }
  
  
   FILE *outfile=fopen("myparticles.txt","w");
   fprintf(outfile,"positions and forces are:\n");
   for (int i=0;i<npts;i++) 
   {
      for (int j=0;j<NDIM;j++)
         fprintf(outfile,"%12.6f ",mydata[i].x[j]);
      fprintf(outfile,"      ");
      for (int j=0;j<NDIM;j++)
         fprintf(outfile,"%12.6f ",mydata[i].f[j]);
      fprintf(outfile,"\n");
   }


   return 0;
}
