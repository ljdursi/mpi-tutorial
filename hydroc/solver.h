#ifndef _SOLVER_H

#define _SOLVER_H 1

#define NVARS 4
#define NGUARD 2

#define IDENS 0
#define IMOMY 1
#define IMOMX 2
#define IENER 3

void timestep(float ***u, const int nx, const int ny, float *dt);
void initialconditions(float ***u,const int nx, const int ny);

#endif
