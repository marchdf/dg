#ifndef LIMITER_FCT_H
#define LIMITER_FCT_H
#include <scalar_def.h>
#include <math.h>



inline int signum(scalar val){return val>0? 1 : (val<0? -1 : 0);}
inline scalar minabs(scalar* c, int n);
scalar minmod  (scalar* c, int n);             // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
scalar minmod2 (scalar* c, int n);             // eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
scalar cminmod (scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
scalar cminmod2(scalar* c, int n, scalar eps); // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
int factorial(int n);

#endif
