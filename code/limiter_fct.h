#ifndef LIMITER_FCT_H
#define LIMITER_FCT_H
#include <scalar_def.h>
#include <math.h>



inline int signum(scalar val){return val>0? 1 : (val<0? -1 : 0);}
inline scalar minabs(scalar* c, int n);
scalar minmod(scalar* c, int n);
int factorial(int n);

#endif
