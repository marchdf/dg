#include <limiter_fct.h>

inline scalar minabs(scalar* c, int n){
  scalar minabs = fabs(c[0]);
  for(int i=1;i<n;i++) if (minabs>fabs(c[i])) minabs = fabs(c[i]);
  return minabs;
}

scalar minmod(scalar* c, int n){
  // Generalized minmod function
  // eq 2.19 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  int sign = signum(c[0]);
  for(int i=1; i<n; i++){
    if (sign!=signum(c[i])) return 0;
  }
  return sign*minabs(c,n);
}

scalar minmod2(scalar* c, int n){
  // eq 2.20 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  scalar min = c[0];
  for(int i=1; i<n; i++) if(fabs(c[i])<fabs(min)) min = c[i];
  return min;
}

scalar cminmod(scalar* c, int n, scalar eps){
  // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  // using minmod
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = minmod(cc,2);
  delete[] cc;
  return m;    
}

scalar cminmod2(scalar* c, int n, scalar eps){
  // eq 2.21 of "Hierarchical reconstruction for discontinuous Galerkin methods..."
  // using minmod2
  scalar* cc = new scalar[2];
  scalar sum = 0;
  for(int i=0;i<n;i++) sum += c[i];
  cc[0] =(1+eps)*minmod(c,n);
  cc[1] =(scalar)sum/n;
  scalar m = minmod2(cc,2);
  delete[] cc;
  return m;    
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
