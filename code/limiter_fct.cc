#include <limiter_fct.h>

inline scalar minabs(scalar* c, int n){
  scalar minabs = fabs(c[0]);
  for(int i=1;i<n;i++) if (minabs>fabs(c[i])) minabs = fabs(c[i]);
  return minabs;
}

scalar minmod(scalar* c, int n){
  // Generalize minmod function
  int sign = signum(c[0]);
  for(int i=1; i<n; i++){
    if (sign!=signum(c[i])) return 0;
  }
  return sign*minabs(c,n);
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
