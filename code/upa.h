/*!
  \file upa.h
  \brief Find |u+a|
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \section Description
  These functions find |u+a| for a given rho, u, E, gamma.  For 2D,
  they return the max(|u+a|,|v+a|);
*/

#ifndef UPA_H
#define UPA_H
#include <scalar_def.h>
#include <math.h>
#include <stdio.h>

#ifdef PASSIVE
#ifdef ONED
arch_device scalar oned_passive_upa(scalar rho,
				    scalar u,
				    scalar E,
				    scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*u*u);
  return fabs(u)+sqrt(gamma*p/rho);
}
#elif TWOD
arch_device scalar twod_passive_upa(scalar rho,
				    scalar u,
				    scalar v,
				    scalar E,
				    scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*(u*u+v*v));
  scalar a = sqrt(gamma*p/rho);
  return MAX(fabs(u)+a,fabs(v)+a);
}
				  
#endif // dimensions
#elif MULTIFLUID
#ifdef ONED
arch_device scalar oned_multifluid_upa(scalar rho,
				       scalar u,
				       scalar E,
				       scalar alpha){
#ifdef GAMCONS
  alpha = alpha/rho;
#endif
  scalar gamma = 1.0+1.0/alpha;
  scalar p = (gamma-1)*(E-0.5*rho*u*u);
  return fabs(u)+sqrt(gamma*p/rho);
}
#elif TWOD
arch_device scalar twod_multifluid_upa(scalar rho,
				       scalar u,
				       scalar v,
				       scalar E,
				       scalar alpha){
#ifdef GAMCONS
  alpha = alpha/rho;
#endif
  scalar gamma = 1.0+1.0/alpha;
  scalar p = (gamma-1)*(E-0.5*rho*(u*u+v*v));
  scalar a = sqrt(gamma*p/rho);
  return MAX(fabs(u)+a,fabs(v)+a);
}

#endif //dimensions

#elif STIFFENED
#ifdef ONED
arch_device scalar oned_stiffened_upa(scalar rho,
				      scalar u,
				      scalar E,
				      scalar alpha,
				      scalar beta){
  scalar gamma = 1.0+1.0/alpha;
  scalar pinf = (1-1.0/gamma)*beta;
  scalar p = (gamma-1)*(E - beta - 0.5*rho*u*u);
  return fabs(u)+sqrt(gamma*(p+pinf)/rho);
}
#elif TWOD
arch_device scalar twod_stiffened_upa(scalar rho,
				       scalar u,
				       scalar v,
				       scalar E,
				       scalar alpha,
				       scalar beta){
  scalar gamma = 1.0+1.0/alpha;
  scalar pinf = (1-1.0/gamma)*beta;
  scalar p = (gamma-1)*(E - beta - 0.5*rho*(u*u+v*v));
  scalar a = sqrt(gamma*(p+pinf)/rho);
  return MAX(fabs(u)+a,fabs(v)+a);
}
#endif //dimensions
#endif // problem type
#endif // header file
