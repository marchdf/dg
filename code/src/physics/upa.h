/*!
  \file upa.h
  \brief Find |u+a|
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  These functions find |u+a| for a given rho, u, E, gamma.  For 2D,
  they return the max(|u+a|,|v+a|);
*/

#ifndef UPA_H
#define UPA_H
#include "scalar_def.h"
#include <math.h>
#include <stdio.h>

#ifdef SCALARAD

  //PEJ 05/29/2017:
  //in this case, the maximum physics-based eigenvalue
  //corresponds to max(Vx,Vy,Vz); if we are getting technical,
  //max Delta t depends on the interaction between
  //these three speeds and the element geometry,
  //but for starters, I'm going to take the square
  //root of the squares
#ifdef ONED 
  arch_device scalar oned_scalarad_upa(scalar rho,
				       scalar Vx)
  {
    return fabs(Vx);
  }
#endif //end if for 1D 
  
#ifdef TWOD
  arch_device scalar twod_scalarad_upa(scalar rho,
				       scalar Vx,
				       scalar Vy) 
  {
    return sqrt(Vx*Vx + Vy*Vy);
  }
#endif //end if for 2D
  
#ifdef THREED
arch_device scalar threed_scalarad_upa(scalar rho,
				scalar Vx,
				       scalar Vy,
				       scalar Vz) 
{
  return sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
}
#endif //end if for 3D

#endif //endif for SCALARAD. Hopefully, nothing else here will activate when I set SCALARAD

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

#elif SINGLEFLUID //means we are singlefluid, not passive

#ifdef ONED
arch_device scalar oned_singlefluid_upa(scalar rho,
					scalar u,
					scalar E,
					scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*u*u);
  return fabs(u)+sqrt(gamma*p/rho);
}
#elif TWOD
arch_device scalar twod_singlefluid_upa(scalar rho,
					scalar u,
					scalar v,
					scalar E,
					scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*(u*u+v*v));
  scalar a = sqrt(gamma*p/rho);
  return MAX(fabs(u)+a,fabs(v)+a);
}		
#elif THREED
arch_device scalar threed_singlefluid_upa(scalar rho,
					  scalar u,
					  scalar v,
					  scalar w,
					  scalar E,
					  scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*(u*u + v*v + w*w));
  scalar a = sqrt(gamma*p/rho);
  return MAX(MAX(fabs(u)+a , fabs(v)+a), fabs(w)+a);
}		  
#endif // dimensions

#elif RADSINGLEFLUID //means we are singlefluid with Reisner AD

#ifdef ONED
arch_device scalar oned_radsinglefluid_upa(scalar rho,
					scalar u,
					scalar E,
					scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*u*u);
  return fabs(u)+sqrt(gamma*p/rho);
}
#elif TWOD
arch_device scalar twod_radsinglefluid_upa(scalar rho,
					scalar u,
					scalar v,
					scalar E,
					scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*(u*u+v*v));
  scalar a = sqrt(gamma*p/rho);
  return MAX(fabs(u)+a,fabs(v)+a);
}		
#elif THREED
arch_device scalar threed_radsinglefluid_upa(scalar rho,
					  scalar u,
					  scalar v,
					  scalar w,
					  scalar E,
					  scalar gamma){
  scalar p = (gamma-1)*(E-0.5*rho*(u*u + v*v + w*w));
  scalar a = sqrt(gamma*p/rho);
  return MAX(MAX(fabs(u)+a , fabs(v)+a), fabs(w)+a);
}		  
#endif // dimensions

#elif MULTIFLUID //we are multifluid, not singlefluid or passive

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

#elif STIFFENED //yet another physics type: not passive, singlefluid, or multifluid

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
