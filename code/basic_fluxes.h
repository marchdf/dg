/*!
  \file basic_fluxes.h  
  \brief Functions to calculate simple flux products  
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup fluxes
  \section Description
  Simple flux products often used
*/
#ifndef BASIC_FLUXES_H
#define BASIC_FLUXES_H
#include <scalar_def.h>
#include <math.h>


//
// Generic fluxes
//
arch_device inline scalar flux_ab(scalar rho, scalar u){/*!Return rho*u*/return rho*u;}
arch_device inline scalar flux_apb(scalar a, scalar b){/*!Return a+b*/return a+b;}
arch_device inline scalar flux_ab2pc(scalar rho, scalar u, scalar p){/*!Return rho*u*u+p*/return rho*u*u+p;} 
arch_device inline scalar flux_abc(scalar rho, scalar u, scalar a) {/*!Return rho*u*a*/return rho*u*a;}
arch_device inline scalar fhll(scalar UL, scalar SL, scalar FL, scalar UR, scalar SR, scalar FR){
  /*! Return HLL flux: (SR*FL-SL*FR+SL*SR*(UR-UL))/fabs(SR-SL)*/
  return (SR*FL-SL*FR+SL*SR*(UR-UL))/fabs(SR-SL);
}


#endif 
