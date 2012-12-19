#ifndef BASIC_FLUXES_H
#define BASIC_FLUXES_H
#include <scalar_def.h>
#include <math.h>


//
// Generic fluxes
//
arch_device scalar flux_ab(scalar rho, scalar u){return rho*u;}
arch_device scalar flux_apb(scalar a, scalar b){return a+b;}
arch_device scalar flux_ab2pc(scalar rho, scalar u, scalar p){return rho*u*u+p;} 
arch_device scalar flux_abc(scalar rho, scalar u, scalar a) {return rho*u*a;}
arch_device scalar fhll(scalar UL, scalar SL, scalar FL, scalar UR, scalar SR, scalar FR){
  return (SR*FL-SL*FR+SL*SR*(UR-UL))/fabs(SR-SL);
}


#endif 
