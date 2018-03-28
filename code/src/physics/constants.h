/*!
  \file constants.h
  \brief Global constants
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  Constants such as gamma and gravity
*/
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "scalar_def.h"

namespace constants
{
  const scalar GLOBAL_PI = 3.14159265358979; 

  //calorically perfect gas stuff:
  const scalar GLOBAL_GAMMA = 1.4;
  const scalar GLOBAL_RGAS = 287.15;
  const scalar GLOBAL_PRAN = 0.71; //My hal code uses 4*rsh/(9*rsh-5)=0.737 for rsh=1.4
  const scalar GLOBAL_CPGAS = (GLOBAL_GAMMA*GLOBAL_RGAS) / (GLOBAL_GAMMA-1.0);

  //Sutherland's Law parameters: For the constant-viscosity model, the code will use MEWREF
  //const scalar GLOBAL_MEWREF = 0.00001716; //reference viscosity level
  //const scalar GLOBAL_MEWREF = 1000; //reference viscosity level
  //mew_KushJet = 0.00010635;
  //mew_TGV = 0.01;
  const scalar GLOBAL_MEWREF = 0.000625; //reference viscosity level (ONLY for Navier-Stokes)
  const scalar GLOBAL_TREF = 273.15; //reference viscosity temperature
  const scalar GLOBAL_CVIS = 110.4; //Additional Sutherland parameter

  //Some stuff for scalar advection-diffusion: good for L=0.1 square
  /*
  const scalar GLOBAL_KLIN = 0.1; //diffusion coefficient (ONLY for scalar AD)
  const scalar GLOBAL_VX = 100.0; //scalar advection-diffusion propagation velocity
  const scalar GLOBAL_VY = 0.0;
  const scalar GLOBAL_VZ = 0.0;
  */
  //Some stuff for scalar advection-diffusion: good for 2pi cube
  const scalar GLOBAL_KLIN = 0.000625; //diffusion coefficient (ONLY for scalar AD)
  const scalar GLOBAL_VX = GLOBAL_PI; //1.0 //scalar advection-diffusion propagation velocity
  const scalar GLOBAL_VY = GLOBAL_PI; //0.5
  const scalar GLOBAL_VZ = GLOBAL_PI; //-0.5


  //Marc's original contrubution
  // I don't like the way this is done but I have nothing better for
  // now. It seems that it will initialize to zero by default.
  extern scalar GLOBAL_GX; // gravity x-dir
  extern scalar GLOBAL_GY; // gravity y-dir
  
} // namespace constants

#endif
