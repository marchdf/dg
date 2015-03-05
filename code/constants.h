/*!
  \file constants.h
  \brief Global constants
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  Constants such as gamma and gravity
*/
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <scalar_def.h>

namespace constants
{
  const scalar GLOBAL_GAMMA = 1.4;

  // I don't like the way this is done but I have nothing better for
  // now. It seems that it will initialize to zero by default.
  extern scalar GLOBAL_GX; // gravity x-dir
  extern scalar GLOBAL_GY; // gravity y-dir
  
} // namespace constants

#endif
