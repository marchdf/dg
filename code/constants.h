/*!
  \file constants.h
  \brief Global constants
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
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
