/*!
  \file loop.h
  \brief Preprocessor for loop... shady business
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  This is straight from: http://www.codeproject.com/Tips/444338/Pre-processor-Iteration
  And the code is on github: https://github.com/debdattabasu/pp_Iteration
  Licensed under The Code Project Open License (CPOL)
  
  Right now things max out at 10 iterations but I could increase that easily
   
  You use it by doing, for example: 
  #include "loopstart.h"
  #define LOOP_END 10
  #define MACRO(x) printf("%d\n", x);
  #include "loop.h"
*/
#if(LOOP_END > LOOP_MAX)
#error Loop Counter Too Big
#endif

#if(LOOP_END != 0)
MACRO(LOOP_START) 

#if(LOOP_START == 0)
#undef LOOP_START
#define LOOP_START 1
#elif(LOOP_START == 1)
#undef LOOP_START
#define LOOP_START 2
#elif(LOOP_START == 2)
#undef LOOP_START
#define LOOP_START 3
#elif(LOOP_START == 3)
#undef LOOP_START
#define LOOP_START 4
#elif(LOOP_START == 4)
#undef LOOP_START
#define LOOP_START 5
#elif(LOOP_START == 5)
#undef LOOP_START
#define LOOP_START 6
#elif(LOOP_START == 6)
#undef LOOP_START
#define LOOP_START 7
#elif(LOOP_START == 7)
#undef LOOP_START
#define LOOP_START 8
#elif(LOOP_START == 8)
#undef LOOP_START
#define LOOP_START 9
#elif(LOOP_START == 9)
#undef LOOP_START
#define LOOP_START 10
#endif

#if(LOOP_START < LOOP_END)
#include "loop.h"
#endif

#endif
