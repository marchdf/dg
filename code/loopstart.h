/*!
  \file loopstart.h
  \brief Start definitions fo preprocessor for loop
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  This is straight from: http://www.codeproject.com/Tips/444338/Pre-processor-Iteration
  And the code is on github: https://github.com/debdattabasu/pp_Iteration

  Right now things max out at 10 iterations but I could increase that easily
   
  You use it by doing, for example: 
  #include "loopstart.h"
  #define LOOP_END 10
  #define MACRO(x) printf("%d\n", x);
  #include "loop.h"
*/

#ifdef LOOP_START
#undef LOOP_START
#endif

#ifdef LOOP_MAX
#undef LOOP_MAX
#endif

#define LOOP_START 0
#define LOOP_MAX 10

#ifdef LOOP_END
#undef LOOP_END
#endif

#ifdef MACRO
#undef MACRO
#endif
