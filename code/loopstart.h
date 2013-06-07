/* Preprocessor for loop... shady business
   This is straight from: http://www.codeproject.com/Tips/444338/Pre-processor-Iteration
   And the code is on github: https://github.com/debdattabasu/pp_Iteration

   Right now things max out at 10 iterations but I could increase that easily
   
   You use it by doing, for example: 
   #include "loopstart.h"
   #define LOOP_END 10
   #define MACRO(x) printf("%d\n", x);
   #include "loop.h"
 */

#define LOOP_START 0
#define LOOP_MAX 10

#ifdef LOOP_END
#undef LOOP_END
#endif

#ifdef MACRO
#undef MACRO
#endif
