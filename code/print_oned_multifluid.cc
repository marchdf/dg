#include <print_sol.h>

//===========================================
//
// Output 1D multifluid scalar solution
//
//===========================================
#ifdef ONED
#ifdef MULTIFLUID

// Used to define dynamically variables
#define Y(x) Y ##x 

void print_dg(const int N_s, const int N_E, scalar* U, const simpleMesh m, const int elem_type, const int step, const double time, const int append){

  fullMatrix<scalar> Rho(N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> P(N_s, N_E);
  fullMatrix<scalar> G(N_s, N_E);

  // separate the fields
  scalar rho = 0;
  scalar ux  = 0;
  scalar et  = 0;
  scalar gamma = 0;
  
  // Mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) fullMatrix<scalar> Y(x)(N_s,N_E); 
#include "loop.h"

  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Seperate the fields
      rho = U[(e*N_F+0)*N_s+i];
      ux  = U[(e*N_F+1)*N_s+i]/rho;
      et  = U[(e*N_F+2)*N_s+i];
#ifdef GAMCONS
      gamma = 1+rho/U[(e*N_F+3)*N_s+i];
#elif  GAMNCON
      gamma = 1+1.0/U[(e*N_F+3)*N_s+i];
#endif

      // Check for NaN error
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      Rho(i,e) = rho;
      Ux (i,e) = ux;
      G  (i,e) = gamma;
      P  (i,e) = (gamma-1)*(et - 0.5*ux*ux*rho);

      // Mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Y(x)(i,e) = U[(e*N_F+4+x)*N_s+i]/rho;
#include "loop.h"
    }
  }

 
  // print to the output file
  m.writeSolution( Rho,  elem_type, "rho.pos", "Rho", step, time, append);
  m.writeSolution(  Ux,  elem_type,  "ux.pos",  "Ux", step, time, append);
  m.writeSolution(   G,  elem_type,   "g.pos",   "G", step, time, append);
  m.writeSolution(   P,  elem_type,   "p.pos",   "P", step, time, append);

  // Mass fractions output to file
  char buffer1 [5]; char buffer2 [2];
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) sprintf(buffer1, "y%i.pos", x); sprintf(buffer2, "Y%i", x); m.writeSolution(Y(x), elem_type, buffer1, buffer2, step, time, append);
#include "loop.h"

}
#endif
#endif
