#include <print_sol.h>

//===========================================
//
// Output 2D multifluid scalar solution
//
//===========================================
#ifdef TWOD
#ifdef MULTIFLUID
void print_dg(const int N_s, const int N_E, const int N_F, scalar* U, const simpleMesh m, const int elem_type, const int step, const double time, const int append){

  fullMatrix<scalar> Rho(N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> Uy(N_s, N_E);
  fullMatrix<scalar> Et(N_s, N_E);
  fullMatrix<scalar> P(N_s, N_E);
  fullMatrix<scalar> G(N_s, N_E)  ;

  // separate the fields
  scalar rho = 0;
  scalar ux  = 0;
  scalar uy  = 0;
  scalar et  = 0;
  scalar gamma = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Seperate the fields
      rho = U[(e*N_F+0)*N_s+i];
      ux  = U[(e*N_F+1)*N_s+i]/rho;
      uy  = U[(e*N_F+2)*N_s+i]/rho;
      et  = U[(e*N_F+3)*N_s+i];
#ifdef GAMCONS
      gamma = 1+rho/U[(e*N_F+4)*N_s+i];
#elif  GAMNCON
      gamma = 1+1.0/U[(e*N_F+4)*N_s+i];
#endif

      // Check for NaN error
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      Rho(i,e) = rho;
      Ux (i,e) = ux;
      Uy (i,e) = uy;
      Et (i,e) = et;
      G  (i,e) = gamma;
      P  (i,e) = (gamma-1)*(et - 0.5*(ux*ux+uy*uy)/rho);
    }
  }

  // print to the output file
  m.writeSolution( Rho,  elem_type, "rho.pos", "Rho", step, time, append);
  m.writeSolution(  Ux,  elem_type,  "ux.pos",  "Ux", step, time, append);
  m.writeSolution(  Uy,  elem_type,  "uy.pos",  "Uy", step, time, append);
  m.writeSolution(  Et,  elem_type,  "et.pos",  "Et", step, time, append);
  m.writeSolution(   G,  elem_type,   "g.pos",   "G", step, time, append);
  m.writeSolution(   P,  elem_type,   "p.pos",   "P", step, time, append);
}
#endif
#endif
