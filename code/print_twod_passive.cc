#include <print_sol.h>

//===========================================
//
// Output 2D passive scalar solution
//
//===========================================
#ifdef TWOD
#ifdef PASSIVE
void print_dg(const int N_s, const int N_E, scalar* U, const simpleMesh m, const int elem_type, const int step, const double time, const bool append){

  fullMatrix<scalar> Rho(N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> Uy(N_s, N_E);
  fullMatrix<scalar> P(N_s, N_E);
  scalar gamma = constants::GLOBAL_GAMMA;
  fullMatrix<scalar> PhiC(N_s, N_E)  ;
  fullMatrix<scalar> PhiNC(N_s, N_E)  ;

  // separate the fields
  scalar rho = 0;
  scalar ux  = 0;
  scalar uy  = 0;
  scalar et  = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Seperate the fields
      rho = U[(e*N_F+0)*N_s+i];
      ux  = U[(e*N_F+1)*N_s+i]/rho;
      uy  = U[(e*N_F+2)*N_s+i]/rho;
      et  = U[(e*N_F+3)*N_s+i];
      
      // Check for NaN error
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      Rho(i,e) = rho;
      Ux (i,e) = ux;
      Uy (i,e) = uy;
      PhiC  (i,e) = U[(e*N_F+4)*N_s+i]/rho;
      PhiNC (i,e) = U[(e*N_F+5)*N_s+i];
      P     (i,e) = (gamma-1)*(et - 0.5*(ux*ux+uy*uy)*rho);
    }
  }

  // print to the output file
  m.writeSolution(  Rho,  elem_type,  "rho.pos",  "Rho", step, time, append);
  m.writeSolution(   Ux,  elem_type,   "ux.pos",   "Ux", step, time, append);
  m.writeSolution(   Uy,  elem_type,   "uy.pos",   "Uy", step, time, append);
  m.writeSolution( PhiC,  elem_type, "phic.pos", "PhiC", step, time, append);
  m.writeSolution(PhiNC,  elem_type,"phinc.pos","PhiNC", step, time, append);
  m.writeSolution(    P,  elem_type,    "p.pos",    "P", step, time, append);
}
#endif
#endif
