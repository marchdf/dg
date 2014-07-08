/*!
  \file limiting.cc
  \brief Function definitions for Limiting class.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup limiting
*/
#include <limiting.h>

void Limiting::HRlimiting(COMMUNICATOR &communicator, scalar* U){
  /*!
    \brief HR limiting function
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */

  _timers.start_timer(20);
  
#ifdef ONED

  // Go from lagrange to monomial representation
  blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
  // Limit the solution according to Liu
  Lhrl1D(_N_s, _N_E, N_F, _N_N, 1, _neighbors, 0, _A, _Alim);
  // Go back to lagrange representation
  blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

#elif TWOD
  if(_cartesian){
    // Go from lagrange to monomial representation wrt x (assume communication has happened)
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts)*N_F, _N_s, 1, _Lag2MonoX, _N_s, U, _N_s, 0.0, _A, _N_s);

    // Limit the solution according to Liu (for each x slice)
    Lhrl1D(_N_s1D, _N_E, N_F, _N_N, _N_s1D, _neighbors, 0, _A, _Alim);
      
    // Go to the monomial representation wrt y
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _MonoX2MonoY, _N_s, _Alim, _N_s, 0.0, _A, _N_s);

    // Communicate the elements on different partitions if necessary
    communicator.CommunicateGhosts(N_F, _A);

    // Limit the solution according to Liu (for each y slice)
    Lhrl1D(_N_s1D, _N_E, N_F, _N_N, _N_s1D, _neighbors, 2, _A, _Alim);
      
    // Go back to lagrange
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _MonoY2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);
  }
    
  if(!_cartesian){
    /* for(int e=0; e<_N_E;e++){ */
    /*   for(int fc=0;fc<N_F;fc++){ */
    /* 	U[(e*N_F+fc)*_N_s+0] = 2; */
    /* 	U[(e*N_F+fc)*_N_s+1] = 3; */
    /* 	U[(e*N_F+fc)*_N_s+2] = 3; */
    /* 	U[(e*N_F+fc)*_N_s+3] = 3; */
    /* 	U[(e*N_F+fc)*_N_s+4] = 3; */
    /* 	U[(e*N_F+fc)*_N_s+5] = 3; */
    /*   } */
    /* } */
    /* printf("Before limiting:\n"); */
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("U(e=%i,i=%i) = %f\n",8,i,U[(8*N_F+0)*_N_s+i]); */
    /* } */

    // Go from lagrange to monomial representation
    //LChangeBasis(_L2Msize1, _L2Msize2, _N_E, _Lag2Mono, U, _A);
      
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("A(e=%i,i=%i) = %.12f\n",8,i,_A[(8*N_F+0)*_N_s+i]); */
    /* } */
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("A(e=%i,i=%i) = %.12f\n",5,i,_A[(5*N_F+0)*_N_s+i]); */
    /* } */
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("A(e=%i,i=%i) = %.12f\n",13,i,_A[(13*N_F+0)*_N_s+i]); */
    /* } */
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("A(e=%i,i=%i) = %.12f\n",1,i,_A[(1*N_F+0)*_N_s+i]); */
    /* } */
  
    // Limit the solution according to Liu
    // Lhrl2D(_L2Msize1, _N_E, _N_G, _N_N, _order, _XYZCen, _powersXYZG, _neighbors, _TaylorDxIdx, _TaylorDyIdx, _weight, _refArea, _A, _Alim);

    /* printf("After limiting:\n"); */
    /* for(int i = 0; i < _N_s; i++){ */
    /* 	printf("Alim(e=%i,i=%i) = %f\n",8,i,_Alim[(8*N_F+0)*_N_s+i]); */
    /* } */
    // Go back to lagrange representation
    //LChangeBasis(_L2Msize2, _L2Msize1, _N_E, _Mono2Lag, _Alim, U);
    /* for(int i = 0; i < _N_s; i++){ */
    /*   printf("U(e=%i,i=%i) = %f\n",8,i,U[(8*N_F+0)*_N_s+i]); */
    /* } */
  }
#endif

  _timers.stop_timer(20);
}

void Limiting::PRlimiting(COMMUNICATOR &communicator, scalar* U){
  /*!
    \brief PR limiting function (like HRlimiting but limit the primitive variables)
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */

  // Transform to primitive variables
  _timers.start_timer(20);
  LCons2Prim(_N_s, _N_E, U);
  _timers.stop_timer(20);

  // Limit primitive variables  
  HRlimiting(communicator, U);

  // Transform to conserved variables
  _timers.start_timer(20);
  LPrim2Cons(_N_s, _N_E, U);
  _timers.stop_timer(20);
}

void Limiting::M2limiting(COMMUNICATOR &communicator, scalar* U){
  /*!
    \brief Modified limiting function (see jcp2014)
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */
  
  _timers.start_timer(20);
  
#ifdef ONED
#ifdef MULTIFLUID  //===========================================================
#ifdef GAMNCON

  // Get the density field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,0,0,U,_rho); // copy from U into rho
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _rhoMono, _rhoLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);
    
  // Get the momentum field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,_N_s,0,U,_rhou); // copy from U into rhou
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _rhouMono, _rhouLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

  // Get the energy field, transform to monomial basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,2*_N_s,0,U,_E); // copy from U into E
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _E, _N_s, 0.0, _EMono, _N_s);
    
  // Get the 1/gamma-1 field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,3*_N_s,0,U,_gamma); // copy from U into gamma
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _gammaMono, _gammaLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);
    
  // Get the pressure field, transform to monomial basis, limit
  Lpressure(_N_s, _N_E, U, _pressure);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _pressureMono, _pressureLim);

  // Get the limited kinetic energy by using the limited rho and limited rhou
  // Start in Lagrange formulation then transform to monomial basis
  Lkinetic_energy1D(_N_s,_N_E,_rho,_rhou,_K); // K = 0.5*rhou^2/rho;
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _K, _N_s, 0.0, _KLim, _N_s);
    
  // Get the limited internal energy by using the limited pressure and limited 1/(gamma-1)
  Linternal_energy_multifluid(_N_s,_N_E,1,_pressureLim,_gammaLim,_rhoeLim);

  // Reconstruct the limited energy coefficients, go back to lagrange basis
  Lreconstruct_energy(_N_s,_N_E,1,_rhoeLim,_KLim,_EMono,_ELim);
  // If we did the following, it would be normal HRL: 
  //Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _EMono, _ELim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);

  // Copy the variables back into U (now limited)
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,0     ,_rho  ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,_N_s  ,_rhou ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,2*_N_s,_E    ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,3*_N_s,_gamma,U);

  // Get the mass fraction field, transform to monomial basis, limit, transform to Lagrange basis
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,(4+x)*_N_s,0,U,_Y(x)); \
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _YMono(x), _YLim(x)); \
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s); \
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,(4+x)*_N_s,_Y(x),U);
#include "loop.h"    
    
#else 
  printf("M2L is only implemented for gamncon. Exiting...\n");
  exit(1);
#endif

#elif STIFFENED  //===========================================================

  // Get the density field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,0,0,U,_rho); // copy from U into rho
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _rhoMono, _rhoLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);
    
  // Get the momentum field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,_N_s,0,U,_rhou); // copy from U into rhou
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _rhouMono, _rhouLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

  // Get the energy field, transform to monomial basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,2*_N_s,0,U,_E); // copy from U into E
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _E, _N_s, 0.0, _EMono, _N_s);
    
  // Get the 1/gamma-1 field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,3*_N_s,0,U,_gamma); // copy from U into gamma
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _gammaMono, _gammaLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);

  // Get the gamma pinf/gamma-1 field, transform to monomial basis, limit, transform to Lagrange basis
  Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,4*_N_s,0,U,_beta); // copy from U into beta
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _beta, _N_s, 0.0, _betaMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _betaMono, _betaLim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _betaLim, _N_s, 0.0, _beta, _N_s);

  // Get the pressure field, transform to monomial basis, limit
  Lpressure(_N_s, _N_E, U, _pressure);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _pressureMono, _pressureLim);

  // Get the limited kinetic energy by using the limited rho and limited rhou
  // Start in Lagrange formulation then transform to monomial basis
  Lkinetic_energy1D(_N_s,_N_E,_rho,_rhou,_K); // K = 0.5*rhou^2/rho;
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _K, _N_s, 0.0, _KLim, _N_s);
    
  // Get the limited internal energy by using the limited pressure, limited 1/(gamma-1) and limited beta
  Linternal_energy_stiffened(_N_s,_N_E,1,_pressureLim,_gammaLim,_betaLim,_rhoeLim);

  // Reconstruct the limited energy coefficients, go back to lagrange basis
  Lreconstruct_energy(_N_s,_N_E,1,_rhoeLim,_KLim,_EMono,_ELim);
  // If we did the following, it would be normal HRL: 
  //Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _EMono, _ELim);
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);
    
  // Copy the variables back into U (now limited)
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,0     ,_rho  ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,_N_s  ,_rhou ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,2*_N_s,_E    ,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,3*_N_s,_gamma,U);
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,4*_N_s,_beta ,U);

  // Get the mass fraction field, transform to monomial basis, limit, transform to Lagrange basis
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,(5+x)*_N_s,0,U,_Y(x)); \
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
  Lhrl1D(_N_s, _N_E, 1, _N_N, 1, _neighbors, 0, _YMono(x), _YLim(x)); \
  blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Mono2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s); \
  Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,(5+x)*_N_s,_Y(x),U);
#include "loop.h"    


#else 
  printf("M2L is only implemented for multifluid and stiffened. Exiting...\n");
  exit(1);
#endif // problem type
    
#elif TWOD
#ifdef MULTIFLUID  //===========================================================
#ifdef GAMNCON
  if(_cartesian){

    //
    // Limit wrt x
    // 
      
    // Get the density field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy((_N_E+_N_ghosts),_N_s,_N_s*N_F,_N_s,0,0,U,_rho); // copy from U into rho
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhoMono, _rhoLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);
      
    // Get the momentum x field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy((_N_E+_N_ghosts),_N_s,_N_s*N_F,_N_s,_N_s,0,U,_rhou); // copy from U into rhou
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhouMono, _rhouLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

    // Get the momentum y field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy((_N_E+_N_ghosts),_N_s,_N_s*N_F,_N_s,2*_N_s,0,U,_rhov); // copy from U into rhov
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _rhov, _N_s, 0.0, _rhovMono, _N_s);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhovMono, _rhovLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhovLim, _N_s, 0.0, _rhov, _N_s);
    
    // Get the 1/gamma-1 field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy((_N_E+_N_ghosts),_N_s,_N_s*N_F,_N_s,4*_N_s,0,U,_gamma); // copy from U into gamma
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _gammaMono, _gammaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);

    // Get the pressure field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lpressure(_N_s, (_N_E+_N_ghosts), U, _pressure);
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _pressureMono, _pressureLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _pressureLim, _N_s, 0.0, _pressure, _N_s);

    // Get the limited kinetic energy by using the limited rho and limited rhou and limited rhov
    // Start in Lagrange formulation then transform to monomial basis in x
    Lkinetic_energy2D(_N_s,_N_E,_rho,_rhou,_rhov,_K); // K = 0.5*(rhou^2+rhov^2)/rho;
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _K, _N_s, 0.0, _KLim, _N_s); 

    // Get the limited internal energy by using the limited pressure and limited 1/(gamma-1) in x
    Linternal_energy_multifluid(_N_s1D,_N_E,_N_s1D,_pressureLim,_gammaLim,_rhoeLim); 

    // Get the energy field, transform to monomial basis in x
    // Reconstruct the limited energy coefficients in x, back to lagrange
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,3*_N_s,0,U,_E); // copy from U into E
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _E, _N_s, 0.0, _EMono, _N_s); 
    Lreconstruct_energy(_N_s1D,_N_E,_N_s1D, _rhoeLim,_KLim,_EMono,_ELim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);

    // Get the mass fraction field, transform to monomial basis, limit in x, transform to Lagrange basis
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Lstridedcopy((_N_E+_N_ghosts),_N_s,_N_s*N_F,_N_s,(5+x)*_N_s,0,U,_Y(x)); \
    blasGemm('N','N', _N_s, (_N_E+_N_ghosts), _N_s, 1, _Lag2MonoX, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _YMono(x), _YLim(x)); \
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s);
#include "loop.h"    

    //
    // Limit wrt y
    //

    // Get the density field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
    communicator.CommunicateGhosts(1, _rhoMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhoMono, _rhoLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);

    // Get the momentum x field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
    communicator.CommunicateGhosts(1, _rhouMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhouMono, _rhouLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

    // Get the momentum y field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rhov, _N_s, 0.0, _rhovMono, _N_s);
    communicator.CommunicateGhosts(1, _rhovMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhovMono, _rhovLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhovLim, _N_s, 0.0, _rhov, _N_s);

    // Get the 1/gamma-1 field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
    communicator.CommunicateGhosts(1, _gammaMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _gammaMono, _gammaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);

    // Get the pressure field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    communicator.CommunicateGhosts(1, _pressureMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _pressureMono, _pressureLim);

    // Get the limited kinetic energy by using the limited rho and limited rhou and limited rhov
    // Start in Lagrange formulation then transform to monomial basis in y
    Lkinetic_energy2D(_N_s,_N_E,_rho,_rhou,_rhov,_K); // K = 0.5*(rhou^2+rhov^2)/rho;
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _K, _N_s, 0.0, _KLim, _N_s);

    // Get the limited internal energy by using the limited pressure and limited 1/(gamma-1) in x
    Linternal_energy_multifluid(_N_s1D,_N_E,_N_s1D,_pressureLim,_gammaLim,_rhoeLim);

    // Get the energy field, transform to monomial basis in y
    // Reconstruct the limited energy coefficients in y, back to lagrange
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _E, _N_s, 0.0, _EMono, _N_s);
    Lreconstruct_energy(_N_s1D,_N_E,_N_s1D, _rhoeLim,_KLim,_EMono,_ELim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);

    // Copy the variables back into U (now limited)
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,0     ,_rho  ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,_N_s  ,_rhou ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,2*_N_s,_rhov ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,3*_N_s,_E    ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,4*_N_s,_gamma,U);

    // Get the mass fraction field, transform to monomial basis, limit in y, transform to Lagrange basis, copy to U
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
    communicator.CommunicateGhosts(1, _YMono(x)); \
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _YMono(x), _YLim(x)); \
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s); \
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,(5+x)*_N_s,_Y(x),U);
#include "loop.h"
      
  } // end if cartesian
#else 
  printf("M2L is only implemented for gamncon. Exiting...\n");
  exit(1);
#endif // gamma model

#elif STIFFENED  //===========================================================
  if(_cartesian){

    //
    // Limit wrt x
    // 
      
    // Get the density field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,0,0,U,_rho); // copy from U into rho
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
    communicator.CommunicateGhosts(1, _rhoMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhoMono, _rhoLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);
      
    // Get the momentum x field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,_N_s,0,U,_rhou); // copy from U into rhou
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
    communicator.CommunicateGhosts(1, _rhouMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhouMono, _rhouLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

    // Get the momentum y field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,2*_N_s,0,U,_rhov); // copy from U into rhov
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _rhov, _N_s, 0.0, _rhovMono, _N_s);
    communicator.CommunicateGhosts(1, _rhovMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _rhovMono, _rhovLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _rhovLim, _N_s, 0.0, _rhov, _N_s);
    
    // Get the 1/gamma-1 field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,4*_N_s,0,U,_gamma); // copy from U into gamma
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
    communicator.CommunicateGhosts(1, _gammaMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _gammaMono, _gammaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);

    // Get the gamma pinf/gamma-1 field, transform to monomial basis, limit, transform to Lagrange basis
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,5*_N_s,0,U,_beta); // copy from U into beta
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _beta, _N_s, 0.0, _betaMono, _N_s);
    communicator.CommunicateGhosts(1, _betaMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _betaMono, _betaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _betaLim, _N_s, 0.0, _beta, _N_s);

    // Get the pressure field, transform to monomial basis, limit in x, transform to Lagrange basis
    Lpressure(_N_s, _N_E, U, _pressure);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    communicator.CommunicateGhosts(1, _pressureMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _pressureMono, _pressureLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _pressureLim, _N_s, 0.0, _pressure, _N_s);

    // Get the limited kinetic energy by using the limited rho and limited rhou and limited rhov
    // Start in Lagrange formulation then transform to monomial basis in x
    Lkinetic_energy2D(_N_s,_N_E,_rho,_rhou,_rhov,_K); // K = 0.5*(rhou^2+rhov^2)/rho;
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _K, _N_s, 0.0, _KLim, _N_s); 

    // Get the limited internal energy by using the limited pressure and limited 1/(gamma-1) in x
    Linternal_energy_stiffened(_N_s1D,_N_E,_N_s1D,_pressureLim,_gammaLim,_betaLim,_rhoeLim); 

    // Get the energy field, transform to monomial basis in x
    // Reconstruct the limited energy coefficients in x, back to lagrange
    Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,3*_N_s,0,U,_E); // copy from U into E
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _E, _N_s, 0.0, _EMono, _N_s); 
    Lreconstruct_energy(_N_s1D,_N_E,_N_s1D, _rhoeLim,_KLim,_EMono,_ELim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);

    // Get the mass fraction field, transform to monomial basis, limit in x, transform to Lagrange basis
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) Lstridedcopy(_N_E,_N_s,_N_s*N_F,_N_s,(6+x)*_N_s,0,U,_Y(x)); \
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoX, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
    communicator.CommunicateGhosts(1, _YMono(x)); \
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 0, _YMono(x), _YLim(x)); \
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoX2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s);
#include "loop.h"    

    //
    // Limit wrt y
    //

    // Get the density field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rho, _N_s, 0.0, _rhoMono, _N_s);
    communicator.CommunicateGhosts(1, _rhoMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhoMono, _rhoLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhoLim, _N_s, 0.0, _rho, _N_s);

    // Get the momentum x field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rhou, _N_s, 0.0, _rhouMono, _N_s);
    communicator.CommunicateGhosts(1, _rhouMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhouMono, _rhouLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhouLim, _N_s, 0.0, _rhou, _N_s);

    // Get the momentum y field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _rhov, _N_s, 0.0, _rhovMono, _N_s);
    communicator.CommunicateGhosts(1, _rhovMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _rhovMono, _rhovLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _rhovLim, _N_s, 0.0, _rhov, _N_s);

    // Get the 1/gamma-1 field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _gamma, _N_s, 0.0, _gammaMono, _N_s);
    communicator.CommunicateGhosts(1, _gammaMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _gammaMono, _gammaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _gammaLim, _N_s, 0.0, _gamma, _N_s);

    // Get the gamma*pinf/gamma-1 field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _beta, _N_s, 0.0, _betaMono, _N_s);
    communicator.CommunicateGhosts(1, _betaMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _betaMono, _betaLim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _betaLim, _N_s, 0.0, _beta, _N_s);

    // Get the pressure field, transform to monomial basis, limit in y, transform to Lagrange basis
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    communicator.CommunicateGhosts(1, _pressureMono);
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _pressureMono, _pressureLim);

    // Get the limited kinetic energy by using the limited rho and limited rhou and limited rhov
    // Start in Lagrange formulation then transform to monomial basis in y
    Lkinetic_energy2D(_N_s,_N_E,_rho,_rhou,_rhov,_K); // K = 0.5*(rhou^2+rhov^2)/rho;
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _K, _N_s, 0.0, _KLim, _N_s);

    // Get the limited internal energy by using the limited pressure and limited 1/(gamma-1) in x
    Linternal_energy_stiffened(_N_s1D,_N_E,_N_s1D,_pressureLim,_gammaLim,_betaLim,_rhoeLim);

    // Get the energy field, transform to monomial basis in y
    // Reconstruct the limited energy coefficients in y, back to lagrange
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _E, _N_s, 0.0, _EMono, _N_s);
    Lreconstruct_energy(_N_s1D,_N_E,_N_s1D, _rhoeLim,_KLim,_EMono,_ELim);
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _ELim, _N_s, 0.0, _E, _N_s);

    // Copy the variables back into U (now limited)
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,0     ,_rho  ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,_N_s  ,_rhou ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,2*_N_s,_rhov ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,3*_N_s,_E    ,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,4*_N_s,_gamma,U);
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,5*_N_s,_beta,U);

    // Get the mass fraction field, transform to monomial basis, limit in y, transform to Lagrange basis, copy to U
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2MonoY, _N_s, _Y(x), _N_s, 0.0, _YMono(x), _N_s); \
    communicator.CommunicateGhosts(1, _YMono(x)); \
    Lhrl1D(_N_s1D, _N_E, 1, _N_N, _N_s1D, _neighbors, 2, _YMono(x), _YLim(x)); \
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _MonoY2Lag, _N_s, _YLim(x), _N_s, 0.0, _Y(x), _N_s); \
    Lstridedcopy(_N_E,_N_s,_N_s,_N_s*N_F,0,(6+x)*_N_s,_Y(x),U);
#include "loop.h"

  } // end if cartesian

#endif // problem type
#endif // dimensions

  _timers.stop_timer(20);  
} // end m2limiting

void Limiting::HRIlimiting(COMMUNICATOR &communicator, SENSOR &sensor, scalar* U){
  /*!
    \brief HR limiting function for limiting each element individually
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */

  _timers.start_timer(20);
  
  // Calculate the sensor
  sensor.sensing(_neighbors,U);

#ifdef ONED

  // Call limiting function
  Lhri1D(_N_s, _N_E, _N_N, _neighbors, _N_s, 1, 0, _Lag2Mono, _Mono2Lag, sensor.getSensors(), U, _Utmp);
  sensor.copy_detected_elements(_Utmp, U);
  
#elif TWOD
  if(_cartesian){

    // Call limiting function in x
    Lhri1D(_N_s, _N_E, _N_N, _neighbors, _N_s1D, _N_s1D, 0, _Lag2MonoX, _MonoX2Lag, sensor.getSensors(), U, _Utmp);
    sensor.copy_detected_elements(_Utmp, U);
    
    // Communicate the elements on different partitions if necessary
    communicator.CommunicateGhosts(N_F, U);

    // Call limiting function in y
    Lhri1D(_N_s, _N_E, _N_N, _neighbors, _N_s1D, _N_s1D, 2, _Lag2MonoY, _MonoY2Lag, sensor.getSensors(), U, _Utmp);
    sensor.copy_detected_elements(_Utmp, U);
  }
#endif
  
  _timers.stop_timer(20);
}

void Limiting::PRIlimiting(COMMUNICATOR &communicator, SENSOR &sensor, scalar* U){
  /*!
    \brief PR limiting function for limiting each element individually (primitive variables)
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */

  // Transform to primitive variables
  _timers.start_timer(20);
  LCons2Prim(_N_s, _N_E, U);
  _timers.stop_timer(20);
  
  // Limit primitive variables  
  HRIlimiting(communicator, sensor, U);

  // Transform to conserved variables
  _timers.start_timer(20);
  LPrim2Cons(_N_s, _N_E, U);
  _timers.stop_timer(20);
}

void Limiting::M2Ilimiting(COMMUNICATOR &communicator, SENSOR &sensor, scalar* U){
  /*!
    \brief modified limiting function (see jcp2014) for limiting each element individually
    \param[in] communicator communicator object for MPI communications
    \param[out] U solution to be limited
  */
  
  _timers.start_timer(20);
  
  // Calculate the sensor
  sensor.sensing(_neighbors,U);

#ifdef ONED

  // Call limiting function
  Lm2i1D(_N_s, _N_E, _N_N, _neighbors, _N_s, 1, 0, _Lag2Mono, _Mono2Lag, sensor.getSensors(), U, _Utmp);
  sensor.copy_detected_elements(_Utmp, U);
  
#elif TWOD
  if(_cartesian){

    // Call limiting function in x
    Lm2i1D(_N_s, _N_E, _N_N, _neighbors, _N_s1D, _N_s1D, 0, _Lag2MonoX, _MonoX2Lag, sensor.getSensors(), U, _Utmp);
    sensor.copy_detected_elements(_Utmp, U);
    
    // Communicate the elements on different partitions if necessary
    communicator.CommunicateGhosts(N_F, U);

    // Call limiting function in y
    Lm2i1D(_N_s, _N_E, _N_N, _neighbors, _N_s1D, _N_s1D, 2, _Lag2MonoY, _MonoY2Lag, sensor.getSensors(), U, _Utmp);
    sensor.copy_detected_elements(_Utmp, U);
  }
#endif

  _timers.stop_timer(20);
}
