/*!
  \file rk.h
  \brief Runge-Kutta time integration class
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup rk Runge-Kutta
*/
#ifndef RK_H
#define RK_H

#include <macros.h>
#include <rk_kernels.h>
#include <kernels.h>
#include <limiting.h>
#include <dg_solver.h>
#include <communicator.h>
#include <printer.h>
#include <sensor.h>
#ifdef USE_MPI
#include "mpi.h"
#endif

class RK
{
 private:
  int     _order;  // RK order (only implemented RK4)
  scalar* _beta;
  scalar* _gamma;
  
 public:
  /*!
    \brief Constructor
    \param[in] order DG polynomial order
  */
  RK(int order) : _order(order){
    switch (_order){
    case 4:
      _beta = new scalar[4];
      _beta[0] = 0.0; _beta[1] = 0.5; _beta[2] = 0.5; _beta[3] = 1.0;
      _gamma = new scalar[4];
      _gamma[0] = 1.0/6.0; _gamma[1] = 2.0/6.0; _gamma[2] = 2.0/6.0; _gamma[3] = 1.0/6.0;
      break;
    default:
      printf("Invalid RK order. Defaulting to RK4.\n");
      _beta = new scalar[4];
      _beta[0] = 0.0; _beta[1] = 0.5; _beta[2] = 0.5; _beta[3] = 1.0;
      _gamma = new scalar[4];
      _gamma[0] = 1.0/6.0; _gamma[1] = 2.0/6.0; _gamma[2] = 2.0/6.0; _gamma[3] = 1.0/6.0;
    }
  };

  /*! Destructor*/
  ~RK(){
    if(_beta)        delete[] _beta;
    if(_gamma)       delete[] _gamma;
  };

  void RK_integration(double DtOut, double Tf, scalar CFL,
		      int N_E, int N_s, int N_G, int M_T, int M_s, int N_ghosts,
		      scalar* h_Minv, 
		      scalar* h_U,
		      Limiting &Limiter, bool order0, DG_SOLVER &dgsolver, COMMUNICATOR &communicator, PRINTER &printer, SENSOR &sensor);
    
  scalar DtFromCFL(const int N_s, const int N_E, const scalar CFL, scalar* U, scalar* UPA);
  
}; //end the RK clcass

#endif // RK_H

