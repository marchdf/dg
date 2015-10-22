/*!
  \file rk.h
  \brief Runge-Kutta time integration class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef RK_H
#define RK_H

#include <vector>
#include "macros.h"
#include "rk_kernels.h"
#include "kernels.h"
#include "limiting.h"
#include "dg_solver.h"
#include "communicator.h"
#include "printer.h"
#include "sensor.h"
#include "timers.h"
#include "mem_counter.h"
#include "lagrange_particles.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

class RK
{
 private:
  int     _order;  // RK order (only implemented RK4)
  scalar* _beta;
  scalar* _gamma;
  std::vector<double> _output_time_array;
  double _Tf;
  
 public:
  /*!
    \brief Constructor
    \param[in] order DG polynomial order
    \param[in] DtOut output time step
    \param[in] Tf final time
  */
  RK(int order, double DtOut, double Tf, const std::vector<double> &output_time_array = std::vector<double>()) : _order(order), _Tf(Tf){
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

    // Calculate the output time array
    // if the array was not specified in the deck, default to constant DtOut
    if ((output_time_array.size() == 0) && (DtOut > 0)){
      scalar next_time = 0;
      while (next_time <= Tf){
	_output_time_array.push_back(next_time);
	next_time = next_time+DtOut;
      }
    }
    // if it was specified in the deck, use it
    else if (output_time_array.size() != 0){
      _output_time_array = output_time_array;
      _Tf = output_time_array.back();
    }
    else{
      printf("Could not figure out the output times.\nCheck the deck: either specify the output array OR the output delta t\n");
      exit(1);
    }

      
  };

  /*! Destructor*/
  ~RK(){
    if(_beta)        delete[] _beta;
    if(_gamma)       delete[] _gamma;
  };

  void RK_integration(scalar CFL, int restart_step,
		      int N_E, int N_s, int N_G, int M_T, int M_s, int N_ghosts,
		      scalar* h_Minv, 
		      scalar* h_U,
		      Limiting &Limiter, bool order0, DG_SOLVER &dgsolver, COMMUNICATOR &communicator, PRINTER &printer, SENSOR &sensor, TIMERS &timers, MEM_COUNTER &mem_counter, LAGRANGE_PARTICLES &particles);
    
  scalar DtFromCFL(const int N_s, const int N_E, const scalar CFL, scalar* U, scalar* UPA);
  
}; //end the RK clcass

#endif // RK_H

