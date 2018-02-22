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
  //PEJ 11/09/2017: ALtering rk routines to allow for non-diagonal butcher table
 private:
  int     _order;  // RK order (only implemented RK4)
  int     _steps; //steps in the RK integration scheme; 1 for RK1, 4 for RK1, 13 for RK8
  scalar* _beta;  //substep time (Tau) at which each RK sub-residual is calculated
  scalar* _gamma; //how much each substep influences master solution
  scalar* _alpha; //Updates predictor solution from known substep residuals
  std::vector<double> _output_time_array;
  double _Tf;
  //scalar _Tf;
  
 public:
  /*!
    \brief Constructor
    \param[in] order RK scheme order
    \param[in] steps RK scheme step count
    \param[in] DtOut output time step
    \param[in] Tf final time
  */
 RK(int order, int steps, double DtOut, double Tf, const std::vector<double> &output_time_array = std::vector<double>()) : _order(order), _steps(steps), _Tf(Tf){
    switch (_order){
    case 1:
      //I can get away without using alpha, but it's he
      _beta = new scalar[1];
      _beta[0] = 0.0;
      _gamma = new scalar[1];
      _gamma[0] = 1.0;
      _alpha = new scalar[1];
      _alpha[0] = 0.0;
      break;

    case 2:
      // This is a case for optimal RK2-TVD method as suggested by
      // Gottlieb and Shu (1998), Mathematics of Computation 
      _beta = new_scalar[2];
      _gamma = new_scalar[2];
      -beta[0] = 0.0; -beta[1] = 1.0;
      _gamma[0] = 1.0/2.0; _gamma[1] = 1.0/2.0;
      break;

    case 3:
      // This is a case for optimal RK3-TVD method as suggested by
      // Gottlieb and Shu (1998), Mathematics of Computation
      _beta = new_scalar[3];
      _gamma = new_scalar[3];
      -beta[0] = 0.0; -beta[1] = 1.0; _beta[2] = 0.5;
      _gamma[0] = 1.0/6.0; _gamma[1] = 1.0/6.0; _gamma[3] = 2.0/3.0;
      break;
   
    case 4:
      //I can get away without using alpha, but it's here to help me debug the RK[8/13] case
      _beta = new scalar[4];
      _beta[0] = 0.0; _beta[1] = 0.5; _beta[2] = 0.5; _beta[3] = 1.0;
      _gamma = new scalar[4];
      _gamma[0] = 1.0/6.0; _gamma[1] = 2.0/6.0; _gamma[2] = 2.0/6.0; _gamma[3] = 1.0/6.0;
      _alpha = new scalar[4*4];
      for (int i1 = 0; i1 < 4; i1 = i1 + 1){
	    for (int i2 = 0; i2 < 4; i2 = i2 + 1){
	      _alpha[i1*4 + i2] = 0.0;
	    }}
      _alpha[0*4 + 0] = 0.5;
      _alpha[1*4 + 0] = 0.0;
      _alpha[1*4 + 1] = 0.5;
      _alpha[2*4 + 0] = 0.0;
      _alpha[2*4 + 1] = 0.0;
      _alpha[2*4 + 2] = 1.0;
      _alpha[3*4 + 0] = 0.0;
      _alpha[3*4 + 1] = 0.0;
      _alpha[3*4 + 2] = 0.0;
      _alpha[3*4 + 3] = 1.0;
      break;
    case 8:
      {
	//this is different from rk1 and rk4 in the fact that we need to account for the alpha matrix
	_beta = new scalar[13];
	_gamma = new scalar[13];
	_alpha = new scalar[13*13];
	//not sure how to write a subroutine to populate these coefficients, so it's gonna happen here.
	_beta[0] = 0.0;
	_beta[1] = 1.0 / 18.0;
	_beta[2] = 1.0 / 12.0;
	_beta[3] = 1.0 / 8.0;
	_beta[4] = 5.0 / 16.0;
	_beta[5] = 3.0 / 8.0;
	_beta[6] = 59.0 / 400.0;
	_beta[7] = 93.0 / 200.0;
	_beta[8] = 5490023248.0 / 9719169821.0;
	_beta[9] = 13.0 / 20.0;
	_beta[10] = 1201146811.0 / 1299019798;
	_beta[11] = 1.0;
	_beta[12] = 1.0;

	_gamma[0] = 14005451.0 / 335480064;
	_gamma[1] = 0.0;
	_gamma[2] = 0.0;
	_gamma[3] = 0.0;
	_gamma[4] = 0.0;
	_gamma[5] = -59238493.0 / 1068277825.0;
	_gamma[6] = 181606767.0 / 758867731.0;
	_gamma[7] = 561292985.0 / 797845732.0;
	_gamma[8] = -1041891430.0 / 1371343529.0;
	_gamma[9] = 760417239.0 / 1151165299.0;
	_gamma[10] = 118820643.0 / 751138087.0;
	_gamma[11] = -528747749.0 / 2220607170.0;
	_gamma[12] = 1.0 / 4.0;
	 
	for (int i1 = 0; i1 < 13; i1 = i1 + 1){
	    for (int i2 = 0; i2 < 13; i2 = i2 + 1){
		_alpha[i1*13 + i2] = 0.0;
	    }}
	_alpha[0*13 + 0] = 1.0 / 18.0;
	_alpha[1*13 + 0] = 1.0 / 48.0;
	_alpha[2*13 + 0] = 1.0 / 32.0;
	_alpha[3*13 + 0] = 5.0 / 16.0;
	_alpha[4*13 + 0] = 3.0 / 80.0;
	_alpha[5*13 + 0] = 29443841.0 / 614563906.0;
	_alpha[6*13 + 0] = 16016141.0 / 946692911.0;
	_alpha[7*13 + 0] = 39632708.0 / 573591083.0;
	_alpha[8*13 + 0] = 246121993.0 / 1340847787.0;
	_alpha[9*13 + 0] = -1028468189.0 / 846180014.0;
	_alpha[10*13 + 0] = 185892177.0 / 718116043.0;
	_alpha[11*13 + 0] = 403863854.0 / 491063109.0;
	
	_alpha[1*13 + 1] = 1.0 / 16.0;
	//lots of zeros in the 1 column
	
	_alpha[2*13 + 2] = 3.0 / 32.0;
	_alpha[3*13 + 2] = -75.0 / 64.0;
	_alpha[4*13 + 2] = 0.0;
	//lots of zeros in the 2 column
	
	_alpha[3*13 + 3] = 75.0 / 64.0;
	_alpha[4*13 + 3] = 3.0 / 16.0;
	_alpha[5*13 + 3] = 77736538.0 / 692538347.0;
	_alpha[6*13 + 3] = 61564180.0 / 158732637.0;
	_alpha[7*13 + 3] = -433636366.0 / 683701615.0;
	_alpha[8*13 + 3] = -37695042795.0 / 15268766246.0;
	_alpha[9*13 + 3] = 8478235783.0 / 508512852.0;
	_alpha[10*13 + 3] = -3185094517.0 / 667107341.0;
	_alpha[11*13 + 3] = -5068492393.0 / 434740067.0;
	
	_alpha[4*13 + 4] = 3.0 / 20.0;
	_alpha[5*13 + 4] = -28693883.0 / 1125000000.0;
	_alpha[6*13 + 4] = 22789713.0 / 633445777.0;
	_alpha[7*13 + 4] = -421739975.0 / 2616292301.0;
	_alpha[8*13 + 4] = -309121744.0 / 1061227803.0;
	_alpha[9*13 + 4] = 1311729495.0 / 1432422823.0;
	_alpha[10*13 + 4] = -477755414.0 / 1098053517.0;
	_alpha[11*13 + 4] = -411421997.0 / 543043805.0;
	
	_alpha[5*13 + 5] = 23124283.0 / 1800000000.0;
	_alpha[6*13 + 5] = 545815736.0 / 2771057229.0;
	_alpha[7*13 + 5] = 100302831.0 / 723423059.0;
	_alpha[8*13 + 5] = -12992083.0 / 490766935.0;
	_alpha[9*13 + 5] = -10304129995.0 / 1701304382.0;
	_alpha[10*13 + 5] = -703635378.0 / 230739211.0;
	_alpha[11*13 + 5] = 652783627.0 / 914296604.0;
	
	_alpha[6*13 + 6] = -180193667.0 / 1043307555.0;
	_alpha[7*13 + 6] = 790204164.0 / 839813087.0;
	_alpha[8*13 + 6] = 6005943493.0 / 2108947869.0;
	_alpha[9*13 + 6] = -48777925059.0 / 3047939560.0;
	_alpha[10*13 + 6] = 5731566787.0 / 1027545527.0;
	_alpha[11*13 + 6] = 11173962825.0 / 925320556.0;
	
	_alpha[7*13 + 7] = 800635310.0 / 3783071287.0;
	_alpha[8*13 + 7] = 393006217.0 / 1396673457.0;
	_alpha[9*13 + 7] = 15336726248.0 / 1032824649.0;
	_alpha[10*13 + 7] = 5232866602.0 / 850066563.0;
	_alpha[11*13 + 7] = -13158990841.0 / 6184727034.0;
	
	_alpha[8*13 + 8] = 123872331.0 / 1001029789.0;
	_alpha[9*13 + 8] = -45442868181.0 / 3398467696.0;
	_alpha[10*13 + 8] = -4093664535.0 / 808688257.0;
	_alpha[11*13 + 8] = 3936647629.0 / 1978049680.0;
	
	_alpha[9*13 + 9] = 3065993473.0 / 597172653.0;
	_alpha[10*13 + 9] = 3962137247.0 / 1805957418.0;
	_alpha[11*13 + 9] = -160528059 / 685178525.0;
	
	_alpha[10*13 + 10] = 65686358.0 / 487910083.0;
	_alpha[11*13 + 10] = 248638103.0 / 1413531060.0;
	
	_alpha[11*13 + 11] = 0.0;
	
	break;
      }
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
      printf("In RK constructor: currently, output time array size is 0 and DTout>0\n");
      // Number of ouputs

      int num_outputs = int(Tf/DtOut)+1; 
      printf("Tf=%f, DtOut=%f, num_outputs=%d\n",Tf,DtOut,num_outputs);
      // Populate the array
      for(int k=0; k < num_outputs; k ++){
	_output_time_array.push_back(k*DtOut);
      }

      // Enforce _Tf to be identically exactly the same as the final
      // time calculated from the for loop that populated the output
      // time array
      _Tf = _output_time_array.back(); 
      printf("Using output time array, set _Tf=%f\n", _Tf);
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
    if(_alpha)       delete[] _alpha;
  };

  void RK_integration(scalar CFL, scalar VNN, scalar VNNAD, int restart_step,
		      int N_E, int N_s, int N_G, int M_T, int M_s, int N_ghosts, int N_N,
		      scalar* h_Minv, 
		      scalar* h_U, int* neighbors,
		      Limiting &Limiter, bool order0, DG_SOLVER &dgsolver, COMMUNICATOR &communicator, PRINTER &printer, SENSOR &sensor, TIMERS &timers, MEM_COUNTER &mem_counter, LAGRANGE_PARTICLES &particles);
    
  scalar DtFromCFL(const int N_s, const int N_E, const scalar CFL, scalar* U, scalar* UPA);
  scalar DtFromCFL_and_VNN(const int N_s, const int N_E, const scalar CFL, const scalar VNN, scalar* U, scalar* UPA, scalar* Mew);
  scalar DtFromCFL_and_VNN_and_VNNAD(const int N_s, const int N_E, const scalar CFL, const scalar VNN, const scalar VNNAD, scalar* U, scalar* UPA, scalar* Mew, scalar* DivMax, scalar* ADepsMax, scalar Beta_S, scalar Mew_S);
  
}; //end the RK clcass

#endif // RK_H

