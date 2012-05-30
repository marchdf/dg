//
// Runge-Kutta class
//
#ifndef RK_H
#define RK_H

#include <stdio.h>
#include <scalar_def.h>
#include <macros.h>
//#include <common_kernels.h>
#include <cpu_kernels.h>
#include <limiting.h>

class RK
{
 private:
  int     _order;  // RK order (only implemented RK4)
  scalar* _beta;
  scalar* _gamma;

 public:
  // constructor
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
  }

  // destructor
  ~RK(){
    if(_beta)        delete[] _beta;
    if(_gamma)       delete[] _gamma;
  }
  
  // main RK integration function
  void RK_integration(scalar Dt, int Nt, int output_factor, scalar* h_U, int blas, 
		      int N_s, int N_E, int M_T, int N_G, int N_F, Limiting &Limiter){

  
    // Initialize some vars
    double T = 0;
    double Tstar = 0;
    int count = 1;
#ifdef USE_CPU
    scalar* h_Us    = new scalar[N_s*N_E*N_F];
    scalar* h_Ustar = new scalar[N_s*N_E*N_F];
    scalar* h_DU    = new scalar[N_s*N_E*N_F];
#elif USE_GPU
    scalar* d_U, *d_Us, *d_Ustar, *d_DU, *d_UF;
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_U,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMemcpy(d_U, h_U, N_s*N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Us,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Ustar,N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_DU,N_s*N_E*N_F*sizeof(scalar)));
#endif

    // ATTENTION: need to initialize the first step with euler.
    // Need to do a solve operation (and average the solution for p=0)
    
    // Time integration
    for (int n = 1; n <= Nt; n++){

      // Us = U
      if (blas==1) {blasCopy(N_F*N_s*N_E, arch(U), 1, arch(Us), 1);}
      else Lcpu_equal(N_s, N_E, N_F, arch(Us), arch(U));

      for(int k = 0; k < _order; k++){
	// Ustar = Us + beta*DU
	if (blas==1) {blasCopy(N_F*N_s*N_E, arch(Us), 1, arch(Ustar), 1);}
	else Lcpu_equal(N_s, N_E, N_F, arch(Ustar), arch(Us)); // make Ustar = Us;
	if (blas==1) {blasAxpy(N_s*N_F*N_E, _beta[k], arch(DU), 1, arch(Ustar), 1);}
	else Lcpu_add(N_s, N_E, N_F, arch(Ustar), arch(DU), _beta[k]);// do Ustar.add(DU,beta[k]);
	Tstar = T + _beta[k]*Dt;

	//Limit the solution if you so want to do so
	if(k>0){if (Limiter.getLimitingStatus()) Limiter.HRlimiting(arch(Ustar));}

	// Solve: DU = Dt*Minv*f(Ustar)
	// solve: requires Q, F, S, Dt, Minv, DU 
	// Lcpu_solve(N_s, N_E, N_F, arch(DU), arch(S), arch(F), arch(Q), arch(Minv), Dt);

	// if 0-order average the solution in the cells
	//if (order0){Lcpu_average_cell_p0(N_s, N_E, N_F, arch(DU));}

	// U = U + gamma*DU
	/* if (blas==1) {blasAxpy(N_s*N_F*N_E, _gamma[k], arch(DU), 1, arch(U), 1);}       */
	/* else Lcpu_add(N_s, N_E, N_F, arch(U), arch(DU), _gamma[k]); // do U.add(DU,gamma[k]) */

      }// end loop on k
      
      T = T + Dt;

      // Limit solution
      //if (Limiter.getLimitingStatus()) Limiter.HRlimiting(arch(U));

      // Output the solution
      if(n % (Nt/output_factor) == 0){

	// Get the solution to the CPU
#ifdef  USE_GPU
	CUDA_SAFE_CALL(cudaMemcpy(h_U, d_U, N_s*N_F*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
#endif
		     
	//printf("Solution written to output file at step %i and time %f.\n",n,n*Dt);
	/* if(multifluid)print_dg_multifluid(N_s, N_E, N_F, model, h_U, m, msh_lin, count, n*Dt, 1,-1); */
	/* if(passive)   print_dg_passive(N_s, N_E, N_F, gamma0, h_U, m, msh_lin, count, n*Dt, 1,-1); */
	count++;
      }


    }// end loop on time


    // Free some stuff
#ifdef USE_CPU
    delete[] h_Us;
    delete[] h_Ustar;
    delete[] h_DU;
#elif USE_GPU
    CUDA_SAFE_CALL(cudaFree(d_U));
    CUDA_SAFE_CALL(cudaFree(d_Us));
    CUDA_SAFE_CALL(cudaFree(d_Ustar));
    CUDA_SAFE_CALL(cudaFree(d_DU));
#endif
  };

  
};

#endif // RK_H
