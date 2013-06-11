//
// Runge-Kutta class
//
#ifndef RK_H
#define RK_H

#include <macros.h>
#include <rk_kernels.h>
#include <cpu_kernels.h>
#include <limiting.h>
#include <dg_solver.h>
#include <print_sol.h>
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
  };

  // destructor
  ~RK(){
    if(_beta)        delete[] _beta;
    if(_gamma)       delete[] _gamma;
  };
  
  // main RK integration function
  void RK_integration(double DtOut, double Tf, scalar CFL,
		      int N_E, int N_s, int N_G, int M_T, int M_s,
		      scalar* h_Minv, 
		      scalar* h_U,
		      Limiting &Limiter, bool order0, DG_SOLVER &dgsolver,
		      int elem_type, simpleMesh &m){

    // Initialize some vars
    double T = 0;          // current time
    double Tstar = 0;
    double Tout = T+DtOut; // next output time
    scalar Dt = 0;
    scalar DtCFL = 0;
    int n = 0;             // counts the time steps
    int count = 1;         // counts the output steps
    bool output = false;
    bool done = false;
    
    scalar* _Us;
    scalar* _Ustar;
    scalar* _DU;
    scalar* _UPA; 
    scalar* _f;
    scalar* _Minv;    
#ifdef USE_CPU
    _Us    = new scalar[N_s*N_E*N_F];  makeZero(_Us   ,N_s*N_E*N_F);	 
    _Ustar = new scalar[N_s*N_E*N_F];  makeZero(_Ustar,N_s*N_E*N_F);
    _DU    = new scalar[N_s*N_E*N_F];  makeZero(_DU   ,N_s*N_E*N_F);
    _UPA   = new scalar[N_s*N_E];      makeZero(_UPA  ,N_s*N_E);
    _f     = new scalar[N_s*N_E*N_F];  makeZero(_f    ,N_s*N_E*N_F);
    _Minv  = new scalar[N_s*N_s*N_E];  memcpy(_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar));
#elif USE_GPU
    scalar* d_U;
    // Allocate on device
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_U   , N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_Us   , N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_Ustar, N_s*N_E*N_F*sizeof(scalar)));  
    CUDA_SAFE_CALL(cudaMalloc((void**) &_DU   , N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_UPA  , N_s*N_E*sizeof(scalar)));     
    CUDA_SAFE_CALL(cudaMalloc((void**) &_f    , N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_Minv , N_s*N_s*N_E*sizeof(scalar)));

    // Set to zero
    CUDA_SAFE_CALL(cudaMemset(_Us   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMemset(_Ustar, (scalar)0.0, N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMemset(_DU   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMemset(_UPA  , (scalar)0.0, N_s*N_E*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMemset(_f    , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar)));
    
    // Copy info to device
    CUDA_SAFE_CALL(cudaMemcpy(d_U, h_U, N_s*N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar), cudaMemcpyHostToDevice));
#endif
  
    /* // Limit solution */
    /* if      (Limiter.getLimitingMethod()==1) Limiter.HRlimiting(arch(U)); */
    /* else if (Limiter.getLimitingMethod()==2) Limiter.MYlimiting(arch(U)); */
    /* else if (Limiter.getLimitingMethod()==3) Limiter.M2limiting(arch(U)); */

    // print the initial condition to the file
    printf("Initial condition written to output file.\n");
    print_dg(N_s, N_E, h_U, m, elem_type, 0, 0, 0);
    
    // Output conservation of the fields
    dgsolver.conservation(h_U,0.0);
	   
    // Time integration
    while (!done){
      // Find new Dt
      Dt = DtFromCFL(N_s, N_E, CFL, arch(U),_UPA); output = false;
      printf("helloS Dt = %e\n",Dt);
#ifdef USE_MPI // Make sure everyone is at the same time
      MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
      MPI_Allreduce(MPI_IN_PLACE, &Dt, 1, MPI_SCALAR, MPI_MIN, MPI_COMM_WORLD);
#endif
      if(Dt<1e-14){ printf("Next time step is too small (%e<1e-14). Exiting.\n",Dt); exit(1);}
      if     (Dt>(Tf  -T)){ DtCFL = Dt; Dt = Tf  -T; output = true; done = true;}
      else if(Dt>(Tout-T)){ DtCFL = Dt; Dt = Tout-T; output = true;}
      /* printf("current time=%f, this Dt=%f, next output at %f\n",T+Dt,Dt,Tout); */
      /* Dt = 1e-7; */
      /* if ((n+1)%100==0){output=true;} */
      /* else {output=false;} */
      /* if ((n+1)==1000) {done = true;} */

      // Us = U
#ifdef HAVE_BLAS
      blasCopy(N_F*N_s*N_E, arch(U), 1, _Us, 1);
#else
      Lcpu_equal(N_s, N_E, N_F, _Us, arch(U));
#endif
      for(int k = 0; k < _order; k++){
	// Ustar = Us + beta*DU
#ifdef HAVE_BLAS
	blasCopy(N_F*N_s*N_E, _Us, 1, _Ustar, 1);
	blasAxpy(N_s*N_F*N_E, _beta[k], _DU, 1, _Ustar, 1);
#else
	Lcpu_equal(N_s, N_E, N_F, _Ustar, _Us); // make Ustar = Us;
	Lcpu_add(N_s, N_E, N_F, _Ustar, _DU, _beta[k]);// do Ustar.add(DU,beta[k]);
#endif

	Tstar = T + _beta[k]*Dt;

	//Limit the solution if you so want to do so
	if(k>0){
	  if      (Limiter.getLimitingMethod()==1) Limiter.HRlimiting(_Ustar);
	  else if (Limiter.getLimitingMethod()==2) Limiter.MYlimiting(_Ustar);
	  else if (Limiter.getLimitingMethod()==3) Limiter.M2limiting(_Ustar);
	}

	// Now you have to calculate f(Ustar)
	dgsolver.dg_solver(_Ustar,_f);
	
	// Solve: DU = Dt*Minv*f(Ustar)
	Lsolver(N_s, N_E, Dt, _Minv, _f, _DU);

	// if 0-order average the solution in the cells
	if (order0){Laverage_cell_p0(N_s, N_E, _DU);}

	// U = U + gamma*DU
#ifdef HAVE_BLAS
	blasAxpy(N_s*N_F*N_E, _gamma[k], _DU, 1, arch(U), 1);
#else
	Lcpu_add(N_s, N_E, N_F, arch(U), _DU, _gamma[k]); // do U.add(DU,gamma[k])
#endif

      }// end loop on k
      
      // Limit solution
      if      (Limiter.getLimitingMethod()==1) Limiter.HRlimiting(arch(U));
      else if (Limiter.getLimitingMethod()==2) Limiter.MYlimiting(arch(U));
      else if (Limiter.getLimitingMethod()==3) Limiter.M2limiting(arch(U));


      T = T + Dt; // update current time
      n++;        // update the time step counter

      // Output the solution
      if(output){

	// Get the solution to the CPU
#ifdef  USE_GPU
	CUDA_SAFE_CALL(cudaMemcpy(h_U, d_U, N_s*N_F*N_E*sizeof(scalar), cudaMemcpyDeviceToHost));
#endif

	printf("Solution written to output file at step %7i and time %e (current CFL time step:%e).\n",n,T,DtCFL);
	print_dg(N_s, N_E, h_U, m, elem_type, count, T, 1);
	Tout = T + DtOut; // update the new output time
	count++;

	// Output conservation of the fields
	dgsolver.conservation(h_U,T);
	
      }// end output steps
      printf("helloE\n");
    }// end loop on time
    
    // Free some stuff
    del(_Us);
    del(_Ustar);
    del(_DU);
    del(_UPA);
    del(_f);
    del(_Minv);
#ifdef USE_GPU
    CUDA_SAFE_CALL(cudaFree(d_U));
#endif
  }; //end main RK function

  scalar DtFromCFL(const int N_s, const int N_E, const scalar CFL, scalar* U, scalar* UPA){
    LfindUPA(N_s, N_E, U, UPA);
    int maxUPAIdx;
    maxUPAIdx = blasIamax(N_E*N_s,UPA,1)-1; // Fortran starts numbering at 1
#ifdef USE_CPU
    scalar maxUPA = UPA[maxUPAIdx];
    scalar Dt = CFL/maxUPA;
#elif USE_GPU
    scalar* maxUPA = new scalar[1];
    CUDA_SAFE_CALL(cudaMemcpy(maxUPA, &UPA[maxUPAIdx], sizeof(scalar), cudaMemcpyDeviceToHost));
    scalar Dt = CFL/maxUPA[0];
    delete[] maxUPA;
#endif
    return Dt;
  }

  
}; //end the RK clcass

#endif // RK_H

