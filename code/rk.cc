/*!
  \file rk.cc
  \brief Runge-Kutta time integration function definitions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup rk
*/
#include "rk.h"

void RK::RK_integration(double DtOut, double Tf, scalar CFL, int restart_step,
			int N_E, int N_s, int N_G, int M_T, int M_s, int N_ghosts,
			scalar* h_Minv, 
			scalar* h_U,
			Limiting &Limiter, bool order0, DG_SOLVER &dgsolver, COMMUNICATOR &communicator, PRINTER &printer, SENSOR &sensor, TIMERS &timers, MEM_COUNTER &mem_counter, LAGRANGE_PARTICLES &particles){
  /*!
    \brief Main RK integration function
    \param[in] DtOut output time step
    \param[in] Tf final time
    \param[in] CFL CFL number
    \param[in] restart_step output step for a restart
    \param[in] N_E number of elements
    \param[in] N_s number of nodes per element
    \param[in] N_G number of gaussian nodes per element
    \param[in] M_T number of interfaces
    \param[in] M_s number of nodes per interface
    \param[in] N_ghosts number of ghost elements
    \param[in] h_Minv host array containing the inverse mass matrices of each element
    \param[out] h_U solution to integrate in time
    \param[in] Limiter limiter object
    \param[in] order0 true if DG p=0
    \param[in] dgsolver solver object
    \param[in] communicator communicator object
    \param[in] printer printer object
    \param[in] sensor sensor object
    \param[in] timers timers object
    \param[in] mem_counter memory counter object
    \param[in] particles lagrange particles object
  */

  // Arrays
  scalar* _Us;
  scalar* _Ustar;
  scalar* _DU;
  scalar* _UPA; 
  scalar* _f;
  scalar* _Minv;
  
#ifdef USE_CPU
  _Us    = new scalar[N_s*N_E*N_F];  makeZero(_Us   ,N_s*N_E*N_F);	                  mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _Ustar = new scalar[N_s*(N_E+N_ghosts)*N_F];  makeZero(_Ustar,N_s*(N_E+N_ghosts)*N_F);  mem_counter.addToCPUCounter(N_s*(N_E+N_ghosts)*N_F*sizeof(scalar));
  _DU    = new scalar[N_s*N_E*N_F];  makeZero(_DU   ,N_s*N_E*N_F);                        mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _UPA   = new scalar[N_s*N_E];      makeZero(_UPA  ,N_s*N_E);                            mem_counter.addToCPUCounter(N_s*N_E*sizeof(scalar));
  _f     = new scalar[N_s*N_E*N_F];  makeZero(_f    ,N_s*N_E*N_F);                        mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _Minv  = new scalar[N_s*N_s*N_E];  memcpy(_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar));   mem_counter.addToCPUCounter(N_s*N_s*N_F*sizeof(scalar));
#elif USE_GPU
  scalar* d_U;
  // Allocate on device
  cudaMalloc((void**) &d_U   , N_s*N_E*N_F*sizeof(scalar));                               mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));
  cudaMalloc((void**) &_Us   , N_s*N_E*N_F*sizeof(scalar));                               mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_Ustar, N_s*N_E*N_F*sizeof(scalar));  				  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar)); 
  cudaMalloc((void**) &_DU   , N_s*N_E*N_F*sizeof(scalar));				  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_UPA  , N_s*N_E*sizeof(scalar));     				  mem_counter.addToGPUCounter(N_s*N_E*sizeof(scalar));		      
  cudaMalloc((void**) &_f    , N_s*N_E*N_F*sizeof(scalar));				  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_Minv , N_s*N_s*N_E*sizeof(scalar));				  mem_counter.addToGPUCounter(N_s*N_s*N_F*sizeof(scalar));
  
  // Set to zero
  cudaMemset(_Us   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_Ustar, (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_DU   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_UPA  , (scalar)0.0, N_s*N_E*sizeof(scalar));
  cudaMemset(_f    , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  
  // Copy info to device
  cudaMemcpy(d_U, h_U, N_s*N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice);
  cudaMemcpy(_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar), cudaMemcpyHostToDevice);
#endif

  // Get cpu id
  int myid = 0;
#ifdef USE_MPI 
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

  // Initialize some vars
  double T = 0;             // current time
  int count = restart_step; // counts the output steps
  double Tstar = 0;
  double Tout = 0;          // next output time
  scalar Dt = 0;
  scalar DtCFL = 0;
  int n = 0;                // counts the time steps
  bool output = false;
  bool done = false;

  // If count is not equal to zero, read output files for restart
  if (count!=0){
    printer.read(count,T,arch(U));
  }
  else{ // print the initial condition to the file
    if(CFL>=0){
      timers.start_timer(2);
      
      if(myid==0){printf("Initial condition written to output file.\n");}
      printer.print(arch(U), particles, count, T);
      
      // Limit the initial solution before integrating to avoid problems
      if      (Limiter.getLimitingMethod()==1){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRlimiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==2){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2limiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==3){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRIlimiting(communicator, sensor, arch(U));}
      else if (Limiter.getLimitingMethod()==4){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2Ilimiting(communicator, sensor, arch(U));}
      else if (Limiter.getLimitingMethod()==5){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRlimiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==6){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRIlimiting(communicator, sensor, arch(U));}
      printer.print_sensor(sensor, count, T);
    
      // Output conservation of the fields (if wanted). Works only
      // when data is on host (really this is just for small
      // validation runs)
#ifdef CONS
      dgsolver.conservation(h_U,T);
#endif
      
      timers.stop_timer(2);
    }
  }
  Tout = T+DtOut;
  count++;

  // Time integration
  timers.start_timer(1);
  while (!done){

    // Give the deck a negative CFL for fixed time step and no output
    timers.start_timer(4);
    if(CFL<0){
      Dt = DtOut; output = false;
      if(Dt>(Tf  -T)){Dt = Tf  -T; done = true;}
    }
    else{
      // Find new Dt
      Dt = DtFromCFL(N_s, N_E, CFL, arch(U),_UPA); output = false;
#ifdef USE_MPI // Make sure everyone is at the same time
      MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
      MPI_Allreduce(MPI_IN_PLACE, &Dt, 1, MPI_SCALAR, MPI_MIN, MPI_COMM_WORLD);
#endif
      if(Dt<1e-14){ printf("Next time step is too small (%e<1e-14). Exiting at step %7i and time %e.\n",Dt,n,T); exit(1);}
      if(Dt!=Dt){ printf("Time step is NaN. Exiting at step %7i and time %e.\n",n,T); exit(1);}
      if     (Dt>(Tf  -T)){ DtCFL = Dt; Dt = Tf  -T; output = true; done = true;}
      else if(Dt>(Tout-T)){ DtCFL = Dt; Dt = Tout-T; output = true;}
      //printf("current time=%e, this Dt=%e, next output at %e\n",T+Dt,Dt,Tout);
      /* Dt = 1e-7; */
      /* if ((n+1)%100==0){output=true;} */
      /* else {output=false;} */
      /* if ((n+1)==1000) {done = true;} */
    }
    timers.stop_timer(4);

    // Advect the particles over that delta t step using the velocity at
    // the current time step
    if(particles.haveParticles()){particles.advectParticles(Dt, arch(U));}
    
    // Us = U
    blasCopy(N_F*N_s*N_E, arch(U), 1, _Us, 1);
    for(int k = 0; k < _order; k++){
      // Ustar = Us + beta*DU
      blasCopy(N_F*N_s*N_E, _Us, 1, _Ustar, 1);           // make Ustar = Us;
      if(k>0){blasAxpy(N_s*N_F*N_E, _beta[k], _DU, 1, _Ustar, 1);} // do Ustar.add(DU,beta[k]) if k=0, beta is zero so this isn't necessary

      Tstar = T + _beta[k]*Dt;

      // Communications for Ustar
      communicator.CommunicateGhosts(N_F, _Ustar);

      //Limit the solution if you so want to do so
      if(k>0){
	if      (Limiter.getLimitingMethod()==1) Limiter.HRlimiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==2) Limiter.M2limiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==3) Limiter.HRIlimiting(communicator, sensor, _Ustar);
	else if (Limiter.getLimitingMethod()==4) Limiter.M2Ilimiting(communicator, sensor, _Ustar);
	else if (Limiter.getLimitingMethod()==5) Limiter.PRlimiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==6) Limiter.PRIlimiting(communicator, sensor, _Ustar);
      }

      // Now you have to calculate f(Ustar)
      dgsolver.dg_solver(_Ustar,_f);
	
      // Solve: DU = Dt*Minv*f(Ustar)
      timers.start_timer(3);
      Lsolver(N_s, N_E, Dt, _Minv, _f, _DU);
      timers.stop_timer(3);
      
      // if 0-order average the solution in the cells
      if (order0){Laverage_cell_p0(N_s, N_E, _DU);}

      // U = U + gamma*DU
      blasAxpy(N_s*N_F*N_E, _gamma[k], _DU, 1, arch(U), 1); // do U.add(DU,gamma[k])

    }// end loop on k
      
    // Communicate and limit solution
    if      (Limiter.getLimitingMethod()==1){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRlimiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==2){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2limiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==3){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRIlimiting(communicator, sensor, arch(U));}
    else if (Limiter.getLimitingMethod()==4){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2Ilimiting(communicator, sensor, arch(U));}
    else if (Limiter.getLimitingMethod()==5){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRlimiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==6){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRIlimiting(communicator, sensor, arch(U));}
        
    T = T + Dt; // update current time
    n++;        // update the time step counter

    // Output the solution
    if(output){
      timers.start_timer(2);
      
      if(myid==0){printf("Solution written to file at step %7i and time %e (current CFL time step:%e).\n",n,T,DtCFL);}
      printer.print(arch(U), particles, count, T);
      printer.print_sensor(sensor, count, T);

      Tout = T + DtOut; // update the new output time
      count++;

      // Output conservation of the fields (if wanted). Works only
      // when data is on host (really this is just for small
      // validation runs)
#ifdef CONS
      dgsolver.conservation(h_U,T);
#endif

      timers.stop_timer(2);
    }// end output steps
  }// end loop on time
  timers.stop_timer(1);
    
  // Free some stuff
  del(_Us);
  del(_Ustar);
  del(_DU);
  del(_UPA);
  del(_f);
  del(_Minv);
#ifdef USE_GPU
  cudaFree(d_U);
#endif
} //end main RK function


scalar RK::DtFromCFL(const int N_s, const int N_E, const scalar CFL, scalar* U, scalar* UPA){
  /*!
    \brief Get the next time step
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] CFL CFL number
    \param[in] U solution to integrate in time
    \param[out] UPA array of |u+a|
    \return next time step
  */
  LfindUPA(N_s, N_E, U, UPA);
  int maxUPAIdx = blasIamax(N_E*N_s,UPA,1)-1; // Fortran starts numbering at 1
#ifdef USE_CPU
  scalar maxUPA = UPA[maxUPAIdx];
  scalar Dt = CFL/maxUPA;
#elif USE_GPU
  scalar* maxUPA = new scalar[1];
  cudaMemcpy(maxUPA, &UPA[maxUPAIdx], sizeof(scalar), cudaMemcpyDeviceToHost);
  scalar Dt = CFL/maxUPA[0];
  delete[] maxUPA;
#endif
  return Dt;
}
