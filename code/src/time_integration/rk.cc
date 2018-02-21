/*!
  \file rk.cc
  \brief Runge-Kutta time integration function definitions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup rk
*/
#include "rk.h"

void RK::RK_integration(scalar CFL, scalar VNN, scalar VNNAD, int restart_step,
			int N_E, int N_s, int N_G, int M_T, int M_s, int N_ghosts, int N_N,
			scalar* h_Minv, 
			scalar* h_U, int* neighbors,
			Limiting &Limiter, bool order0, DG_SOLVER &dgsolver, COMMUNICATOR &communicator, PRINTER &printer, SENSOR &sensor, TIMERS &timers, MEM_COUNTER &mem_counter, LAGRANGE_PARTICLES &particles){
  /*!
    \brief Main RK integration function
    \param[in] CFL CFL number
    \param[in] VNN VNN number for viscous physics
    \param[in] VNNAD The VNN number for AD treatment
    \param[in] restart_step output step for a restart
    \param[in] N_E number of elements
    \param[in] N_s number of nodes per element
    \param[in] N_G number of gaussian nodes per element
    \param[in] M_T number of interfaces
    \param[in] M_s number of nodes per interface
    \param[in] N_ghosts number of ghost elements
    \param[in] N_N number of sides (neighbors) per element
    \param[in] h_Minv host array containing the inverse mass matrices of each element
    \param[out] h_U solution to integrate in time
    \param[in] neighbors the neighbor array used for sensor and limiting
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
  scalar* _superK; //PEJ 11/09/2017
  scalar* _Us;
  scalar* _Ustar;
  scalar* _DU;
  scalar* _UPA; 
  scalar* _Mew; //PEJ 06/01/2017
  scalar* _f;
  scalar* _Minv;
  int* _neighbors; //PEJ 10/24/2017
  int* _SensorTag_Aug; //PEJ 11/29/2017: Holds sensor values for all elements (flesh+ghost)
  scalar* _LamMax_local;
  scalar* _DivMax_local;
  scalar* _CsMax_local;
  scalar* _LamMax; //AD parameter, global maxima
  scalar* _DivMax; //AD parameter, global maxima
  scalar* _CsMax; //AD parameter, global maxima
  
#ifdef USE_CPU
  _superK = new scalar[N_s*N_E*N_F*_steps]; makeZero(_superK, N_s*N_E*N_F*_steps);        mem_counter.addToCPUCounter(N_s*N_E*N_F*_steps*sizeof(scalar));
  _Us    = new scalar[N_s*N_E*N_F];  makeZero(_Us   ,N_s*N_E*N_F);	                  mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _Ustar = new scalar[N_s*(N_E+N_ghosts)*N_F];  makeZero(_Ustar,N_s*(N_E+N_ghosts)*N_F);  mem_counter.addToCPUCounter(N_s*(N_E+N_ghosts)*N_F*sizeof(scalar));
  _DU    = new scalar[N_s*N_E*N_F];  makeZero(_DU   ,N_s*N_E*N_F);                        mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _UPA   = new scalar[N_s*N_E];      makeZero(_UPA  ,N_s*N_E);                            mem_counter.addToCPUCounter(N_s*N_E*sizeof(scalar));
  _Mew   = new scalar[N_s*N_E];      makeZero(_Mew  ,N_s*N_E);                            mem_counter.addToCPUCounter(N_s*N_E*sizeof(scalar));
  _f     = new scalar[N_s*N_E*N_F];  makeZero(_f    ,N_s*N_E*N_F);                        mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));
  _Minv  = new scalar[N_s*N_s*N_E];  memcpy(_Minv, h_Minv, N_s*N_s*N_E*sizeof(scalar));   mem_counter.addToCPUCounter(N_s*N_s*N_F*sizeof(scalar));
  _neighbors = new int[N_E*N_N]; memcpy(_neighbors, neighbors, N_N*N_E*sizeof(int)); mem_counter.addToCPUCounter(N_E*N_N*sizeof(int));
  _SensorTag_Aug = new int[N_E+N_ghosts];                                  mem_counter.addToCPUCounter((N_E+N_ghosts)*sizeof(int));
  _LamMax = new scalar[D]; makeZero(_LamMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
  _DivMax = new scalar[D]; makeZero(_DivMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
  _CsMax = new scalar[D]; makeZero(_CsMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
   _LamMax_local = new scalar[D]; makeZero(_LamMax_local, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
  _DivMax_local = new scalar[D]; makeZero(_DivMax_local, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
  _CsMax_local = new scalar[D]; makeZero(_CsMax_local, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
#elif USE_GPU
  scalar* d_U;
  // Allocate on device
  cudaMalloc((void**) &_superK   , N_s*N_E*N_F*_steps*sizeof(scalar));            mem_counter.addToGPUCounter(N_s*N_E*N_F*_steps*sizeof(scalar));
  cudaMalloc((void**) &d_U   , N_s*N_E*N_F*sizeof(scalar));                       mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));
  cudaMalloc((void**) &_Us   , N_s*N_E*N_F*sizeof(scalar));                       mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_Ustar, N_s*N_E*N_F*sizeof(scalar));  		       	  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar)); 
  cudaMalloc((void**) &_DU   , N_s*N_E*N_F*sizeof(scalar));		       	  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_UPA  , N_s*N_E*sizeof(scalar));     			  mem_counter.addToGPUCounter(N_s*N_E*sizeof(scalar));		      
  cudaMalloc((void**) &_Mew  , N_s*N_E*sizeof(scalar));     		      	  mem_counter.addToGPUCounter(N_s*N_E*sizeof(scalar));		      
  cudaMalloc((void**) &_f    , N_s*N_E*N_F*sizeof(scalar));			  mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));	      
  cudaMalloc((void**) &_Minv , N_s*N_s*N_E*sizeof(scalar));		      	  mem_counter.addToGPUCounter(N_s*N_s*N_E*sizeof(scalar));
  cudaMalloc((void**) &_neighbors , N_N*N_E*sizeof(int));			  mem_counter.addToGPUCounter(N_N*N_E*sizeof(int));
  cudaMalloc((void**) &_SensorTag_Aug , (N_E+N_ghosts)*sizeof(int));		  mem_counter.addToGPUCounter((N_E+N_ghosts)*sizeof(int));
  cudaMalloc((void**) &_LamMax, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));
  cudaMalloc((void**) &_DivMax, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));
  cudaMalloc((void**) &_CsMax, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));
   cudaMalloc((void**) &_LamMax_local, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));
  cudaMalloc((void**) &_DivMax_local, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));
  cudaMalloc((void**) &_CsMax_local, D*sizeof(scalar));                                mem_counter.addToGPUCounter(D*sizeof(scalar));

  // Set to zero
  cudaMemset(_superK   , (scalar)0.0, N_s*N_E*N_F*_steps*sizeof(scalar));
  cudaMemset(_Us   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_Ustar, (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_DU   , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_UPA  , (scalar)0.0, N_s*N_E*sizeof(scalar));
  cudaMemset(_f    , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
  cudaMemset(_SensorTag_Aug    , (int)0, (N_E+N_ghosts)*sizeof(int));
  cudaMemset(_LamMax, (scalar)0.0, D*sizeof(scalar));
  cudaMemset(_DivMax, (scalar)0.0, D*sizeof(scalar));
  cudaMemset(_CsMax, (scalar)0.0 , D*sizeof(scalar));
  cudaMemset(_LamMax_local, (scalar)0.0, D*sizeof(scalar));
  cudaMemset(_DivMax_local, (scalar)0.0, D*sizeof(scalar));
  cudaMemset(_CsMax_local, (scalar)0.0 , D*sizeof(scalar));

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
  double Tout = _output_time_array[0];   // next output time
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
      if(particles.haveParticles()){particles.advectParticles(0.0, arch(U));} // get the particle locations by advecting them by dt = 0
      printf("particle print call\n"); fflush(stdout);
      printer.print(arch(U), particles, count, T);

      //PEJ 10/24/2017: initialize sensor values and global Cparam maxima
#ifdef RADSINGLEFLUID
      printf("pr=%d: Attempting initial sensing in RK\n",myid); fflush(stdout);
      communicator.CommunicateGhosts(N_F, arch(U));
      printf("Calling sensor for pr=%d\n",myid);
      sensor.sensing(_neighbors, arch(U));

      //Put sensor results in the SensorTag_Aug array:
      for (int e = 0; e < N_E; e++) {
	_SensorTag_Aug[e] = sensor.get1Sensor(e); }
      //Treat ghost elements in SensorTag_Aug array:
      communicator.CommunicateSensor(_SensorTag_Aug);

      /*
      printf("The initial sensor array:\n");
      for (int e = 0; e < N_E; e++)
	{
	  printf("sensor[%d] = %d\n", e, sensor.getSensors()[e]);
	}
      */
#endif //end RADSINGLEFLUID
      //exit(1);
      
      // Limit the initial solution before integrating to avoid problems
      printf("pr=%d: Limiting initial condition if limiter prescribed\n",myid); fflush(stdout);
      if      (Limiter.getLimitingMethod()==1){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRlimiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==2){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2limiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==3){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRIlimiting(communicator, sensor, arch(U));}
      else if (Limiter.getLimitingMethod()==4){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2Ilimiting(communicator, sensor, arch(U));}
      else if (Limiter.getLimitingMethod()==5){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRlimiting(communicator, arch(U));}
      else if (Limiter.getLimitingMethod()==6){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRIlimiting(communicator, sensor, arch(U));}
      else if (Limiter.getLimitingMethod()==7){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.P0Ilimiting(communicator, sensor, dgsolver, arch(U));}
      printer.print_sensor(sensor, count, T);
  


      printf("pr=%d: The output time array:\n",myid);
      for (int j = 0; j < _output_time_array.size(); j++)
	{
	  printf("output time[%d] = %f\n", j, _output_time_array[j]);
	}

      // Output conservation of the fields (if wanted). Works only
      // when data is on host (really this is just for small
      // validation runs)
#ifdef CONS
      dgsolver.conservation(h_U,T);
#endif

      //PEJ 10/01/2017: output taylor-green vortex statistics.
      //Only input from here is time, it takes the rest from
      //the dgsolver class
#ifdef TGVSTATS
      dgsolver.TGV_statistics(arch(U),T);
#endif
      
      timers.stop_timer(2);
    }
  }
  count++;
  Tout = _output_time_array[count];

  // Time integration
  timers.start_timer(1);
  printf("pr=%d: Entering the thunderdome for time integration:\n",myid); fflush(stdout);
  while (!done){
  //  while (n < 1){
    
    // Give the deck a negative CFL for fixed time step and no output
    timers.start_timer(4);
    if(CFL<0){
      Dt = _output_time_array[1] - _output_time_array[0]; output = false;
      if(Dt>(_Tf  -T)){Dt = _Tf  -T; done = true;}
    }
    else{
      // Find new Dt
      //Dt = DtFromCFL(N_s, N_E, CFL, arch(U),_UPA); output = false;
#ifndef RADSINGLEFLUID
      Dt = DtFromCFL_and_VNN(N_s, N_E, CFL, VNN, arch(U),_UPA,_Mew); output = false;
#endif
#ifdef RADSINGLEFLUID
      Dt = DtFromCFL_and_VNN_and_VNNAD(N_s, N_E, CFL, VNN, VNNAD, arch(U),_UPA,_Mew, dgsolver.get_DivMax(), dgsolver.get_ADepsMax(), dgsolver.get_Beta_S(), dgsolver.get_Mew_S()); output = false;
#endif
#ifdef USE_MPI // Make sure everyone is at the same time
      MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
      MPI_Allreduce(MPI_IN_PLACE, &Dt, 1, MPI_SCALAR, MPI_MIN, MPI_COMM_WORLD); //identify global Dt value
#endif
      if(Dt<1e-14){ printf("Next time step is too small (%e<1e-14). Exiting at step %7i and time %e.\n",Dt,n,T); exit(1);}
      //if(Dt<1e-17){ printf("Next time step is too small (%e<1e-14). Exiting at step %7i and time %e.\n",Dt,n,T); exit(1);}
      if(Dt!=Dt  ){ printf("Time step is NaN. Exiting at step %7i and time %e.\n",n,T); exit(1);}
      if     (Dt>(_Tf  -T)){ DtCFL = Dt; Dt = _Tf  -T; output = true; done = true; }
      else if(Dt>(Tout -T)){ DtCFL = Dt; Dt = Tout -T; output = true; }
      // printf("current time=%e, this Dt=%e, next output at %e\n",T+Dt,Dt,Tout); fflush(stdout);
      /* Dt = 1e-7; */
      /* if ((n+1)%100==0){output=true;} */
      /* else {output=false;} */
      /* if ((n+1)==1000) {done = true;} */
    }
    timers.stop_timer(4);

    // Advect the particles over that delta t step using the velocity at
    // the current time step
    //printf("Preparing to advect particles, maybe\n"); fflush(stdout);
    if(particles.haveParticles()){particles.advectParticles(Dt, arch(U));}
    
    // Us = U
    //printf("Preparing blasCopy between arch(U) and _Us\n"); fflush(stdout);
    blasCopy(N_F*N_s*N_E, arch(U), 1, _Us, 1); //Inside k loop, _Us does not change; it is the starter value used to form the predictor, _Ustar
    //for(int k = 0; k < _order; k++){
    for(int k = 0; k < _steps; k++){
      //  printf("Entered RK loop, k=%d/%d\n", k, _steps); fflush(stdout);
      // Ustar = Us + beta*DU
      //printf("Another blasCopy, this one between _Us and _Ustar\n"); fflush(stdout);
      /*
      printf("Copying Us to Ustar, kRK=%d\n",k);
      for (int e = 0; e < N_E; e++)
	{
	  for (int i = 0; i < N_s; i++)
	    {
	      printf("e=%d,node=%d:, present _Us = %f, present _Ustar = %f\n",e,i, _Us[(e*N_F+0)*N_s+i], _Ustar[(e*N_F+0)*N_s+i]);
	    }
	}
      */
      //Update the predictor solution, _Ustar:
      if (_order == 1 || _order == 4)
	{ //Only need _beta[k] for the Ustar update
	  blasCopy(N_F*N_s*N_E, _Us, 1, _Ustar, 1);           // make Ustar = Us;
	  if(k>0){
	    //printf("For k>0: performing blas Axpy\n"); fflush(stdout);
	    blasAxpy(N_s*N_F*N_E, _beta[k], _DU, 1, _Ustar, 1);
	  } // do Ustar.add(DU,beta[k]) if k=0, beta is zero so this isn't necessary
	}
      else
	{
	  //RK[8/13] or other methods with non-diagonal alpha matrix
	  //_Us tracks the master solution, while _Ustar is the predictor that dgsolver works with
	  blasCopy(N_F*N_s*N_E, _Us, 1, _Ustar, 1);           // make Ustar = Us;
	  //general RK routine: need the full _alpha matrix multiplied by previous superK values:
	  if(k>0)
	    {
	      //Update the predictor, Ustar. This is not happening in parallel-friendly format.
	      //Set Ustar = 0:
	      for (int slot = 0; slot < N_E*N_s*N_F; slot++)
		{ 
		  //_Ustar[slot] = 0.0;
		  //Dt multplication happens in the Lsolver routine for Marc's approach,
		  //so it is already contained in superK
		  for (int j = 0; j < k; j++) //previous RK contribution
		    {
		      _Ustar[slot] += _alpha[(k-1)*_steps + j] * _superK[slot*_steps + j];
		    }
		}
	      /*
	      //scalar* _alphaLocal = new scalar[_steps];
	      for (int j = 0; j < _steps; j++)
		{
		  //alphaLocal[j] = _alpha[k-1][j];
		  scalar alphaLocal = _alpha[k-1][j];
		  blasAxpy(N_s*N_F*N_E, alphaLocal, _DU, 1, _Ustar, 1);
		}
	      //With alphaLocal populated, use blasAxpy to
	      //add contribution to predictor solution
	      blasAxpy(N_s*N_F*N_E, _beta[k], _DU, 1, _Ustar, 1);
	      //delete[] _alphaLocal;
	      */
	    } //end k>0 case for _order outside {1,4}
	    
	}
      Tstar = T + _beta[k]*Dt;

      // Communications for Ustar
      communicator.CommunicateGhosts(N_F, _Ustar);

      //Limit the predictor solution if you so want to do so
      if(k>0){
	//PEJ 10/24/2017: get sensor values.
#ifdef RADSINGLEFLUID
	//	printf("Again attempting sensing withn RK substepping, on Ustar\n"); fflush(stdout);
	sensor.sensing(_neighbors, _Ustar);
	//Re-populate _SensorTag_Aug and treat ghost elements
	for (int e = 0; e < N_E; e++) {
	  _SensorTag_Aug[e] = sensor.get1Sensor(e); }
	//Treat ghost elements in SensorTag_Aug array:
	//communicator.CommunicateSensor(_SensorTag_Aug);
	//Don't grab max AD parameters here because I need to make a trip through residual calculation first
#endif
	if      (Limiter.getLimitingMethod()==1) Limiter.HRlimiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==2) Limiter.M2limiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==3) Limiter.HRIlimiting(communicator, sensor, _Ustar);
	else if (Limiter.getLimitingMethod()==4) Limiter.M2Ilimiting(communicator, sensor, _Ustar);
	else if (Limiter.getLimitingMethod()==5) Limiter.PRlimiting(communicator, _Ustar);
	else if (Limiter.getLimitingMethod()==6) Limiter.PRIlimiting(communicator, sensor, _Ustar);
	else if (Limiter.getLimitingMethod()==7) Limiter.P0Ilimiting(communicator, sensor, dgsolver, _Ustar);
      }

      // Now you have to calculate f(Ustar)
      //k is the RK substep we are presently on
      //   printf("Calling dg_solver from within RK loop\n"); fflush(stdout);
	// dgsolver.dg_solver(_Ustar,_f, k);
      //dgsolver.dg_solver(_Ustar,_f, k, sensor.getSensors());
      dgsolver.dg_solver(_Ustar,_f, k, _SensorTag_Aug, _LamMax, _DivMax, _CsMax, myid);

      // Solve: DU = Dt*Minv*f(Ustar)
      timers.start_timer(3);
      //printf("Performing Lsolver call\n"); fflush(stdout);
      Lsolver(N_s, N_E, Dt, _Minv, _f, _DU);
      timers.stop_timer(3);
      
      // if 0-order average the solution in the cells
      if (order0){Laverage_cell_p0(N_s, N_E, _DU);}

      // U = U + gamma*DU (this is the master solution update)
      //   printf("Another Axpy for _Du and arch(U)\n"); fflush(stdout);
      blasAxpy(N_s*N_F*N_E, _gamma[k], _DU, 1, arch(U), 1); // do U.add(DU,gamma[k])

      if (_order == 8)
	{
	  //RK8 case, and other cases with non-diagonal alpha: need
	  //to store _DU as _superK
	  for (int slot = 0; slot < N_E*N_s*N_F; slot++){
	    _superK[slot*_steps + k] = _DU[slot]; }
	}

    }// end loop on k
      
    // Communicate and limit solution
    //12/11/2017: For the RADSINGLEFLUID case,
    //there is no need to calculate the sensor here;
    //sensor is populated before each residual evaluation
    //in the k-loop.
    if      (Limiter.getLimitingMethod()==1){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRlimiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==2){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2limiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==3){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.HRIlimiting(communicator, sensor, arch(U));}
    else if (Limiter.getLimitingMethod()==4){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.M2Ilimiting(communicator, sensor, arch(U));}
    else if (Limiter.getLimitingMethod()==5){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRlimiting(communicator, arch(U));}
    else if (Limiter.getLimitingMethod()==6){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.PRIlimiting(communicator, sensor, arch(U));}
    else if (Limiter.getLimitingMethod()==7){ communicator.CommunicateGhosts(N_F, arch(U)); Limiter.P0Ilimiting(communicator, sensor, dgsolver, arch(U));}
        
    T = T + Dt; // update current time
    n++;        // update the time step counter
    //PEJ 10/01/2017: output taylor-green vortex statistics.
    //Only input from here is time, it takes the rest from
    //the dgsolver class
    //Take this into case structure below to only run whenver
    //solution is output to .pos files
#ifdef TGVSTATS
      dgsolver.TGV_statistics(arch(U),T);
#endif
    // Output the solution
    if(output){
      timers.start_timer(2);
      //printf("Attempting output\n"); fflush(stdout);
      //if(myid==0){printf("Solution written to file at step %7i and time %e (current CFL time step:%e).\n",n,T,DtCFL);}
#ifdef RADSINGLEFLUID
      //need to calculate and communicate sensor
      //for the printer's sake
      communicator.CommunicateGhosts(N_F, arch(U));
      sensor.sensing(_neighbors, arch(U));
      //Re-populate _SensorTag_Aug and treat ghost elements
      for (int e = 0; e < N_E; e++) {
	_SensorTag_Aug[e] = sensor.get1Sensor(e); }
      //For printer case, don't worry about ghost elements
#endif
      if(myid==0){printf("Solution written to file at step %7i and time %e (current CFL-VNN-VNNAD time step:%e).\n",n,T,DtCFL);}
      printer.print(arch(U), particles, count, T);
      printer.print_sensor(sensor, count, T);

      count++;
      Tout = _output_time_array[count]; // update the new output time
      
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
  del(_superK); //PEJ 11/09/2017
  del(_Us);
  del(_Ustar);
  del(_DU);
  del(_UPA);
  del(_Mew); //PEJ 06/01/2017
  del(_f);
  del(_Minv);
  del(_neighbors); //PEJ 10/24/2017
  del(_SensorTag_Aug);
  del(_LamMax);
  del(_DivMax);
  del(_CsMax);
  del(_LamMax_local);
  del(_DivMax_local);
  del(_CsMax_local);

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
//PEJ Edit 06/01/2017: Subroutine to get the timestep for advection-diffusion problem
scalar RK::DtFromCFL_and_VNN(const int N_s, const int N_E, const scalar CFL, scalar VNN, scalar* U, scalar* UPA, scalar* Mew){
  /*!
    \brief Get the next time step
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] CFL CFL number
    \param[in] VNN Von-Neumann number
    \param[in] U solution to integrate in time
    \param[out] UPA array of |u+a|
    \param[out] array of viscosity at each location
    \return next time step
  */
  //Note: in main, the VNN and CFL have been calculated based on mesh spacing, so there
  //is no need to account for delx here
  //printf("Entered
  LfindUPA(N_s, N_E, U, UPA);
  LfindMew(N_s, N_E, U, Mew); 
  int maxUPAIdx = blasIamax(N_E*N_s,UPA,1)-1; // Fortran starts numbering at 1
  int maxMewIdx = blasIamax(N_E*N_s,Mew,1)-1; // Fortran starts numbering at 1
#ifdef USE_CPU
  scalar maxUPA = UPA[maxUPAIdx];
  scalar maxMew = Mew[maxMewIdx];
  //  printf("max mew = %f, max UPA = %f\n",maxMew,maxUPA);
  //printf("maxUPA = %f, maxMew = %f, CFL/maxUPA = %f, VNN/maxMew = %f\n",maxUPA,maxMew,CFL/maxUPA,VNN/maxMew);
  scalar Dt = 1.0 / (maxUPA/CFL + maxMew/VNN);
#elif USE_GPU
  scalar* maxUPA = new scalar[1];
  scalar* maxMew = new scalar[1];
  cudaMemcpy(maxUPA, &UPA[maxUPAIdx], sizeof(scalar), cudaMemcpyDeviceToHost);
  cudaMemcpy(maxMew, &Mew[maxMewIdx], sizeof(scalar), cudaMemcpyDeviceToHost);
  //scalar Dt = CFL/maxUPA[0];
  scalar Dt = 1.0 / (maxUPA[0]/CFL + maxMew[0]/VNN);
  printf("maxUPA[0]=%f, maxMew[0]=%f, Dt=%f\n", maxUPA[0], maxMew[0], Dt);
  delete[] maxUPA;
  delete[] maxMew;
#endif
  return Dt;
}

//PEJ Edit 06/01/2017: Subroutine to get the timestep for advection-diffusion problem with artificial dissipation
scalar RK::DtFromCFL_and_VNN_and_VNNAD(const int N_s, const int N_E, const scalar CFL, scalar VNN, scalar VNNAD, scalar* U, scalar* UPA, scalar* Mew,  scalar* DivMax, scalar* ADepsMax, scalar Beta_S, scalar Mew_S){
  /*!
    \brief Get the next time step
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] CFL CFL number
    \param[in] VNN Von-Neumann number
    \param[in] Artificial dissipation VNN Von-Neumann number
    \param[in] U solution to integrate in time
    \param[out] UPA array of |u+a|
    \param[out] Mew array of viscosity at each location
    \param[in] DivMax the maximum super-divergence, relayed from DGsolver
    \param[in] APepsMax Sent from main, stored by rk solver, applied here
    \param[in] Beta_S AD global strength parameter, relayed in from DgSolver
    \param[in] Mew_S AD global spread parameter, relayed in from DgSolver
    \return next time step
  */
  //Note: in main, the VNN and CFL have been calculated based on mesh spacing, so there
  //is no need to account for delx here
  //printf("Entered
  LfindUPA(N_s, N_E, U, UPA);
  LfindMew(N_s, N_E, U, Mew); 
  int maxUPAIdx = blasIamax(N_E*N_s,UPA,1)-1; // Fortran starts numbering at 1
  int maxMewIdx = blasIamax(N_E*N_s,Mew,1)-1; // Fortran starts numbering at 1
#ifdef USE_CPU
  scalar maxUPA = UPA[maxUPAIdx];
  scalar maxMew = Mew[maxMewIdx];
  //  printf("max mew = %f, max UPA = %f\n",maxMew,maxUPA);
  //printf("maxUPA = %f, maxMew = %f, CFL/maxUPA = %f, VNN/maxMew = %f\n",maxUPA,maxMew,CFL/maxUPA,VNN/maxMew);
  //Timestep size: combine advective timestep limit, viscous physics timestep, and AD timestep

  //Identify maximum ADeps value (will be same each timestep, could streamline)
  scalar term1 = ADepsMax[0];
  for (int a = 0; a < D; a++)
    {
      term1 = fmax(ADepsMax[a], term1);
    }
  
  //Diffusivity of C equation is first on the list.
  //Throwing in a 2 here for extra safety
  scalar CvisDT = 2.0*Mew_S * term1 * maxUPA; //C diffusivity coefficient includes max wavespeed
  //Now, get the maximum Kappa diffusion for AD on conserved variables
#ifdef ONED
  scalar KappaMax = Beta_S * (term1*term1) * DivMax[0];
#endif
#ifdef TWOD
  scalar KappaMax = Beta_S * (term1*term1) * fmax(DivMax[0], DivMax[1]);
#endif
#ifdef THREED
  scalar KappaMax = Beta_S * (term1*term1) * fmax(fmax(DivMax[0], DivMax[1]), DivMax[2]);
#endif
  //Maximum between KappaMa and C equation diffusity goeverns AD timestep limit.
  //Why can I take minimum? Because Kappa and Cvis act on different field variables
  scalar ADFactor = fmax(KappaMax, CvisDT);
  
  //Now combine upa, mew, and ADFactor to determine the timestep limit
  scalar Dt = 1.0 / (maxUPA/CFL + maxMew/VNN);
#ifdef RADSINGLEFLUID
  Dt = 1.0 / (maxUPA/CFL + maxMew/VNN + ADFactor/VNNAD);
  //printf("maxUPA = %f, maxMew = %f, CvisDT = %f, KappaMax = %f, ADepsMax=%f, Beta_S=%f, DivMax=%f, CFL/maxUPA = %f, VNN/maxMew = %f, VNNAD/ADFactor=%f, Dt=%f\n",maxUPA,maxMew,CvisDT,KappaMax,term1,Beta_S, DivMax[0], CFL/maxUPA,VNN/maxMew,VNNAD/ADFactor,Dt);
#endif
  
#elif USE_GPU
  scalar* maxUPA = new scalar[1];
  scalar* maxMew = new scalar[1];
  cudaMemcpy(maxUPA, &UPA[maxUPAIdx], sizeof(scalar), cudaMemcpyDeviceToHost);
  cudaMemcpy(maxMew, &Mew[maxMewIdx], sizeof(scalar), cudaMemcpyDeviceToHost);
  //scalar Dt = CFL/maxUPA[0];
  scalar Dt = 1.0 / (maxUPA[0]/CFL + maxMew[0]/VNN);
  scalar Dt = 1.0 / (maxUPA[0]/CFL + maxMew[0]/VNN);
  printf("maxUPA[0]=%f, maxMew[0]=%f, Dt=%f\n", maxUPA[0], maxMew[0], Dt);
  delete[] maxUPA;
  delete[] maxMew;
#endif
  return Dt;
}
