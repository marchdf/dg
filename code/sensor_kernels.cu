/*!
  \file sensor_kernels.cu
  \brief Kernels used by the SENSOR class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup sensor
*/
#include <sensor_kernels.h>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void calc_sensors(int N_E, int N_N, bool sensor1, scalar thresh1, bool sensor2, scalar thresh2, bool sensor3, scalar thresh3, int* neighbors, scalar* Uavg, int* sensors){
  /*!
    \brief Calculate the sensors kernel
    \param[in] N_E number of elements
    \param[in] N_N number of neighbors
    \param[in] sensor1 true if using first sensor
    \param[in] thresh1 threshold value for first sensor
    \param[in] sensor2 true if using second sensor
    \param[in] thresh2 threshold value for second sensor
    \param[in] sensor3 true if using third sensor
    \param[in] thresh3 threshold value for third sensor
    \param[in] neighbors array of neighbor element indices
    \param[in] Uavg average of solution in element
    \param[out] sensors array to hold the sensor
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
#endif

    // Check only if it hasn't been flagged yet
    if(sensors[e] == 0){

      // Initialize
      scalar rhoL,vxL,vyL=0,alphaL,gammaL,betaL=0,pinfL=0,EtL,pL,iL,aL,HL;
      scalar rhoR,vxR,vyR=0,alphaR,gammaR,betaR=0,pinfR=0,EtR,pR,iR,aR,HR;
      scalar RT, vx, vy=0, alpha, gamma, a, i, H, drho, dp, dV2,G,PHI,PSI,W;

      // Get variables for this element (call it left)
      int fc=0;
      rhoL   = Uavg[e*N_F+fc];       fc++;
      vxL    = Uavg[e*N_F+fc]/rhoL;  fc++;
#ifdef TWOD
      vyL    = Uavg[e*N_F+fc]/rhoL;  fc++;
#endif
      EtL    = Uavg[e*N_F+fc];       fc++;
#ifdef PASSIVE
      alphaL = 1/(constants::GLOBAL_GAMMA-1);
#else
      alphaL = Uavg[e*N_F+fc];       fc++;
#endif
#ifdef STIFFENED
      betaL  = Uavg[e*N_F+fc];       fc++;
#endif 
#ifdef GAMCONS
      alphaL = alphaL/rhoL;
#endif
      gammaL = 1.0+1.0/alphaL;
      pinfL = (1-1.0/gammaL)*betaL;
      pL = (gammaL-1)*(EtL - betaL - 0.5*rhoL*(vxL*vxL+vyL*vyL));
      iL = pL*alphaL;
      aL = sqrt((gammaL*(pL+pinfL))/rhoL);
      HL = (EtL + pL)/rhoL;
    
      // Loop on neighbors
      bool done_detecting = false;
      for(int nn = 0; (nn < N_N) && (!done_detecting) ; nn++){
	int right = neighbors[e*N_N+nn]; // neighbor index
	if (right >= 0){

	  // Get variables for this neighbor (call it right)
	  fc=0;
	  rhoR   = Uavg[right*N_F+fc];       fc++;
	  vxR    = Uavg[right*N_F+fc]/rhoR;  fc++;
#ifdef TWOD
	  vyR    = Uavg[right*N_F+fc]/rhoR;  fc++;
#endif
	  EtR    = Uavg[right*N_F+fc];       fc++;
#ifdef PASSIVE
	  alphaR = 1/(constants::GLOBAL_GAMMA-1);
#else
	  alphaR = Uavg[right*N_F+fc];       fc++;
#endif
#ifdef STIFFENED
	  betaR  = Uavg[right*N_F+fc];       fc++;
#endif 
#ifdef GAMCONS
	  alphaR = alphaR/rhoR;
#endif
	  gammaR = 1.0+1.0/alphaR;
	  pinfR = (1-1.0/gammaR)*betaR;
	  pR = (gammaR-1)*(EtR - betaR - 0.5*rhoR*(vxR*vxR+vyR*vyR));
	  iR = pR*alphaR;
	  aR = sqrt((gammaR*(pR+pinfR))/rhoR);
	  HR = (EtR + pR)/rhoR;

	  // Roe averages
	  RT    = sqrt(rhoR/rhoL);
	  vx    = (vxL+RT*vxR)/(1+RT);
	  vy    = (vyL+RT*vyR)/(1+RT);
	  alpha = (alphaL+RT*alphaR)/(1+RT);
	  gamma = 1+1.0/alpha;
	  H     = (HL+RT*HR)/(1+RT);
	  a     = sqrt((gamma-1)*(H-0.5*(vx*vx+vy*vy)));
	  i     =  (iL+RT*iR)/(1+RT);
      
	  // First sensor (contact sensor from Sreenivas)
	  if(sensor1){
	    // Wave strength for contact
	    drho = rhoR - rhoL;
	    dp   = (gamma-1)*(gamma-1)*(alpha*(iR-iL) - i*(alphaR-alphaL));
	    dV2 = drho - dp/(a*a);

	    G = fabs(dV2)/(rhoL+rhoR);
	    PHI = 2*G/((1+G)*(1+G));
	    if(PHI>thresh1){
	      sensors[e]     = 1;
	      sensors[right] = 1;  // this is basically a race condition (another thread could be messing with this too)
	      done_detecting = true;
	    }
	  } // sensor 1

	  // Second sensor (shock sensor from Sreenivas)
	  if((sensor2)&&(!done_detecting)){
	    if(((vxL-aL > vx-a) && (vx-a > vxR-aR)) || ((vyL-aL > vy-a) && (vy-a > vyR-aR))){
	      PSI = fabs(pR-pL)/(pL+pR);
	      W   = 2*PSI/((1+PSI)*(1+PSI));
	      if(W>thresh2){
		sensors[e]     = 2;
		sensors[right] = 2; 
		done_detecting = true;
	      }
	    }
	  } // sensor 2

	  // Third sensor (contact sensor based on gamma)
	  if((sensor3)&&(!done_detecting)){
	    G   = fabs(gammaR-gammaL)/(gammaL+gammaR);
	    PHI = 2*G/((1+G)*(1+G));
	    if(PHI>thresh3){
	      sensors[e]     = 3;
	      sensors[right] = 3; 
	      done_detecting = true;
	    }
	  } // sensor 3
	}

	// If right < 0 then we are at a interesting boundary (and we
	// turn on the sensor)
	else{sensors[e] = 1; done_detecting = true;}
      
      } // loop on neighbors
    } // if sensor is zero
  } // loop on elements
}

arch_global void copy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U){
  /*!
    \brief Copy only the elements which were detected with a sensor and limited.
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] sensors array to hold the sensor
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
    \section Description
    In GPU mode, launches N_E blocks of N_s x N_F x 1 threads. 
  */

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    int sen = sensors[e];
    if (sen != 0){
      for(int i = 0; i < N_s; i++){
	for(int fc = 0; fc < N_F; fc++){
#elif USE_GPU
  int e = blockIdx.x;
  if (e < N_E){
    int sen = sensors[e];
    if (sen != 0){
      int i = threadIdx.x;{
	int fc= threadIdx.y;{
#endif  
	  U[(e*N_F+fc)*N_s+i] = Uold[(e*N_F+fc)*N_s+i];
	} // loop on fields
      } // loop on nodes
    } // if condition on sensor
  } // loop on element
}

  
//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C"
void Lcalc_sensors(int N_E, int N_N, bool sensor1, scalar thresh1, bool sensor2, scalar thresh2, bool sensor3, scalar thresh3, int* neighbors, scalar* Uavg, int* sensors){
  /*!
    \brief Host C function to launch calc_sensors kernel.
    \param[in] N_E number of elements
    \param[in] N_N number of neighbors
    \param[in] sensor1 true if using first sensor
    \param[in] thresh1 threshold value for first sensor
    \param[in] sensor2 true if using second sensor
    \param[in] thresh2 threshold value for second sensor
    \param[in] sensor3 true if using third sensor
    \param[in] thresh3 threshold value for third sensor
    \param[in] neighbors array of neighbor element indices
    \param[in] Uavg average of solution in element
    \param[out] sensors array to hold the sensor
    \section Description
    In GPU mode, launches N_E/blkE blocks of 1 x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(1,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  calc_sensors arch_args (N_E,N_N,sensor1,thresh1,sensor2,thresh2,sensor3,thresh3,neighbors,Uavg,sensors);
};

extern "C"
void Lcopy_detected(int N_s, int N_E, int* sensors, scalar* Uold, scalar* U){
  /*!
    \brief Host C function to launch copy_detected kernel.
    \param[in] N_s number of nodes per element.
    \param[in] N_E number of elements
    \param[in] sensors array to hold the sensor
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
    \section Description
    In GPU mode, launches N_E blocks of N_s x N_F x 1 threads. 
  */
#ifdef USE_GPU
  dim3 dimBlock(N_s,N_F,1);
  dim3 dimGrid(N_E,1);
#endif

  copy_detected arch_args (N_s, N_E, sensors, Uold, U);
};
