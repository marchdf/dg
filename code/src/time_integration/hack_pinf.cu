/*!
  \file hack_pinf.cu
  \brief Kernels used by hack_pinf
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup rk
*/
#include "hack_pinf.h"
#include <cstdlib>
#include <stdio.h>


//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void hack_pinf_20150727(int N_s, int N_E, scalar* U){
  /*!
    \brief Solve for pinf using dubious methods. See notes 2015/07/27
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[out] U the solution (we will modify the pinf field)
  */


#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

    // field properties
    scalar rho   = U[(e*N_F+0)*N_s+i];
    scalar u     = U[(e*N_F+1)*N_s+i]/rho;  // (rho u / rho) = u
    scalar v     = U[(e*N_F+2)*N_s+i]/rho;  // (rho v / rho) = v
    scalar Et    = U[(e*N_F+3)*N_s+i];
    scalar G     = U[(e*N_F+4)*N_s+i];
    scalar gamma = 1+1.0/G;
    scalar beta  = U[(e*N_F+5)*N_s+i];
    scalar pinf  = beta*(gamma-1)/gamma;

    // ND quantities
    scalar rho_air = 1.1765;
    scalar gamma_air = 1.4;
    scalar patm = 101325;
    scalar cs_air = sqrt(gamma_air*patm/rho_air);
    scalar p_ND   = rho_air*cs_air*cs_air;
    
    // mixture
    scalar rv = 0.8;
    scalar alpha_g = rv/(1+rv);
    scalar alpha_l = 1-alpha_g;

    // material properties (non-dimensionalized)
    scalar gamma_l = 5.5;
    scalar pinf_l = 492115000/p_ND;
    scalar gamma_g = 1.4;

    // Newton solver
    scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v)) - gamma*pinf;
    scalar E = Et;
    scalar eps = 1e-16;
    for(int k=0; k<100; k++){

      // evaluate function and derivative
      scalar get_f_p = Et - 0.5*rho*(u*u+v*v) - p/(gamma - 1) - 1/(gamma - 1)*(1.0/(alpha_l/(gamma_l*(p + pinf_l)) + alpha_g/(gamma_g*p)) - gamma*p);
      scalar get_fp_p = 1.0/(gamma-1) + 1.0/(gamma-1) * ( ( -alpha_g/(p*p * gamma_g) - alpha_l/(gamma_l*(p+pinf_l)*(p+pinf_l)) )/( (alpha_g/(p*gamma_g) - alpha_l/(gamma_l*(p+pinf_l))) * (alpha_g/(p*gamma_g) - alpha_l/(gamma_l*(p+pinf_l))))  - gamma);
      
      // get the next guess
      scalar pnew = p - get_f_p/get_fp_p;

      // save data first
      p = pnew;
      pinf = (1.0/(alpha_l/(gamma_l*(p+pinf_l)) + alpha_g/(gamma_g*p)) - gamma*p)/gamma;
      E = p/(gamma-1) + gamma*pinf/(gamma-1) + 0.5*rho*(u*u+v*v) ;
      // print k, p, pinf, E, Ef, '{0:.5e}'.format(fabs(E-Ef))  
      
      // Test for convergence and exit if converged
      scalar get_f_pnew = Et - 0.5*rho*(u*u+v*v) - pnew/(gamma - 1) - 1/(gamma - 1)*(1.0/(alpha_l/(gamma_l*(pnew + pinf_l)) + alpha_g/(gamma_g*pnew)) - gamma*pnew);
      if (fabs(get_f_p-get_f_pnew) < eps){break;}

    }

    // Update the energy and pinf field
    U[(e*N_F+3)*N_s+i] = E;
    U[(e*N_F+5)*N_s+i] = gamma*pinf/(gamma-1);
    
#ifdef USE_CPU
    }
#endif
  }
  
}

//==========================================================================
arch_global void hack_pinf_20150807(int N_s, int N_E, scalar pm, scalar* U){
  /*!
    \brief Solve for pinf using dubious methods. See notes 2015/08/07
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] pm  initial mixture pressure (global variable from initial condition)
    \param[out] U the solution (we will modify the pinf field)
  */


#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

    // field properties
    scalar rho   = U[(e*N_F+0)*N_s+i];
    scalar u     = U[(e*N_F+1)*N_s+i]/rho;  // (rho u / rho) = u
    scalar v     = U[(e*N_F+2)*N_s+i]/rho;  // (rho v / rho) = v
    scalar Et    = U[(e*N_F+3)*N_s+i];
    scalar G     = U[(e*N_F+4)*N_s+i];
    scalar gamma = 1+1.0/G;
    scalar beta  = U[(e*N_F+5)*N_s+i];
    scalar pinf  = beta*(gamma-1)/gamma;

    // ND quantities
    // scalar rho_air = 1.1765;
    // scalar gamma_air = 1.4;
    // scalar patm = 101325;
    // scalar cs_air = sqrt(gamma_air*patm/rho_air);
    scalar p_ND   = 141855.000000; // rho_air*cs_air*cs_air;
    
    // mixture
    scalar rv = 0.8;
    scalar alpha_g = rv/(1+rv);
    scalar alpha_l = 1-alpha_g;

    // material properties (non-dimensionalized)
    scalar gamma_l = 5.5;
    scalar pinf_l = 492115000/p_ND;
    scalar gamma_g = 1.4;

    // redefine the volume fractions to account for the pressure change
    // initial pressure in flow
    //scalar pm = 1.306875;    // Ms = 2
    //scalar pm = 0.580833;    // Ms = 3
    //scalar pm = 0.326719;    // Ms = 4
    //scalar pm = 0.209100;    // Ms = 5
    //scalar pm = 0.145208;    // Ms = 6
    //scalar pm = 0.106684;    // Ms = 7
    //scalar pm = 0.081680;    // Ms = 8
    //scalar pm = 0.064537;    // Ms = 9
    //scalar pm = 0.052275;    // Ms = 10
    scalar gamma_i = 1; // maybe use gamma_g?
    
    // Newton solver
    scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v)) - gamma*pinf;
    scalar E = Et;
    scalar reltol = 1e-13;
    scalar tol    = 1e-13;
    scalar get_f_p, get_f_pnew, get_fp_p, pnew,K,delta = 0;
    for(int k=0; k<1000; k++){
     
      // // evaluate function and derivative
      // get_f_p = Et - 0.5*rho*(u*u+v*v) + (gamma_g*p*((-1 + gamma - gamma_l)*p - gamma_l*pinf_l) + (-1 + gamma - gamma_g)*gamma_l*p*(p + pinf_l)*pow(p/pm,gamma_i)*rv)/((-1 + gamma)*(gamma_g*p + gamma_l*(p + pinf_l)*pow(p/pm,gamma_i)*rv));
      // get_fp_p = (gamma*gamma*(-1 + gamma - gamma_l)*p*p - gamma_g*gamma_l*((2 - 2*gamma + gamma_g + gamma_g*gamma_i + gamma_l - gamma_i*gamma_l)*p*p + (2 - 2*gamma + gamma_g*gamma_i + 2*gamma_l - 2*gamma_i*gamma_l)*p*pinf_l - (-1 + gamma_i)*gamma_l*pinf*pinf_l)*pow(p/pm,gamma_i)*rv + (-1 + gamma - gamma_g)*gamma_l*gamma_l*(p+pinf_l)*(p+pinf_l)* pow(p/pm,2*gamma_i)*rv*rv)/((-1 + gamma)*(gamma_g*p + gamma_l*(p + pinf_l)*pow(p/pm,gamma_i)*rv)*(gamma_g*p + gamma_l*(p + pinf_l)*pow(p/pm,gamma_i)*rv));

      // simplified formula for function and derivative
      get_f_p = gamma_g*p*(Et*(-1 + gamma) + (0.5*rho*(u*u+v*v)) - p + gamma*(-(0.5*rho*(u*u+v*v)) + p) - gamma_l*(p + pinf_l)) - gamma_l*(Et - Et*gamma + (-1 + gamma)*(0.5*rho*(u*u+v*v)) + (1 - gamma + gamma_g)*p)*(p + pinf_l)*pow(p/pm,gamma_i)*rv;
      get_fp_p = (gamma_g*p*(Et*(-1 + gamma) + (0.5*rho*(u*u+v*v)) - gamma*(0.5*rho*(u*u+v*v)) + 2*(-1 + gamma - gamma_l)*p - gamma_l*pinf_l) + gamma_l*(Et*(-1 + gamma)*(p + gamma_i*p + gamma_i*pinf_l) - (-1 + gamma)*(0.5*rho*(u*u+v*v))*(p + gamma_i*p + gamma_i*pinf_l) + (-1 + gamma - gamma_g)*p*((2 + gamma_i)*p + pinf_l + gamma_i*pinf_l))*pow(p/pm,gamma_i)*rv)/p;
	
      // get the next guess
      delta = get_f_p/get_fp_p;
      pnew = p - delta;
     
      // save data first
      p = pnew;

      // get pinf and energy
      K = rv*pow(p/pm,gamma_i); 
      alpha_g = K/(1.0+K);
      alpha_l = 1-alpha_g;
      pinf = (1.0/(alpha_l/(gamma_l*(p+pinf_l)) + alpha_g/(gamma_g*p)) - gamma*p)/gamma;
      E = p/(gamma-1) + gamma*pinf/(gamma-1) + 0.5*rho*(u*u+v*v) ;
      // print k, p, pinf, E, Ef, '{0:.5e}'.format(fabs(E-Ef))  
      
      // Test for convergence and exit if converged
      get_f_pnew = gamma_g*pnew*(Et*(-1 + gamma) + (0.5*rho*(u*u+v*v)) - pnew + gamma*(-(0.5*rho*(u*u+v*v)) + pnew) - gamma_l*(pnew + pinf_l)) - gamma_l*(Et - Et*gamma + (-1 + gamma)*(0.5*rho*(u*u+v*v)) + (1 - gamma + gamma_g)*pnew)*(pnew + pinf_l)*pow(pnew/pm,gamma_i)*rv;

      if ((fabs(delta) < reltol) || (fabs(get_f_pnew)<tol)){break;}

      if (k>990){
	printf("not converging: delta = %20.16f (f=%20.16f and f'=%20.16f)\n",fabs(delta),get_f_p,get_fp_p);
      }
      
    }

    // Update the energy and pinf field
    U[(e*N_F+3)*N_s+i] = E;
    U[(e*N_F+5)*N_s+i] = gamma*pinf/(gamma-1);
    
#ifdef USE_CPU
    }
#endif
  }
}


//==========================================================================
arch_global void hack_pinf_20150813(int N_s, int N_E, scalar pm, scalar* U){
  /*!
    \brief Solve for pinf using dubious methods. See notes 2015/08/13
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] pm  initial mixture pressure (global variable from initial condition)
    \param[out] U the solution (we will modify the pinf field)
  */


#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

    // field properties
    scalar rho   = U[(e*N_F+0)*N_s+i];
    scalar u     = U[(e*N_F+1)*N_s+i]/rho;  // (rho u / rho) = u
    scalar v     = U[(e*N_F+2)*N_s+i]/rho;  // (rho v / rho) = v
    scalar Et    = U[(e*N_F+3)*N_s+i];
    scalar G     = U[(e*N_F+4)*N_s+i];
    scalar gamma = 1+1.0/G;
    scalar beta  = U[(e*N_F+5)*N_s+i];
    scalar pinf  = beta*(gamma-1)/gamma;

    // ND quantities
    // scalar rho_air = 1.1765;
    // scalar gamma_air = 1.4;
    // scalar patm = 101325;
    // scalar cs_air = sqrt(gamma_air*patm/rho_air);
    scalar p_ND   = 141855.000000; // rho_air*cs_air*cs_air;
    
    // mixture
    scalar rv = 0.8;
    scalar alpha_g = rv/(1+rv);
    scalar alpha_l = 1-alpha_g;

    // material properties (non-dimensionalized)
    scalar gamma_l = 5.5;
    scalar pinf_l = 492115000/p_ND;
    scalar gamma_g = 1.4;

    // redefine the volume fractions to account for the pressure change
    // initial pressure in flow
    //scalar pm = 1.306875;    // Ms = 2
    //scalar pm = 0.580833;    // Ms = 3
    //scalar pm = 0.326719;    // Ms = 4
    //scalar pm = 0.209100;    // Ms = 5
    //scalar pm = 0.145208;    // Ms = 6
    //scalar pm = 0.106684;    // Ms = 7
    //scalar pm = 0.081680;    // Ms = 8
    //scalar pm = 0.064537;    // Ms = 9
    //scalar pm = 0.052275;    // Ms = 10
    scalar gamma_i = gamma_g; // maybe use gamma_g?
    
    // Newton solver
    scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v)) - gamma*pinf;
    scalar E = Et;
    scalar reltol = 1e-13;
    scalar tol    = 1e-13;
    scalar get_f_p, get_f_pnew, get_fp_p, pnew,K,delta = 0;
    for(int k=0; k<1000; k++){
     
      // evaluate function and derivative
      get_f_p = -((-1 + gamma_g)*gamma_g*p*(Et - Et*gamma_l + (-1 + gamma_l)*(0.5*rho*(u*u+v*v)) + p + gamma_l*pinf_l)) + (-1 + gamma_l)*gamma_l*((-1 + gamma_g)*(Et - (0.5*rho*(u*u+v*v))) - p)*(p + pinf_l)*pow(p/pm,gamma_i)*rv;
      get_fp_p = (-((-1 + gamma_g)*gamma_g*p*(Et - Et*gamma_l + (-1 + gamma_l)*(0.5*rho*(u*u+v*v)) + 2*p + gamma_l*pinf_l)) + (-1 + gamma_l)*gamma_l*(Et*(-1 + gamma_g)*(p + gamma_i*p + gamma_i*pinf_l) - (-1 + gamma_g)*(0.5*rho*(u*u+v*v))*(p + gamma_i*p + gamma_i*pinf_l) - p*((2 + gamma_i)*p + pinf_l + gamma_i*pinf_l))*pow(p/pm,gamma_i)*rv)/p;
	
      // get the next guess
      delta = get_f_p/get_fp_p;
      pnew = p - delta;
     
      // save data first
      p = pnew;

      // get pinf and energy
      K = rv*pow(p/pm,gamma_i); 
      alpha_g = K/(1.0+K);
      alpha_l = 1-alpha_g;
      G = alpha_g/(gamma_g-1) + alpha_l/(gamma_l-1);
      gamma = 1.0+1.0/G;
      pinf = (1.0/(alpha_l/(gamma_l*(p+pinf_l)) + alpha_g/(gamma_g*p)) - gamma*p)/gamma;
      E = p/(gamma-1) + gamma*pinf/(gamma-1) + 0.5*rho*(u*u+v*v) ;
      // print k, p, pinf, E, Ef, '{0:.5e}'.format(fabs(E-Ef))  
      
      // Test for convergence and exit if converged
      get_f_p = -((-1 + gamma_g)*gamma_g*pnew*(Et - Et*gamma_l + (-1 + gamma_l)*(0.5*rho*(u*u+v*v)) + pnew + gamma_l*pinf_l)) + (-1 + gamma_l)*gamma_l*((-1 + gamma_g)*(Et - (0.5*rho*(u*u+v*v))) - pnew)*(pnew + pinf_l)*pow(pnew/pm,gamma_i)*rv;

      if ((fabs(delta) < reltol) || (fabs(get_f_pnew)<tol)){break;}

      if (k>990){
	printf("not converging: delta = %20.16f (f=%20.16f and f'=%20.16f)\n",fabs(delta),get_f_p,get_fp_p);
      }
      
    }

    // Update the energy, 1/(gamma-1) and pinf field
    U[(e*N_F+3)*N_s+i] = E;
    U[(e*N_F+4)*N_s+i] = G;
    U[(e*N_F+5)*N_s+i] = gamma*pinf/(gamma-1);
    
#ifdef USE_CPU
    }
#endif
  }
}


//==========================================================================
//
//  Host C functions
//
//==========================================================================

extern "C" 
void Lhack_pinf(int N_s, int N_E, scalar* U){
  /*!
    \brief Host C function to launch hack_pinf kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[out] U we will modify the pinf field
    \section Description
    In GPU mode, launches N_E/blkE blocks of N_s x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_s,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  //hack_pinf_20150727 arch_args (N_s, N_E, U);
  //hack_pinf_20150807 arch_args (N_s, N_E, constants::GLOBAL_P0_BBLWEDG, U);
  hack_pinf_20150813 arch_args (N_s, N_E, constants::GLOBAL_P0_BBLWEDG, U);
};
