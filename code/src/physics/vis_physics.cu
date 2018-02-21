/*!
  \file physics.cu
  \brief Kernels used for the viscous/diffusive portion of the physics
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnsen <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "vis_physics.h"
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

/*
  08/16/2017: Doing the viscous energy flux as described
  in "Overview of the NASA Glenn Flux Reconstruction Based
  High-Order Unstructured Grid Code" by Spiegel, Debonis, and Huynh.
 */

#ifdef USE_GPU
__device__ double AtomicAdd_vp(double* address, double val)
{
  /*!
    \brief Atomic addition definition for a GPU.
    \param[out] address address of value to add to
    \param[in] val value to be added to address value
  */
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif 


  //=============================================================================
arch_global void evaluate_sf_vis(int N_G, int N_E, int PROC, scalar* f, scalar* Ug, scalar* dUg){//, scalar* xyz){
    /*!
      \brief Evaluates element interior Viscous flux, adds it to interior inviscid flux
      \param[in] N_G number of gaussian nodes per element
      \param[in] N_E number of elements
      \param[in] PROC processor id.
      \param[in and out] f element flux array, incremented by viscous flux
      \param[in] Ug solution (evaluated at gaussian nodes)
      \param[in] dUg gradient of solution wrt x,y (evaluated at gaussian nodes)
    */
    //Gradients are already wrt physical coordinates, no need for correction.
    //template: dUg[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f]; a=dimensional component, g=quadrature node, f=field

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int g = threadIdx.x;
#endif


      //====================================================================
      //
      // 1D problem
      //
      //====================================================================
#ifdef ONED
    //Two choices for 1D: either scalar advection-diffusion or singlefluid
#ifdef SCALARAD
    {
      scalar rho   = Ug[e*N_F*N_G + 0*N_G + g];
      scalar rho_x = dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0];
      
      //viscosity may be zero or constant:
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif
      //add the viscous adjustment to existing flux:
      //      printf("e=%d,g=%d: rho=%f, gradrho=%f, f_input = %f,",e,g,rho,rho_x,f[((e*N_F+0)*N_G+g)*D+0]);
      //Adding something to f, so I need to use the atomic functionality in GPU case
#ifdef USE_CPU
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - mew*rho_x;
#endif
#ifdef USE_GPU
      scalar addend = (-mew)*rho_x;
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend);
#endif
      //      printf("f_output=%f\n",f[((e*N_F+0)*N_G+g)*D+0]);
	     
    } //end if for scalar ad
#endif
    

#ifdef SINGLEFLUID //=========================================================
    {
    /* Marc's array technique, quoted for reference:
       scalar rho   = Ug[(e*N_F+0)*N_G+g];   
       scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
       scalar Et    = Ug[(e*N_F+2)*N_G+g];
    */
      //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants:: GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)

      //Get the conserved variables:
      scalar rho   = Ug[e*N_F*N_G + 0*N_G + g];
      scalar mox   = Ug[e*N_F*N_G + 1*N_G + g];
      scalar Et    = Ug[e*N_F*N_G + 2*N_G + g];

      //Get some primitive variables:
      scalar u     = mox/rho;  // (rho u / rho) = u
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
      scalar rho_sq = rho*rho;
      scalar rho_cub = rho_sq*rho;
      scalar u_sq = u*u;

      //Viscosity Section:
      scalar Temp = p / (rho*Rgas); //temperature, using ideal gas law
      scalar mew;
#ifdef NOVIS
      mew = 0.0;
#elif CONSTANTVIS
      mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;
      mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //end if on viscosity model
      
      //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;
      
      //Import gradients: Subscript means derivative here
      scalar rho_x = dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0];
      scalar mox_x = dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 1];
      scalar Et_x =  dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 2];

      //u_x is veclocity gradient
      //dpr_dx is gradient of (pressure/density) term
      scalar u_x = (rho*mox_x - mox*rho_x)/rho_sq;
      scalar p_x = (gamma-1.0) * (Et_x - 0.5*(u*mox_x + mox*u_x));
      scalar dpr_dx = (rho*p_x - p*rho_x)/rho_sq; 
      scalar dT_dx = dpr_dx / Rgas; //temperature gradient, assuming ideal gas law
         
      //Calculate viscous stress/diffusion values
      //Subscript no longer means derivative
      scalar Tau_xx = 4.0/3.0 * mew * u_x;
      scalar heatflux_x = -LamT * dT_dx;
      
      // Source term not altered

      // Flux: the viscous contribution is subtracted from the previously calculated inviscid contribution
#ifdef USE_CPU
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - 0; //no change to rho flux       
      f[((e*N_F+1)*N_G+g)*D+0] = f[((e*N_F+1)*N_G+g)*D+0] - Tau_xx; //momentum flux - Tau_xx (remember sign change)
      f[((e*N_F+2)*N_G+g)*D+0] = f[((e*N_F+2)*N_G+g)*D+0] - (u*Tau_xx - heatflux_x); //heatflux_x is q_x in my code
#endif
#ifdef USE_GPU
      scalar addend_rho = 0.0;
      scalar addend_mox = -(Tau_xx);
      scalar addend_en = -(u*Tau_xx - heatflux_x);
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_rho);
      AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+0], addend_mox);
      AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+0], addend_en);
#endif
      
    }
#endif // end ifs on physics
    

      //===================================================================
      //
      // 2D problem
      //
      //===================================================================
#elif TWOD
    //Two choices for 2D: either scalar advection-diffusion or singlefluid
#ifdef SCALARAD
    {
      //sticking to laplacian diffusion for now
      scalar rho   = Ug[e*N_F*N_G + 0*N_G + g];
      scalar rho_x = dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0];
      scalar rho_y = dUg[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 0];

      //Get viscosity based on viscosity definition
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif

      //add the viscous adjustment to existing flux:
#ifdef USE_CPU
      //Flux in x direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - mew*rho_x;
      // Flux in y direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+1] = f[((e*N_F+0)*N_G+g)*D+1] - mew*rho_y;
#endif
#ifdef USE_GPU
      scalar addend_x = -mew*rho_x;
      scalar addend_y = -mew*rho_y;
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_x);
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+1], addend_y);
#endif
     
    } //end if for scalar ad
#endif
#ifdef SINGLEFLUID //========================================================
    {
      //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants:: GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)

      //Get conserved variables
      scalar rho   = Ug[(e*N_F+0)*N_G+g];
      scalar mox   = Ug[(e*N_F+1)*N_G+g];
      scalar moy   = Ug[(e*N_F+2)*N_G+g];
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      
      //Get some primitive variables:
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v));
      scalar rho_sq = rho*rho;

      //Viscosity Section:
      scalar Temp = p / (rho*Rgas);
      scalar mew;
#ifdef NOVIS
      mew = 0.0;
#elif CONSTANTVIS
      mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;
      mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //End if on viscosity model
     
       //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;
      
      

      //The input is the gradient of conserved variables wrt physical coordinates
      int fetch_x = e*D*N_G*N_F + 0*N_G*N_F + g*N_F;
      int fetch_y = e*D*N_G*N_F + 1*N_G*N_F + g*N_F;
      scalar rho_x = dUg[fetch_x + 0];
      scalar mox_x = dUg[fetch_x + 1];
      scalar moy_x = dUg[fetch_x + 2];
      scalar Et_x  = dUg[fetch_x + 3];
      scalar rho_y = dUg[fetch_y + 0];
      scalar mox_y = dUg[fetch_y + 1];
      scalar moy_y = dUg[fetch_y + 2];
      scalar Et_y  = dUg[fetch_y + 3];

      //Get velocity gradient wrt physical coordinates
      scalar u_x = (rho*mox_x - mox*rho_x)/rho_sq; 
      scalar u_y = (rho*mox_y - mox*rho_y)/rho_sq;	
      scalar v_x = (rho*moy_x - moy*rho_x)/rho_sq;	
      scalar v_y = (rho*moy_y - moy*rho_y)/rho_sq;	

      //pressure gradient:
      scalar p_x = (gamma-1)*(Et_x - 0.5*rho_x*(u*u+v*v) - 0.5*rho*(2*u*u_x + 2*v*v_x));
      scalar p_y = (gamma-1)*(Et_y - 0.5*rho_y*(u*u+v*v) - 0.5*rho*(2*u*u_y + 2*v*v_y));

      //Gradient of (pressure/rho)
      scalar dpr_dx = (rho*p_x - p*rho_x)/rho_sq; 
      scalar dpr_dy = (rho*p_y - p*rho_y)/rho_sq;

      //Get temperature gradient with ideal gas law:
      scalar dT_dx = dpr_dx / Rgas;
      scalar dT_dy = dpr_dy / Rgas;

      //Viscous stress/diffusion components
      scalar Tau_xx = 2.0*mew/3.0 * (2.0*u_x - v_y);
      scalar Tau_yy = 2.0*mew/3.0 * (2.0*v_y - u_x);
      scalar Tau_xy = mew * (u_y + v_x);
      scalar heatflux_x = - LamT * dT_dx;
      scalar heatflux_y = - LamT * dT_dy;

      // Source term: NO change to source term
      //for CPU case, I can subtract the viscous contribution from flux.
      //Need to use atomic add in the GPU case
#ifdef USE_CPU
      //Flux in x direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - 0;      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = f[((e*N_F+1)*N_G+g)*D+0] - Tau_xx; // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = f[((e*N_F+2)*N_G+g)*D+0] - Tau_xy;   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = f[((e*N_F+3)*N_G+g)*D+0] - (u*Tau_xx + v*Tau_xy - heatflux_x);  // u(E+p)
      
      // Flux in y direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+1] = f[((e*N_F+0)*N_G+g)*D+1] - 0;      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = f[((e*N_F+1)*N_G+g)*D+1] - Tau_xy;   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = f[((e*N_F+2)*N_G+g)*D+1] - Tau_yy; // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = f[((e*N_F+3)*N_G+g)*D+1] - (v*Tau_yy + u*Tau_xy - heatflux_y);  // v(E+p)
#endif
#ifdef USE_GPU
      scalar addend_x_rho = 0.0;
      scalar addend_x_mox = -Tau_xx;
      scalar addend_x_moy = -Tau_xy;
      scalar addend_x_en = -(u*Tau_xx + v*Tau_xy - heatflux_x);
      scalar addend_y_rho = 0.0;
      scalar addend_y_mox = -Tau_xy;
      scalar addend_y_moy = -Tau_yy;
      scalar addend_y_en = -(v*Tau_yy + u*Tau_xy - heatflux_y); 
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_x_rho);
      AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+0], addend_x_mox);
      AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+0], addend_x_moy);
      AtomicAdd_vp(&f[((e*N_F+3)*N_G+g)*D+0], addend_x_en);
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+1], addend_y_rho);
      AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+1], addend_y_mox);
      AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+1], addend_y_moy);
      AtomicAdd_vp(&f[((e*N_F+3)*N_G+g)*D+1], addend_y_en);
#endif

    } //end bracket for if  singleflid
#endif // end ifs on physics

#elif THREED
    //Only one choice for 3D: scalar advection-diffusion
#ifdef SCALARAD
    {
      //sticking to laplacian diffusion for now
      scalar rho   = Ug[e*N_F*N_G + 0*N_G + g];
      scalar rho_x = dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0];
      scalar rho_y = dUg[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 0];
      scalar rho_z = dUg[e*D*N_G*N_F + 2*N_G*N_F + g*N_F + 0];

      //Get the viscostity
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif
      //      printf("Evaluate_sf_vis, pr=%d: Identified rho_x=%f  in e=%d, g=%d\n", PROC, rho_x, e, g);
      //add the viscous adjustment to existing flux:
      //Some debugging stuff:
      /*
      if (fabs(rho_x)>1.1)
	{
	  printf("Evaluate_sf_vis, pr=%d: Identified rho_x=%f too big in e=%d, g=%d\n", PROC, rho_x, e, g);
	}
      if (fabs(rho_y)>0.1)
	{
	  printf("Evaluate_sf_vis, pr=%d: Identified rho_y=%f too big in e=%d, g=%d\n", PROC, rho_y, e, g);
	}
      if (fabs(rho_z)>0.1)
	{
	  printf("Evaluate_sf_vis, pr=%d: Identified rho_z=%f too big in e=%d, g=%d\n", PROC, rho_z, e, g);
	}
      */
#ifdef USE_CPU
      //Flux in x direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - mew*rho_x;
      // Flux in y direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+1] = f[((e*N_F+0)*N_G+g)*D+1] - mew*rho_y;
      // Flux in z direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+2] = f[((e*N_F+0)*N_G+g)*D+2] - mew*rho_z;
#endif
#ifdef USE_GPU
      scalar addend_x = -mew*rho_x;
      scalar addend_y = -mew*rho_y;
      scalar addend_z = -mew*rho_z;
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_x);
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+1], addend_y);
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+2], addend_z);
#endif

    }
#endif //end if for scalar ad
#ifdef SINGLEFLUID
    {
      //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants::GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)

      //Get conserved variables
      scalar rho   = Ug[(e*N_F+0)*N_G+g];
      scalar mox   = Ug[(e*N_F+1)*N_G+g];
      scalar moy   = Ug[(e*N_F+2)*N_G+g];
      scalar moz   = Ug[(e*N_F+3)*N_G+g];
      scalar Et    = Ug[(e*N_F+4)*N_G+g];
      
      //Get some primitive variables:
      scalar u     = mox / rho;  // (rho u / rho) = u
      scalar v     = moy / rho;  // (rho v / rho) = v
      scalar w     = moz / rho;  // (rho w / rho) = w
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u + v*v + w*w));
      scalar rho_sq = rho*rho;
      scalar InvRho_sq = 1.0 / rho_sq;
      scalar InvRho = 1.0 / rho;

      //Viscosity Section:
      scalar Temp = (p / Rgas) * InvRho;
      scalar mew;
#ifdef NOVIS
      mew = 0.0;
#elif CONSTANTVIS
      mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;
      mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //End if on viscosity model
      //      printf("Calling evaluate_sf_vis: mew=%f\n", mew);
       //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;

      //The input is the gradient of conserved variables wrt physical coordinates
      int fetch_x = e*D*N_G*N_F + 0*N_G*N_F + g*N_F;
      int fetch_y = e*D*N_G*N_F + 1*N_G*N_F + g*N_F;
      int fetch_z = e*D*N_G*N_F + 2*N_G*N_F + g*N_F;
      //derivatives wrt x:
      scalar rho_x = dUg[fetch_x + 0];
      scalar mox_x = dUg[fetch_x + 1];
      scalar moy_x = dUg[fetch_x + 2];
      scalar moz_x = dUg[fetch_x + 3];
      scalar Et_x  = dUg[fetch_x + 4];
      //derivatives wrt y:
      scalar rho_y = dUg[fetch_y + 0];
      scalar mox_y = dUg[fetch_y + 1];
      scalar moy_y = dUg[fetch_y + 2];
      scalar moz_y = dUg[fetch_y + 3];
      scalar Et_y  = dUg[fetch_y + 4];
      //derivatives wrt z:
      scalar rho_z = dUg[fetch_z + 0];
      scalar mox_z = dUg[fetch_z + 1];
      scalar moy_z = dUg[fetch_z + 2];
      scalar moz_z = dUg[fetch_z + 3];
      scalar Et_z  = dUg[fetch_z + 4];
      //Get velocity gradient wrt physical coordinates (uses quotient rule)
      //x velocity component (u):
      scalar u_x = (rho*mox_x - mox*rho_x) * InvRho_sq; 
      scalar u_y = (rho*mox_y - mox*rho_y) * InvRho_sq; 
      scalar u_z = (rho*mox_z - mox*rho_z) * InvRho_sq; 
      //y velocity component (v):
      scalar v_x = (rho*moy_x - moy*rho_x) * InvRho_sq; 	
      scalar v_y = (rho*moy_y - moy*rho_y) * InvRho_sq; 
      scalar v_z = (rho*moy_z - moy*rho_z) * InvRho_sq; 
      //z velocity component (w):
      scalar w_x = (rho*moz_x - moz*rho_x) * InvRho_sq; 	
      scalar w_y = (rho*moz_y - moz*rho_y) * InvRho_sq; 
      scalar w_z = (rho*moz_z - moz*rho_z) * InvRho_sq; 

      //Gradient of (pressure/rho): see notebook 09/25/2017
      scalar Lead_Eg = (gamma-1.0) * InvRho_sq;
      scalar Lead_V = (gamma-1.0);
      scalar dpr_dx = Lead_Eg*(rho*Et_x - Et*rho_x) - Lead_V*(u*u_x + v*v_x + w*w_x);
      scalar dpr_dy = Lead_Eg*(rho*Et_y - Et*rho_y) - Lead_V*(u*u_y + v*v_y + w*w_y);
      scalar dpr_dz = Lead_Eg*(rho*Et_z - Et*rho_z) - Lead_V*(u*u_z + v*v_z + w*w_z);

      //Viscous stress/diffusion components
      //Tnm stands for "tau-no-mew"
      scalar DivV = u_x + v_y + w_z; //divergence of velocity field
      scalar TwoThr = 2.0/3.0;
      scalar Tnm_XX = u_x + u_x - TwoThr * DivV;
      scalar Tnm_XY = u_y + v_x;
      scalar Tnm_XZ = u_z + w_x;

      scalar Tnm_YX = Tnm_XY;
      scalar Tnm_YY = v_y + v_y - TwoThr * DivV;
      scalar Tnm_YZ = v_z + w_y;

      scalar Tnm_ZX = Tnm_XZ;
      scalar Tnm_ZY = Tnm_YZ;
      scalar Tnm_ZZ = w_z + w_z - TwoThr * DivV;
     
      //gradT = grad(p/rho) / Rg;
      scalar Mult = -LamT/Rgas; 
      scalar heatflux_x = Mult * dpr_dx;
      scalar heatflux_y = Mult * dpr_dy;
      scalar heatflux_z = Mult * dpr_dz;

      //We are now ready to adjust the fluxes
#ifdef USE_CPU
      //Flux in x direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - 0;           
      f[((e*N_F+1)*N_G+g)*D+0] = f[((e*N_F+1)*N_G+g)*D+0] - mew * Tnm_XX; 
      f[((e*N_F+2)*N_G+g)*D+0] = f[((e*N_F+2)*N_G+g)*D+0] - mew * Tnm_YX; 
      f[((e*N_F+3)*N_G+g)*D+0] = f[((e*N_F+3)*N_G+g)*D+0] - mew * Tnm_ZX; 
      f[((e*N_F+4)*N_G+g)*D+0] = f[((e*N_F+4)*N_G+g)*D+0] - mew*(u*Tnm_XX + v*Tnm_XY + w*Tnm_XZ) + heatflux_x; 
      
      // Flux in y direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+1] = f[((e*N_F+0)*N_G+g)*D+1] - 0;        
      f[((e*N_F+1)*N_G+g)*D+1] = f[((e*N_F+1)*N_G+g)*D+1] - mew * Tnm_XY; 
      f[((e*N_F+2)*N_G+g)*D+1] = f[((e*N_F+2)*N_G+g)*D+1] - mew * Tnm_YY; 
      f[((e*N_F+3)*N_G+g)*D+1] = f[((e*N_F+3)*N_G+g)*D+1] - mew * Tnm_ZY; 
      f[((e*N_F+4)*N_G+g)*D+1] = f[((e*N_F+4)*N_G+g)*D+1] - mew*(u*Tnm_YX + v*Tnm_YY + w*Tnm_YZ) + heatflux_y; 

      // Flux in z direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+2] = f[((e*N_F+0)*N_G+g)*D+2] - 0;        
      f[((e*N_F+1)*N_G+g)*D+2] = f[((e*N_F+1)*N_G+g)*D+2] - mew * Tnm_XZ; 
      f[((e*N_F+2)*N_G+g)*D+2] = f[((e*N_F+2)*N_G+g)*D+2] - mew * Tnm_YZ; 
      f[((e*N_F+3)*N_G+g)*D+2] = f[((e*N_F+3)*N_G+g)*D+2] - mew * Tnm_ZZ; 
      f[((e*N_F+4)*N_G+g)*D+2] = f[((e*N_F+4)*N_G+g)*D+2] - mew*(u*Tnm_ZX + v*Tnm_ZY + w*Tnm_ZZ) + heatflux_z; 
#endif
#ifdef USE_GPU
      printf("Fuck off, evaluate_sf_vis not ready for GPU implementation\n");
#endif

    }
#endif //end if for singlefluid physics
#endif // end ifs on dimensions  
#ifdef USE_CPU
  }
#endif
    }
  }


//==========================================================================
  arch_global void evaluate_q_vis(int M_G, int M_T, int PROC, scalar* q, scalar* UhCommon, scalar* gradCommon, scalar* normals){//, scalar* xyzf){
    /*!
      \brief Evaluate interface viscous fluxes
      \param[in] M_G number of gaussian nodes per interface
      \param[in] M_T number of interfaces
      \param[in] PROC processor id.
      \param[in and out] q interface flux array, incremented by viscous flux
      \param[in] UhCommon interface solution (evaluated at gaussian nodes)
      \param[in] gradCommon is the gradient along the interface wrt physical coordinates
      \param[in] normals interface normals
    */
    //Viscous interface flux subroutine. Adds a viscous contribution to
    //the inviscid q vector
    
    // Apply the fluxes by subtracting from inviscid Riemann solution.
    //The normal points out of cell 0 (I hope), so sign is reversed for cell 1 contribution
    
    //My interface indexing is different than Marc's: I store 1 common U and 1 common gradU for each interface quadrature point,
    //where Marc stores states of left/right variables. For fluxes, must work with Marc's indexing

    //template: UhCommon[t*M_G*N_F + g*N_F + f]
    //template: gradCommon[t*M_G*N_F*D + g*N_F*D + f*D + a]
    
#ifdef USE_CPU
    int blk = 0;
    for(int t = 0; t < M_T; t++){
      for(int g = 0; g < M_G; g++){
#elif USE_GPU
  int blk = threadIdx.z; // buffer needs room for all the elements in the block
  int t = blockIdx.x*blkT+blk;
  if ( t < M_T){
    int g = threadIdx.x;
#endif

#ifdef ONED
#ifdef SCALARAD
    {
      //Conserved Variables:
    scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
      //viscosity may be zero or constant:
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif

     //Grad of conserved variables
    scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
    scalar nx = -normals[t*D+0];
    //The fluxes on each side of the interface must be averaged
    scalar Qvis_rho = nx*mew*rho_x;

    //add viscous term to interface flux:
#ifdef USE_CPU
    //    printf("t=%d, g=%d, rho=%f, gradrho=%f, q_in(L/R)=(%f,%f), ",t,g, rho,rho_x,q[((t*N_F+0)*2+0)*M_G+g],q[((t*N_F+0)*2+1)*M_G+g]);
    q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
    q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;
    //    printf("q_out(L/R)=(%f,%f)\n", q[((t*N_F+0)*2+0)*M_G+g],q[((t*N_F+0)*2+1)*M_G+g]);
#endif
#ifdef USE_GPU
    scalar addend_A = -Qvis_rho;
    scalar addend_B = Qvis_rho;
    AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A);
    AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B);
#endif
    }
#endif
#ifdef  SINGLEFLUID //=========================================================
    {
    //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants:: GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)
   
    //Get the onserved Variables:
    scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
    scalar mox = UhCommon[t*M_G*N_F + g*N_F + 1];
    scalar Et =  UhCommon[t*M_G*N_F + g*N_F + 2];

    //Get some primitive variables:
    scalar u     = mox/rho;  // (rho u / rho) = u
    scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
    scalar rho_sq = rho*rho;
    scalar rho_cub = rho_sq*rho;
    scalar u_sq = u*u;
    
    //Viscosity Section:
    scalar Temp = p / (rho*Rgas); //temperature, using ideal gas law
    scalar mew;
#ifdef NOVIS
    mew = 0.0;
#elif CONSTANTVIS
      mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;
      mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //end if on viscosity model

      //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;
      
      //Import gradients: Subscript means derivative here
      scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
      scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
      scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
      
      //u_x is veclocity gradient
      //dpr_dx is gradient of (pressure/density) term
      scalar u_x = (rho*mox_x - mox*rho_x)/rho_sq;
      scalar p_x = (gamma-1.0) * (Et_x - 0.5*(u*mox_x + mox*u_x));
      scalar dpr_dx = (rho*p_x - p*rho_x)/rho_sq; 
      scalar dT_dx = dpr_dx / Rgas; //temperature gradient, assuming ideal gas law
      
      //Calculate viscous stress/diffusion values
      //Subscript no longer means derivative
      scalar Tau_xx = 4.0/3.0 * mew * u_x;
      scalar heatflux_x = -LamT * dT_dx;
      
      //interface normal:
      //scalar nx = -normals[t*D+0];
      scalar nx = -normals[t*D+0];

      //viscous flux contributions:
      scalar Qvis_rho = 0;
      scalar Qvis_mox = nx*Tau_xx;
      scalar Qvis_en = nx*(u*Tau_xx - heatflux_x);
    
    //Adjust the Q vector based on your findings:
    //The sign is based on the "Q dot n" term in surface residual
#ifdef USE_CPU
    q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
    q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;

    q[((t*N_F+1)*2+0)*M_G+g] = q[((t*N_F+1)*2+0)*M_G+g] - Qvis_mox;
    q[((t*N_F+1)*2+1)*M_G+g] = q[((t*N_F+1)*2+1)*M_G+g] + Qvis_mox;
    
    q[((t*N_F+2)*2+0)*M_G+g] = q[((t*N_F+2)*2+0)*M_G+g] - Qvis_en;
    q[((t*N_F+2)*2+1)*M_G+g] = q[((t*N_F+2)*2+1)*M_G+g] + Qvis_en;
#endif
#ifdef USE_GPU
    scalar addend_A_rho = -Qvis_rho;
    scalar addend_B_rho = Qvis_rho;
    scalar addend_A_mox = -Qvis_mox;
    scalar addend_B_mox = Qvis_mox;
    scalar addend_A_en = -Qvis_en;
    scalar addend_B_en = Qvis_en;
    
    AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A_rho);
    AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B_rho);
    AtomicAdd_vp(&q[((t*N_F+1)*2+0)*M_G+g], addend_A_mox);
    AtomicAdd_vp(&q[((t*N_F+1)*2+1)*M_G+g], addend_B_mox);
    AtomicAdd_vp(&q[((t*N_F+2)*2+0)*M_G+g], addend_A_en);
    AtomicAdd_vp(&q[((t*N_F+2)*2+1)*M_G+g], addend_B_en);
#endif
    }
#endif // physics if
    
    
#elif TWOD      /////////////////////////////////////////////////////////
#ifdef SCALARAD
    {
      //Conserved Variables:
    scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
      //Get viscosity based on viscosity definition
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif
    //Grad of conserved variables
    scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
    scalar rho_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 1];
    scalar nx = -normals[t*D+0];
    scalar ny = -normals[t*D+1];
 
    scalar Qvis_rho = nx*mew*rho_x + ny*mew*rho_y;
    

    //add viscous term to interface flux:
    //Since this function includes the dot product between
    //local flux and local normal, it only has one component per field variable
#ifdef USE_CPU
    q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
    q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;
#endif
#ifdef USE_GPU
    scalar addend_A = -Qvis_rho;
    scalar addend_B = Qvis_rho;
    AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A);
    AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B);
#endif

    }
#endif
#ifdef  SINGLEFLUID //=========================================================
    {
      //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants:: GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)
      
     //Get conserved variables
      scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
      scalar mox = UhCommon[t*M_G*N_F + g*N_F + 1];
      scalar moy = UhCommon[t*M_G*N_F + g*N_F + 2];
      scalar Et = UhCommon[t*M_G*N_F + g*N_F + 3];
      
      //Get some primitive variables:
      scalar u     = mox/rho;  // (rho u / rho) = u
      scalar v     = moy/rho;  // (rho v / rho) = v
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v));
      scalar rho_sq = rho*rho;
    
      //Viscosity Section:
      scalar Temp = p / (rho*Rgas);
      scalar mew;
#ifdef NOVIS
      mew = 0.0;
#elif CONSTANTVIS
    mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
    scalar mew_ref = constants:: GLOBAL_MEWREF;
    scalar T_ref = constants::GLOBAL_TREF;
    scalar Cvis = constants::GLOBAL_CVIS;
    mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //End if on viscosity model
    
     //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;



    //The input is the gradient of conserved variables wrt physical coordinates
    scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
    scalar rho_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 1];
    
    scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
    scalar mox_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 1];
    
    scalar moy_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
    scalar moy_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 1];
    
    scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 0];
    scalar Et_y  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 1];
    

     //Get velocity gradient wrt physical coordinates
      scalar u_x = (rho*mox_x - mox*rho_x)/rho_sq; 
      scalar u_y = (rho*mox_y - mox*rho_y)/rho_sq;	
      scalar v_x = (rho*moy_x - moy*rho_x)/rho_sq;	
      scalar v_y = (rho*moy_y - moy*rho_y)/rho_sq;	

      //pressure gradient:
      scalar p_x = (gamma-1)*(Et_x - 0.5*rho_x*(u*u+v*v) - 0.5*rho*(2*u*u_x + 2*v*v_x));
      scalar p_y = (gamma-1)*(Et_y - 0.5*rho_y*(u*u+v*v) - 0.5*rho*(2*u*u_y + 2*v*v_y));

      //Gradient of (pressure/rho)
      scalar dpr_dx = (rho*p_x - p*rho_x)/rho_sq; 
      scalar dpr_dy = (rho*p_y - p*rho_y)/rho_sq;

      //Get temperature gradient with ideal gas law:
      scalar dT_dx = dpr_dx / Rgas;
      scalar dT_dy = dpr_dy / Rgas;

      //Viscous stress/diffusion components
      scalar Tau_xx = 2.0*mew/3.0 * (2.0*u_x - v_y);
      scalar Tau_yy = 2.0*mew/3.0 * (2.0*v_y - u_x);
      scalar Tau_xy = mew * (u_y + v_x);
      scalar heatflux_x = - LamT * dT_dx;
      scalar heatflux_y = - LamT * dT_dy;




	//Viscous Flux Contribution, averaged between the two competing cells:
	scalar nx = -normals[t*D+0];
	scalar ny = -normals[t*D+1];

	scalar Qvis_rho, Qvis_mox, Qvis_moy, Qvis_en;

	Qvis_rho = 0;
	Qvis_mox = (nx)*(Tau_xx) + (ny)*(Tau_xy);
	Qvis_moy = (nx)*(Tau_xy) + (ny)*(Tau_yy);
	Qvis_en = (nx)*(u*Tau_xx + v*Tau_xy - heatflux_x) + (ny)*(v*Tau_yy + u*Tau_xy - heatflux_y);
	
	// Source term: NO change to source term
	//#endif // physics if
	//The method of adjusting q depends on how many field variables there are, and also on direction
#ifdef USE_CPU
	q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
	q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;

	q[((t*N_F+1)*2+0)*M_G+g] = q[((t*N_F+1)*2+0)*M_G+g] - Qvis_mox;
	q[((t*N_F+1)*2+1)*M_G+g] = q[((t*N_F+1)*2+1)*M_G+g] + Qvis_mox;

	q[((t*N_F+2)*2+0)*M_G+g] = q[((t*N_F+2)*2+0)*M_G+g] - Qvis_moy;
	q[((t*N_F+2)*2+1)*M_G+g] = q[((t*N_F+2)*2+1)*M_G+g] + Qvis_moy;

	q[((t*N_F+3)*2+0)*M_G+g] = q[((t*N_F+3)*2+0)*M_G+g] - Qvis_en;
	q[((t*N_F+3)*2+1)*M_G+g] = q[((t*N_F+3)*2+1)*M_G+g] + Qvis_en;
#endif
#ifdef USE_GPU
	scalar addend_A_rho = -Qvis_rho;
	scalar addend_B_rho = Qvis_rho;

	scalar addend_A_mox = -Qvis_mox;
	scalar addend_B_mox = Qvis_mox;

	scalar addend_A_moy = -Qvis_moy;
	scalar addend_B_moy = Qvis_moy;

	scalar addend_A_en = -Qvis_en;
	scalar addend_B_en = Qvis_en;
	
	AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A_rho);
	AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B_rho);

	AtomicAdd_vp(&q[((t*N_F+1)*2+0)*M_G+g], addend_A_mox);
	AtomicAdd_vp(&q[((t*N_F+1)*2+1)*M_G+g], addend_B_mox);

	AtomicAdd_vp(&q[((t*N_F+2)*2+0)*M_G+g], addend_A_moy);
	AtomicAdd_vp(&q[((t*N_F+2)*2+1)*M_G+g], addend_B_moy);

	AtomicAdd_vp(&q[((t*N_F+3)*2+0)*M_G+g], addend_A_en);
	AtomicAdd_vp(&q[((t*N_F+3)*2+1)*M_G+g], addend_B_en);
#endif
    }
#endif

#elif THREED
#ifdef SCALARAD
    {
    //Conserved Variables:
    scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
      //Get viscosity based on viscosity definition
#ifdef CONSTANTVIS
      scalar mew = constants::GLOBAL_KLIN;
#endif
#ifdef NOVIS
      scalar mew = 0.0;
#endif
    //Grad of conserved variables
    scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
    scalar rho_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 1];
    scalar rho_z = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 2];
    scalar nx = -normals[t*D+0];
    scalar ny = -normals[t*D+1];
    scalar nz = -normals[t*D+2];
    
    //Some debugging stuff:
    /*
    if (fabs(rho_x)>1.5)
      {
	printf("Evaluate_q_vis, pr=%d: Identified rho_x=%f too big in t=%d, g=%d\n", PROC, rho_x, t, g);
      }
    if (fabs(rho_y)>0.1)
      {
	printf("Evaluate_q_vis, pr=%d: Identified rho_y=%f too big in t=%d, g=%d\n", PROC, rho_y, t, g);
      }
    if (fabs(rho_z)>0.1)
      {
	printf("Evaluate_q_vis, pr=%d: Identified rho_z=%f too big in t=%d, g=%d\n", PROC, rho_z, t, g);
      }
    */
    scalar Qvis_rho = nx*mew*rho_x + ny*mew*rho_y + nz*mew*rho_z;
    
    //add viscous term to interface flux:
    //Since this function includes the dot product between
    //local flux and local normal, it only has one component per field variable
#ifdef USE_CPU
    q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
    q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;
#endif
#ifdef USE_GPU
    scalar addend_A = -Qvis_rho;
    scalar addend_B = Qvis_rho;
    AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A);
    AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B);
#endif
    }
#endif //end if for scalarAD physics
#ifdef SINGLEFLUID
    {
      //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
      scalar Rgas =   constants::GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
      scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)
      
      //Get conserved variables
      scalar rho = UhCommon[t*M_G*N_F + g*N_F + 0];
      scalar mox = UhCommon[t*M_G*N_F + g*N_F + 1];
      scalar moy = UhCommon[t*M_G*N_F + g*N_F + 2];
      scalar moz = UhCommon[t*M_G*N_F + g*N_F + 3];
      scalar Et  = UhCommon[t*M_G*N_F + g*N_F + 4];
      
      //Get some primitive variables:
      scalar u     = mox / rho;  // (rho u / rho) = u
      scalar v     = moy / rho;  // (rho v / rho) = v
      scalar w     = moz / rho;  // (rho w / rho) = w
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u + v*v + w*w));
      scalar rho_sq = rho*rho;
      scalar InvRho_sq = 1.0 / rho_sq;
      scalar InvRho = 1.0 / rho;
      
      //Viscosity Section:
      scalar Temp = (p / Rgas) * InvRho;
      scalar mew;
#ifdef NOVIS
      mew = 0.0;
#elif CONSTANTVIS
      mew = constants:: GLOBAL_MEWREF;
#elif SUTHERLAND
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;
      mew = mew_ref * (T_ref+Cvis) / (Temp+Cvis) * pow(Temp/T_ref,1.5);
#endif //End if on viscosity model
      //printf("Calling evaluate_q_vis, mew=%f\n", mew);
       //Thermal diffusivity:
      scalar LamT = mew * CpGas / Pran;

      //The input is the gradient of conserved variables wrt physical coordinates
      scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
      scalar rho_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 1];
      scalar rho_z = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 2];
      
      scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
      scalar mox_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 1];
      scalar mox_z = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 2];
      
      scalar moy_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
      scalar moy_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 1];
      scalar moy_z = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 2];

      scalar moz_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 0];
      scalar moz_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 1];
      scalar moz_z = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 2];
      
      scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 0];
      scalar Et_y  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 1];
      scalar Et_z  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 2];

      //Get velocity gradient wrt physical coordinates (uses quotient rule)
      //x velocity component (u):
      scalar u_x = (rho*mox_x - mox*rho_x) * InvRho_sq; 
      scalar u_y = (rho*mox_y - mox*rho_y) * InvRho_sq; 
      scalar u_z = (rho*mox_z - mox*rho_z) * InvRho_sq; 
      //y velocity component (v):
      scalar v_x = (rho*moy_x - moy*rho_x) * InvRho_sq; 	
      scalar v_y = (rho*moy_y - moy*rho_y) * InvRho_sq; 
      scalar v_z = (rho*moy_z - moy*rho_z) * InvRho_sq; 
      //z velocity component (w):
      scalar w_x = (rho*moz_x - moz*rho_x) * InvRho_sq; 	
      scalar w_y = (rho*moz_y - moz*rho_y) * InvRho_sq; 
      scalar w_z = (rho*moz_z - moz*rho_z) * InvRho_sq; 

      //Gradient of (pressure/rho): see notebook 09/25/2017
      scalar Lead_Eg = (gamma-1.0) * InvRho_sq;
      scalar Lead_V = (gamma-1.0);
      scalar dpr_dx = Lead_Eg*(rho*Et_x - Et*rho_x) - Lead_V*(u*u_x + v*v_x + w*w_x);
      scalar dpr_dy = Lead_Eg*(rho*Et_y - Et*rho_y) - Lead_V*(u*u_y + v*v_y + w*w_y);
      scalar dpr_dz = Lead_Eg*(rho*Et_z - Et*rho_z) - Lead_V*(u*u_z + v*v_z + w*w_z);

      //Viscous stress/diffusion components
      //Tnm stands for "tau-no-mew"
      scalar DivV = u_x + v_y + w_z; //divergence of velocity field
      scalar TwoThr = 2.0/3.0;
      scalar Tnm_XX = u_x + u_x - TwoThr * DivV;
      scalar Tnm_XY = u_y + v_x;
      scalar Tnm_XZ = u_z + w_x;

      scalar Tnm_YX = Tnm_XY;
      scalar Tnm_YY = v_y + v_y - TwoThr * DivV;
      scalar Tnm_YZ = v_z + w_y;

      scalar Tnm_ZX = Tnm_XZ;
      scalar Tnm_ZY = Tnm_YZ;
      scalar Tnm_ZZ = w_z + w_z - TwoThr * DivV;
     
      //gradT = grad(p/rho) / Rg;
      scalar Mult = -LamT/Rgas; 
      scalar heatflux_x = Mult * dpr_dx;
      scalar heatflux_y = Mult * dpr_dy;
      scalar heatflux_z = Mult * dpr_dz;

      //Viscous Flux Contribution, averaged between the two competing cells:
      scalar nx = -normals[t*D+0];
      scalar ny = -normals[t*D+1];
      scalar nz = -normals[t*D+2];

      scalar Qvis_rho, Qvis_mox, Qvis_moy, Qvis_moz, Qvis_en;

      Qvis_rho = 0;
      Qvis_mox = mew*(nx*Tnm_XX + ny*Tnm_XY + nz*Tnm_XZ);
      Qvis_moy = mew*(nx*Tnm_YX + ny*Tnm_YY + nz*Tnm_YZ);
      Qvis_moz = mew*(nx*Tnm_ZX + ny*Tnm_ZY + nz*Tnm_ZZ);
      //increment energy component piece-by-pice
      Qvis_en =  nx*(mew*(u*Tnm_XX + v*Tnm_XY + w*Tnm_XZ) - heatflux_x);
      Qvis_en += ny*(mew*(u*Tnm_YX + v*Tnm_YY + w*Tnm_YZ) - heatflux_y);
      Qvis_en += nz*(mew*(u*Tnm_ZX + v*Tnm_ZY + w*Tnm_ZZ) - heatflux_z);
            
#ifdef USE_CPU
	q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qvis_rho;
	q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qvis_rho;

	q[((t*N_F+1)*2+0)*M_G+g] = q[((t*N_F+1)*2+0)*M_G+g] - Qvis_mox;
	q[((t*N_F+1)*2+1)*M_G+g] = q[((t*N_F+1)*2+1)*M_G+g] + Qvis_mox;

	q[((t*N_F+2)*2+0)*M_G+g] = q[((t*N_F+2)*2+0)*M_G+g] - Qvis_moy;
	q[((t*N_F+2)*2+1)*M_G+g] = q[((t*N_F+2)*2+1)*M_G+g] + Qvis_moy;

	q[((t*N_F+3)*2+0)*M_G+g] = q[((t*N_F+3)*2+0)*M_G+g] - Qvis_moz;
	q[((t*N_F+3)*2+1)*M_G+g] = q[((t*N_F+3)*2+1)*M_G+g] + Qvis_moz;

	q[((t*N_F+4)*2+0)*M_G+g] = q[((t*N_F+4)*2+0)*M_G+g] - Qvis_en;
	q[((t*N_F+4)*2+1)*M_G+g] = q[((t*N_F+4)*2+1)*M_G+g] + Qvis_en;
#endif
#ifdef USE_GPU
        printf("FUCK OFF, evaluate_q_vis not ready for gpu implementation\n");
#endif
    }
#endif //end if for singlefluid physics
#endif // dimension if
      
	
	

	/*for(int fc = 0; fc < N_F; fc++){
	q[((t*N_F+fc)*2+0)*M_G+g] = q[((t*N_F+fc)*2+0)*M_G+g] - Qvis[fc];
	q[((t*N_F+fc)*2+1)*M_G+g] = q[((t*N_F+fc)*2+1)*M_G+g] + Qvis[fc];
	//q[((t*N_F+fc)*2+0)*M_G+g] =-buffer[Fidx+fc] + buffer[ncidx+fc];
	//q[((t*N_F+fc)*2+1)*M_G+g] = buffer[Fidx+fc] + buffer[ncidx+fc];
      
	}*/
      
#ifdef USE_CPU
    } // end loop on g
#endif
  } // end loop on t
}




//==========================================================================
//
//  Host C functions
//
//==========================================================================


extern "C" 
  void Levaluate_sf_vis(int N_G, int N_E, int PROC, scalar* f, scalar* Ug, scalar* dUg){
  /*!
    \brief Host C function to lauch evaluate_sf_vis kernel.
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of element
    \param[out] f element flux array
    \param[in] Ug solution (evaluated at gaussian nodes)
    \param[in] dUg derivative of solution (evaluated at gaussian nodes)

    \section Description
    In GPU mode, launches N_E/blkE blocks of N_G x 1 x blkE
    threads. blkE controls the number of elements to set on each block
  */
  
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  evaluate_sf_vis arch_args (N_G, N_E, PROC, f, Ug, dUg);
}

extern "C" 
  void Levaluate_q_vis(int M_G, int M_T, int PROC, scalar* q, scalar* UhCommon, scalar* gradCommon, scalar* normals){//, scalar* xyzf){
  /*!
    \brief Host C function to lauch evaluate_q kernel.
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[out] q interface flux array
    \param[in] UgF interface solution (evaluated at gaussian nodes)
    \param[in] normals interface normals
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_G x 1 x blkT
    threads. blkT controls the number of elements to set on each block
  */

#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,1,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  evaluate_q_vis arch_args (M_G, M_T, PROC, q, UhCommon, gradCommon, normals);
}

