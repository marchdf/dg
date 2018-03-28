/*!
  \file physics.cu
  \brief Kernels used for the physics
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "physics.h"
#include "basic_fluxes.h"
#include "oned_scalarad_fluxes.h" //PEJ 05/29/2017
#include "twod_scalarad_fluxes.h" //PEJ 05/29/2017
#include "threed_scalarad_fluxes.h" //PEJ 05/29/2017
#include "threed_singlefluid_fluxes.h" //PEJ 09/26/2017
#include "oned_radsinglefluid_fluxes.h" //PEJ 10/13/2017
#include "twod_radsinglefluid_fluxes.h" //PEJ 10/16/2017
#include "oned_passive_fluxes.h"
#include "twod_passive_fluxes.h"
#include "oned_singlefluid_fluxes.h"
#include "twod_singlefluid_fluxes.h"
#include "oned_multifluid_fluxes.h"
#include "twod_multifluid_fluxes.h"
#include "oned_stiffened_fluxes.h"
#include "twod_stiffened_fluxes.h"
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void evaluate_sf(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar GX, scalar GY){//, scalar* xyz){
  /*!
    \brief Evaluate source and element fluxes
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[out] s source array
    \param[out] f element flux array
    \param[in] Ug solution (evaluated at gaussian nodes)
    \param[in] dUg derivative of solution (evaluated at gaussian nodes)
    \param[in] invJac inverse Jacobian of the elements
    \param[in] GX gravity in x-direction
    \param[in] GY gravity in y-direction
  */

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

#ifdef SCALARAD
    {
      //PEJ 05/29/2017: scalar advection-diffusion
      scalar rho = Ug[(e*N_F + 0)*N_G+g];
      //flux has one component per spatial dimension, so just one in this case
      //Source term
      s[(e*N_F+0)*N_G+g] = 0.0;

      //flux in x direction
      scalar Vx = constants::GLOBAL_VX;
      f[((e*N_F+0)*N_G+g)*D + 0] = flux_ab(rho, Vx); //this is Vx times rho.
    }
#endif //end if for scalarAD. Hopefully, nothing else in the ONED section activates now

#ifdef PASSIVE //===========================================================
    {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);

     
      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;
      s[(e*N_F+4)*N_G+g] = -u*dUg[(e*N_F+4)*N_G+g]*invJac[e*N_G+g];// = -u*dphincdx;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(Et+p,u);
      f[((e*N_F+3)*N_G+g)*D+0] = flux_abc(rho,u,Ug[(e*N_F+3)*N_G+g]/rho); // rho*u*phic
      f[((e*N_F+4)*N_G+g)*D+0] = 0;
    }
#elif SINGLEFLUID //=========================================================
    {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
     
      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(Et+p,u);
    }
#elif RADSINGLEFLUID //=====================================================
    {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
      scalar C     = Ug[(e*N_F+3)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
     
      // Source term (for this source term, I deal with it in evaluate_q_rad routine)
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(Et+p,u);
      f[((e*N_F+3)*N_G+g)*D+0] = 0.0; //C flux dealt with in evaluate_q_rad.
    }
#elif MULTIFLUID //=========================================================
    {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
#ifdef GAMCONS
      scalar gamma=1+rho/Ug[(e*N_F+3)*N_G+g];
#elif  GAMNCON
      scalar gamma=1+1.0/Ug[(e*N_F+3)*N_G+g];
#endif
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
      // Adjust pressure for gravity (flux method)
      // scalar x = xyz[(e*N_G+g)*D+0];
      //p = p -rho*constants::GLOBAL_GX*x; //p = 1e5;
      
      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = rho*GX; 
      s[(e*N_F+2)*N_G+g] = rho*u*GX;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(Et+p,u);

      //printf("p=%f,rho=%f,g=%f,x=%f,p0 = %20.16e\n",p,rho,constants::GLOBAL_GX,x,p-rho*constants::GLOBAL_GX*x);
      
#ifdef GAMCONS
      s[(e*N_F+3)*N_G+g] = 0;
      f[((e*N_F+3)*N_G+g)*D+0] = flux_abc(rho,u,1/(gamma-1));
#elif  GAMNCON
      s[(e*N_F+3)*N_G+g] = -u*dUg[(e*N_F+3)*N_G+g]*invJac[e*N_G+g]; // = -u*dalphadx
      f[((e*N_F+3)*N_G+g)*D+0] = 0;
#endif

      // Mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) s[(e*N_F+4+x)*N_G+g] = 0; \
    f[((e*N_F+4+x)*N_G+g)*D+0] = flux_ab(Ug[(e*N_F+4+x)*N_G+g],u); // Ug contains rho
#include "loop.h"

//       // Mass fractions N-C
// #include "loopstart.h"
// #define LOOP_END N_Y
// #define MACRO(x) s[(e*N_F+4+x)*N_G+g] = -u*dUg[(e*N_F+4+x)*N_G+g]*invJac[e*N_G+g];\
//     f[((e*N_F+4+x)*N_G+g)*D+0] = 0; 
// #include "loop.h"
    }
#elif STIFFENED //=========================================================
    {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
      scalar gamma = 1+1.0/Ug[(e*N_F+3)*N_G+g];
      scalar beta  = Ug[(e*N_F+4)*N_G+g];
      scalar p = (gamma-1)*(Et - beta - 0.5*rho*u*u);

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = rho*GX;
      s[(e*N_F+2)*N_G+g] = rho*u*GX;
      s[(e*N_F+3)*N_G+g] = -u*dUg[(e*N_F+3)*N_G+g]*invJac[e*N_G+g]; // = -u*dalphadx
      s[(e*N_F+4)*N_G+g] = -u*dUg[(e*N_F+4)*N_G+g]*invJac[e*N_G+g]; // = -u*dbetadx
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(Et+p,u);
      f[((e*N_F+3)*N_G+g)*D+0] = 0;
      f[((e*N_F+4)*N_G+g)*D+0] = 0;

      // Mass fractions
      int fcnt = 5;
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) s[(e*N_F+fcnt+x)*N_G+g] = 0; \
    f[((e*N_F+fcnt+x)*N_G+g)*D+0] = flux_ab(Ug[(e*N_F+fcnt+x)*N_G+g],u); fcnt++; // Ug contains rho
#include "loop.h"
    }
#endif // end ifs on physics


      //===================================================================
      //
      // 2D problem
      //
      //===================================================================
#elif TWOD

#ifdef SCALARAD
      {
	//PEJ 05/29/2017: scalar advection-diffusion
	scalar rho = Ug[(e*N_F + 0)*N_G+g];
	//flux has one component per spatial dimension, so two in this case
	
	//Source term is zero
	s[(e*N_F+0)*N_G+g] = 0.0;
	
	//flux in x direction
	scalar Vx = constants::GLOBAL_VX;
	f[((e*N_F+0)*N_G+g)*D + 0] = flux_ab(rho, Vx); //this is Vx times rho.
	//flux in y direction
	scalar Vy = constants::GLOBAL_VY;
	f[((e*N_F+0)*N_G+g)*D + 1] = flux_ab(rho, Vy); //this is Vy times rho.
      }
#endif //end if for scalarAD. Hopefully, nothing else in the TWOD section activates now
#ifdef PASSIVE //==========================================================
      {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar phic  = Ug[(e*N_F+4)*N_G+g]/rho; 
      scalar vdotgradphinc = 0; // = -u*dphincdx-v*dphincdy
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradphinc += 
	  -u*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha]; // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
      }
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v));

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;
      s[(e*N_F+4)*N_G+g] = 0;
      s[(e*N_F+5)*N_G+g] = vdotgradphinc;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(Et+p,u);  // u(E+p)
      f[((e*N_F+4)*N_G+g)*D+0] = flux_abc(rho,u,phic);// rho*u*phic
      f[((e*N_F+5)*N_G+g)*D+0] = 0;
      
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(Et+p,v);  // v(E+p)
      f[((e*N_F+4)*N_G+g)*D+1] = flux_abc(rho,v,phic);// rho*v*phic
      f[((e*N_F+5)*N_G+g)*D+1] = 0; 
      }
#elif SINGLEFLUID //========================================================
      {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v));

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(Et+p,u);  // u(E+p)
      
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(Et+p,v);  // v(E+p)
      }
#elif RADSINGLEFLUID
     {
       //Just like singlefluid except that we add in a few zeros for the 
       //C variable
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u + v*v));

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;
      s[(e*N_F+4)*N_G+g] = 0;
      s[(e*N_F+5)*N_G+g] = 0;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(Et+p,u);  // u(E+p)
      f[((e*N_F+4)*N_G+g)*D+0] = 0.0; //no inviscid flux for C
      f[((e*N_F+5)*N_G+g)*D+0] = 0.0; //no inviscid flux for C
            
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(Et+p,v);  // v(E+p)
      f[((e*N_F+4)*N_G+g)*D+1] = 0.0; //no inviscid flux for C
      f[((e*N_F+5)*N_G+g)*D+1] = 0.0; //no inviscid flux for C
     }
#elif MULTIFLUID //========================================================
      {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
#ifdef GAMCONS
      scalar gamma=1+rho/Ug[(e*N_F+4)*N_G+g];
#elif  GAMNCON
      scalar gamma=1+1.0/Ug[(e*N_F+4)*N_G+g];
#endif
      scalar p = (gamma-1)*(Et - 0.5*rho*(u*u+v*v));

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = rho*GX;
      s[(e*N_F+2)*N_G+g] = rho*GY; 
      s[(e*N_F+3)*N_G+g] = rho*(u*GX+v*GY);

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(Et+p,u);     // u(E+p)

      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(Et+p,v);     // v(E+p)

#ifdef GAMCONS
      s[(e*N_F+4)*N_G+g] = 0;
      f[((e*N_F+4)*N_G+g)*D+0] = flux_abc(rho,u,1/(gamma-1));// flux wrt x
      f[((e*N_F+4)*N_G+g)*D+1] = flux_abc(rho,v,1/(gamma-1));// flux wrt y
#elif  GAMNCON
      scalar vdotgradalpha = 0; // = -u*dalphadx-v*dalphady
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradalpha += 
	  -u*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
      }
      s[(e*N_F+4)*N_G+g] = vdotgradalpha; 
      f[((e*N_F+4)*N_G+g)*D+0] = 0;// flux wrt x
      f[((e*N_F+4)*N_G+g)*D+1] = 0;// flux wrt y
#endif

      // Mass fractions
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) s[(e*N_F+5+x)*N_G+g] = 0;				\
    f[((e*N_F+5+x)*N_G+g)*D+0] = flux_ab(Ug[(e*N_F+5+x)*N_G+g],u);	\
    f[((e*N_F+5+x)*N_G+g)*D+1] = flux_ab(Ug[(e*N_F+5+x)*N_G+g],v);	
#include "loop.h"
      }
#elif STIFFENED //========================================================
  {
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar gamma = 1+1.0/Ug[(e*N_F+4)*N_G+g];
      scalar beta  = Ug[(e*N_F+5)*N_G+g];
      scalar p = (gamma-1)*(Et - beta - 0.5*rho*(u*u+v*v));

      scalar vdotgradalpha = 0; // = -u*dalphadx-v*dalphady
      scalar vdotgradbeta = 0; // = -u*dbetadx-v*dbetady
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradalpha += 
	  -u*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
	vdotgradbeta += 
	  -u*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
      }

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = rho*GX;
      s[(e*N_F+2)*N_G+g] = rho*GY; 
      s[(e*N_F+3)*N_G+g] = rho*(u*GX+v*GY);
      s[(e*N_F+4)*N_G+g] = vdotgradalpha;
      s[(e*N_F+5)*N_G+g] = vdotgradbeta;
      
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(Et+p,u);     // u(E+p)
      f[((e*N_F+4)*N_G+g)*D+0] = 0;// flux wrt x
      f[((e*N_F+5)*N_G+g)*D+0] = 0;// flux wrt x

      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(Et+p,v);     // v(E+p)
      f[((e*N_F+4)*N_G+g)*D+1] = 0;// flux wrt y
      f[((e*N_F+5)*N_G+g)*D+1] = 0;// flux wrt x

      // Mass fractions
      int fcnt = 6;
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) s[(e*N_F+fcnt+x)*N_G+g] = 0;				\
    f[((e*N_F+fcnt+x)*N_G+g)*D+0] = flux_ab(Ug[(e*N_F+fcnt+x)*N_G+g],u);	\
    f[((e*N_F+fcnt+x)*N_G+g)*D+1] = flux_ab(Ug[(e*N_F+fcnt+x)*N_G+g],v); fcnt++;
#include "loop.h"
  }
#endif // end ifs on physics

#elif THREED
      //3D section added by PEJ 05/26/2017
#ifdef SCALARAD
      {
	//PEJ 05/29/2017: scalar advection-diffusion
	scalar rho = Ug[(e*N_F + 0)*N_G+g];
	//flux has one component per spatial dimension, so three in this case
	
	//Source term is zero
	s[(e*N_F+0)*N_G+g] = 0.0;
	
	//flux in x direction
	scalar Vx = constants::GLOBAL_VX;
	f[((e*N_F+0)*N_G+g)*D + 0] = flux_ab(rho, Vx); //this is Vx times rho.
	//flux in y direction
	scalar Vy = constants::GLOBAL_VY;
	f[((e*N_F+0)*N_G+g)*D + 1] = flux_ab(rho, Vy); //this is Vy times rho.
	//flux in z direction
	scalar Vz = constants::GLOBAL_VZ;
	f[((e*N_F+0)*N_G+g)*D + 2] = flux_ab(rho, Vz); //this is Vz times rho.
      }
#endif //end if for scalarAD. Hopefully, nothing else in the TWOD section activates now
#ifdef SINGLEFLUID
      {
	//printf("Activating 3D singlefluid element inviscid flux\n");
	//PEJ 09/25/2017: Inviscid fluxes, Navier-Stokes.
	//Get the conserved variables:
	scalar rho   = Ug[(e*N_F+0)*N_G+g];   
	scalar u     = Ug[(e*N_F+1)*N_G+g] / rho;  // (rho u / rho) = u
	scalar v     = Ug[(e*N_F+2)*N_G+g] / rho;  // (rho v / rho) = v
	scalar w     = Ug[(e*N_F+3)*N_G+g] / rho;  // (rho w / rho) = w
	scalar Et    = Ug[(e*N_F+4)*N_G+g];

	//Specific heat and pressure:
	scalar gamma = constants::GLOBAL_GAMMA;
	scalar p = (gamma-1)*(Et - 0.5*rho*(u*u + v*v + w*w));

	//Source term: For efficiensy, maybe get rid of this dependency eventually
	// Source term
	s[(e*N_F+0)*N_G+g] = 0;
	s[(e*N_F+1)*N_G+g] = 0;
	s[(e*N_F+2)*N_G+g] = 0;
	s[(e*N_F+3)*N_G+g] = 0;
	s[(e*N_F+4)*N_G+g] = 0;
      
	// Flux derive par rapport a x
	f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
	f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p); // rho*u^2 + p
	f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,v,u);   // rho*v*u
	f[((e*N_F+3)*N_G+g)*D+0] = flux_abc(rho,w,u);   // rho*w*u
	f[((e*N_F+4)*N_G+g)*D+0] = flux_ab(Et+p,u);  // u(E+p)
	
	// Flux derive par rapport a y
	f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
	f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);   // rho*u*v
	f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p); // rho*v^2 + p
	f[((e*N_F+3)*N_G+g)*D+1] = flux_abc(rho,w,v);   // rho*w*v
	f[((e*N_F+4)*N_G+g)*D+1] = flux_ab(Et+p,v);  // v(E+p)

	// Flux derive par rapport a z
	f[((e*N_F+0)*N_G+g)*D+2] = flux_ab(rho,w);      // rho*w     
	f[((e*N_F+1)*N_G+g)*D+2] = flux_abc(rho,u,w);   // rho*u*w
	f[((e*N_F+2)*N_G+g)*D+2] = flux_abc(rho,v,w);   // rho*v*w
	f[((e*N_F+3)*N_G+g)*D+2] = flux_ab2pc(rho,w,p); // rho*w^2+p
	f[((e*N_F+4)*N_G+g)*D+2] = flux_ab(Et+p,w);  // w(E+p)
      }
#endif


#endif // end ifs on dimensions
     
#ifdef USE_CPU
    }
#endif
  }
}

//=======================================================================



//==========================================================================
  arch_global void evaluate_q(int M_G, int M_T, scalar* q, scalar* UgF, scalar* normals){//, scalar* xyzf){
  /*!
    \brief Evaluate interface fluxes
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[out] q interface flux array
    \param[in] UgF interface solution (evaluated at gaussian nodes)
    \param[in] normals interface normals
  */
 
#ifdef USE_CPU
  int blk = 0;
  for(int t = 0; t < M_T; t++){
    scalar* buffer = new scalar[M_G*2*N_F]; // buffer holds F and ncterm (M_G x 2*N_F: [F,ncterm])
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
  int blk = threadIdx.z; // buffer needs room for all the elements in the block
  int t = blockIdx.x*blkT+blk;
  if ( t < M_T){
    int g = threadIdx.x;
    extern __shared__ scalar buffer[];    // buffer holds F and ncterm (blkT x M_G x 2*N_F: [F,ncterm])
#endif

      // Initialize the buffer
      int Fidx = (blk*M_G+g)*2*N_F;       // index for start of F
      int ncidx = (blk*M_G+g)*2*N_F+N_F;  // index for start of ncterm
      for(int k = 0; k < 2*N_F; k++) buffer[Fidx+k]  = 0;

     
      // Send the data to the Riemann solvers. These routines are stored in face_fluxes folder
#ifdef ONED
#ifdef SCALARAD //====================================================
      {
	//only doing rusanov flux for scalar AD, it is easier than upwinding
#ifdef RUS
	{
	  oned_scalarad_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif //end of rusanov case
#ifdef UPW
	{
	  //the exact upwind case:
	  oned_scalarad_upwind(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif //end of pure upwind case
#ifdef CEN
	{
	  //the exact upwind case:
	  oned_scalarad_central(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif
      }
#endif //end of scalarAD case
      
#ifdef PASSIVE //=========================================================
      {
#ifdef RUS
      oned_passive_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			   UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
			   UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
			   UgF[((t*N_F+4)*2+0)*M_G+g],                            // phincL
			   UgF[((t*N_F+4)*2+1)*M_G+g],                            // phincR
			   normals[t*D+0],                                        // nx
			   &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      oned_passive_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
		       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
		       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
		       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
		       UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
		       UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
		       UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
		       UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
		       UgF[((t*N_F+4)*2+0)*M_G+g],                            // phincL
		       UgF[((t*N_F+4)*2+1)*M_G+g],                            // phincR
		       normals[t*D+0],                                        // nx
		       &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      oned_passive_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
		       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
		       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
		       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
		       UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
		       UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
		       UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
		       UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
		       UgF[((t*N_F+4)*2+0)*M_G+g],                            // phincL
		       UgF[((t*N_F+4)*2+1)*M_G+g],                            // phincR
		       normals[t*D+0],                                        // nx
		       &buffer[Fidx],&buffer[ncidx]);
#endif // flux if
      }
#elif SINGLEFLUID //=========================================================
      {
#ifdef RUS
      oned_singlefluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      oned_singlefluid_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			   normals[t*D+0],                                        // nx
			   &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      oned_singlefluid_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			   normals[t*D+0],                                        // nx
			   &buffer[Fidx],&buffer[ncidx]);
#elif SLAU
      oned_singlefluid_slau2(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			     UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			     UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			     normals[t*D+0],                                        // nx
			     &buffer[Fidx],&buffer[ncidx]);

#endif // flux if
      }
#elif RADSINGLEFLUID //=================================================
      {
	oned_radsinglefluid_slau2(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				  UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
				  UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
				  UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
				  UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
				  UgF[((t*N_F+3)*2+0)*M_G+g],                            // CL
				  UgF[((t*N_F+3)*2+1)*M_G+g],                            // CR
				  normals[t*D+0],                                        // nx
				  &buffer[Fidx],&buffer[ncidx]);
      }
#elif MULTIFLUID //=========================================================
      {
#ifdef RUS
      // Coordinate of face if you want it...
      //scalar x = xyzf[(t*M_G+g)*D+0];
      oned_multifluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			      UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			      UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			      UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			      UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			      UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			      UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			      UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)              UgF[((t*N_F+4+x)*2+0)*M_G+g], UgF[((t*N_F+4+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			      normals[t*D+0],                                        // nx
			      //x,
			      &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      oned_multifluid_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			  UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			  UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			  UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			  UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			  UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			  UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)          UgF[((t*N_F+4+x)*2+0)*M_G+g], UgF[((t*N_F+4+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			  normals[t*D+0],                                        // nx
			  &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      oned_multifluid_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			  UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			  UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			  UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			  UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			  UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			  UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)          UgF[((t*N_F+4+x)*2+0)*M_G+g], UgF[((t*N_F+4+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			  normals[t*D+0],                                        // nx
			  &buffer[Fidx],&buffer[ncidx]);
#endif // flux if
      }
#elif STIFFENED //=========================================================
      {
#ifdef RUS
      oned_stiffened_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			     UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			     UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			     UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			     UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
			     UgF[((t*N_F+4)*2+0)*M_G+g],                            // betaL
			     UgF[((t*N_F+4)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)             UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			     normals[t*D+0],                                        // nx
			     &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      oned_stiffened_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			 UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			 UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			 UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			 UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			 UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			 UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			 UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
			 UgF[((t*N_F+4)*2+0)*M_G+g],                            // betaL
			 UgF[((t*N_F+4)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)         UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			 normals[t*D+0],                                        // nx
			 &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      oned_stiffened_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			 UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			 UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			 UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			 UgF[((t*N_F+2)*2+0)*M_G+g],                            // EtL
			 UgF[((t*N_F+2)*2+1)*M_G+g],                            // EtR
			 UgF[((t*N_F+3)*2+0)*M_G+g],                            // alphaL
			 UgF[((t*N_F+3)*2+1)*M_G+g],                            // alphaR
			 UgF[((t*N_F+4)*2+0)*M_G+g],                            // betaL
			 UgF[((t*N_F+4)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)         UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			 normals[t*D+0],                                        // nx
			 &buffer[Fidx],&buffer[ncidx]);
#endif // flux if
      }
#endif // physics if
  
#elif TWOD      //////////////////////////////////BEGINNING 2D SECTION
      
#ifdef SCALARAD
      {
#ifdef RUS
	{
	  twod_scalarad_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				normals[t*D+1],                                        // ny
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif
	
#ifdef CEN
	{
	  twod_scalarad_central(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				normals[t*D+1],                                        // ny
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif

#ifdef UPW
	{
	  twod_scalarad_upwind(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				normals[t*D+0],                                        // nx
				normals[t*D+1],                                        // ny
				&buffer[Fidx],&buffer[ncidx]);
	}
#endif	
      }
#endif

#ifdef PASSIVE
      {
#ifdef RUS
      twod_passive_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			   UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			   UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			   UgF[((t*N_F+4)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
			   UgF[((t*N_F+4)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
			   UgF[((t*N_F+5)*2+0)*M_G+g],                            // phincL
			   UgF[((t*N_F+5)*2+1)*M_G+g],                            // phincR
			   normals[t*D+0],                                        // nx
			   normals[t*D+1],                                        // ny
			   &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      twod_passive_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
		       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
		       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
		       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
		       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
		       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
		       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
		       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
		       UgF[((t*N_F+4)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
		       UgF[((t*N_F+4)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
		       UgF[((t*N_F+5)*2+0)*M_G+g],                            // phincL
		       UgF[((t*N_F+5)*2+1)*M_G+g],                            // phincR
		       normals[t*D+0],                                        // nx
		       normals[t*D+1],                                        // ny
		       &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      twod_passive_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
		       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
		       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
		       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
		       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
		       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
		       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
		       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
		       UgF[((t*N_F+4)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // phicL
		       UgF[((t*N_F+4)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // phicR
		       UgF[((t*N_F+5)*2+0)*M_G+g],                            // phincL
		       UgF[((t*N_F+5)*2+1)*M_G+g],                            // phincR
		       normals[t*D+0],                                        // nx
		       normals[t*D+1],                                        // ny
		       &buffer[Fidx],&buffer[ncidx]);
      
#endif // flux if
      }
#elif SINGLEFLUID //=========================================================
      {
#ifdef RUS
      twod_singlefluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       normals[t*D+1],                                        // ny
			       &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      twod_singlefluid_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			   UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			   UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			   normals[t*D+0],                                        // nx
			   normals[t*D+1],                                        // ny
			   &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      twod_singlefluid_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			   UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			   UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			   UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			   UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			   normals[t*D+0],                                        // nx
			   normals[t*D+1],                                        // ny
			   &buffer[Fidx],&buffer[ncidx]);
#elif ROP
      twod_singlefluid_roePHIL(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       normals[t*D+1],                                        // ny
			       &buffer[Fidx],&buffer[ncidx]);
#elif SLAU
      twod_singlefluid_slau2(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			     UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			     UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			     UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			     UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			     normals[t*D+0],                                        // nx
			     normals[t*D+1],                                        // ny
			     &buffer[Fidx],&buffer[ncidx]);
      
#endif // flux if
      }
#elif RADSINGLEFLUID //================================================
      {
#ifdef SLAU
        twod_radsinglefluid_slau2(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       normals[t*D+1],                                        // ny
			       &buffer[Fidx],&buffer[ncidx]);
#elif RUS
	twod_radsinglefluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       normals[t*D+1],                                        // ny
			       &buffer[Fidx],&buffer[ncidx]);
#elif ROE
	twod_radsinglefluid_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			       UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			       UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			       UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			       UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			       UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			       UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			       UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			       normals[t*D+0],                                        // nx
			       normals[t*D+1],                                        // ny
			       &buffer[Fidx],&buffer[ncidx]);

#endif //flux if
      }
#elif MULTIFLUID //=========================================================
      {
#ifdef RUS
      twod_multifluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			      UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			      UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			      UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			      UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			      UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			      UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			      UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			      UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			      UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)              UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			      normals[t*D+0],                                        // nx
			      normals[t*D+1],                                        // ny
			      &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      twod_multifluid_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			  UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			  UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			  UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			  UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			  UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			  UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			  UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			  UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)          UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			  normals[t*D+0],                                        // nx
			  normals[t*D+1],                                        // ny
			  &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      twod_multifluid_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			  UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			  UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			  UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			  UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			  UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			  UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			  UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			  UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)          UgF[((t*N_F+5+x)*2+0)*M_G+g], UgF[((t*N_F+5+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			  normals[t*D+0],                                        // nx
			  normals[t*D+1],                                        // ny
			  &buffer[Fidx],&buffer[ncidx]);
#endif // flux if
      }
#elif STIFFENED //=========================================================
      {
#ifdef RUS
      twod_stiffened_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			     UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			     UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			     UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			     UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			     UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			     UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
			     UgF[((t*N_F+5)*2+0)*M_G+g],                            // betaL
			     UgF[((t*N_F+5)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)             UgF[((t*N_F+6+x)*2+0)*M_G+g], UgF[((t*N_F+6+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			     normals[t*D+0],                                        // nx
			     normals[t*D+1],                                        // ny
			     &buffer[Fidx],&buffer[ncidx]);
#elif HLL
      twod_stiffened_hll(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			 UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			 UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			 UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			 UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			 UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			 UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			 UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			 UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			 UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
			 UgF[((t*N_F+5)*2+0)*M_G+g],                            // betaL
			 UgF[((t*N_F+5)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)         UgF[((t*N_F+6+x)*2+0)*M_G+g], UgF[((t*N_F+6+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			 normals[t*D+0],                                        // nx
			 normals[t*D+1],                                        // ny
			 &buffer[Fidx],&buffer[ncidx]);
#elif ROE
      twod_stiffened_roe(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
			 UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
			 UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
			 UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
			 UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
			 UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
			 UgF[((t*N_F+3)*2+0)*M_G+g],                            // EtL
			 UgF[((t*N_F+3)*2+1)*M_G+g],                            // EtR
			 UgF[((t*N_F+4)*2+0)*M_G+g],                            // alphaL
			 UgF[((t*N_F+4)*2+1)*M_G+g],                            // alphaR
			 UgF[((t*N_F+5)*2+0)*M_G+g],                            // betaL
			 UgF[((t*N_F+5)*2+1)*M_G+g],                            // betaR
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x)         UgF[((t*N_F+6+x)*2+0)*M_G+g], UgF[((t*N_F+6+x)*2+1)*M_G+g],// rhoL*YL, rhoR*YR
#include "loop.h"
			 normals[t*D+0],                                        // nx
			 normals[t*D+1],                                        // ny
			 &buffer[Fidx],&buffer[ncidx]);
#endif // flux if
      }
#endif // physics if

#elif THREED
#ifdef SCALARAD
      {
#ifdef RUS
	{
	  threed_scalarad_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				  normals[t*D+0],                                        // nx
				  normals[t*D+1],                                        // ny
				  normals[t*D+2],                                       //nz
				  &buffer[Fidx],&buffer[ncidx]);
	}
#endif //rusanov if
#ifdef UPW
	{
	  threed_scalarad_upwind(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				 UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				 normals[t*D+0],                                        // nx
				 normals[t*D+1],                                        // ny
				 normals[t*D+2],                                       //nz
				 &buffer[Fidx],&buffer[ncidx]);
	} //end upwind if
#endif
#ifdef CEN
	{
	  threed_scalarad_central(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				  UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				  normals[t*D+0],                                        // nx
				  normals[t*D+1],                                        // ny
				  normals[t*D+2],                                       //nz
				  &buffer[Fidx],&buffer[ncidx]);
	} //end central if
#endif //central flux if
	
      } //end scalarAD if
#endif //scalarAD if
      
#ifdef SINGLEFLUID
      {
	//WELCOME TO THE THUNDERDOME
#ifdef CEN
	{
	  //printf("Activating 3D singlefluid CENTRAL interface inviscid flux\n");
	  threed_singlefluid_central(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
				     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
				     UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
				     UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
				     UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vzL
				     UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vzR
				     UgF[((t*N_F+4)*2+0)*M_G+g],                            // EtL
				     UgF[((t*N_F+4)*2+1)*M_G+g],                            // EtR
				     normals[t*D+0],                                        // nx
				     normals[t*D+1],                                        // ny
				     normals[t*D+2],                                        // nz
				     &buffer[Fidx],&buffer[ncidx]);
	}
#endif //end if for central flux
#ifdef RUS
	{
	  threed_singlefluid_rusanov(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				     UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				     UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
				     UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
				     UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
				     UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
				     UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vzL
				     UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vzR
				     UgF[((t*N_F+4)*2+0)*M_G+g],                            // EtL
				     UgF[((t*N_F+4)*2+1)*M_G+g],                            // EtR
				     normals[t*D+0],                                        // nx
				     normals[t*D+1],                                        // ny
				     normals[t*D+2],                                        // nz
				     &buffer[Fidx],&buffer[ncidx]);
	}
#endif //end if for rusanov flux
#ifdef SLAU
	{
	  threed_singlefluid_slau2(UgF[((t*N_F+0)*2+0)*M_G+g],                            // rhoL
				   UgF[((t*N_F+0)*2+1)*M_G+g],                            // rhoR
				   UgF[((t*N_F+1)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vxL
				   UgF[((t*N_F+1)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vxR
				   UgF[((t*N_F+2)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vyL
				   UgF[((t*N_F+2)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vyR
				   UgF[((t*N_F+3)*2+0)*M_G+g]/UgF[((t*N_F+0)*2+0)*M_G+g], // vzL
				   UgF[((t*N_F+3)*2+1)*M_G+g]/UgF[((t*N_F+0)*2+1)*M_G+g], // vzR
				   UgF[((t*N_F+4)*2+0)*M_G+g],                            // EtL
				   UgF[((t*N_F+4)*2+1)*M_G+g],                            // EtR
				   normals[t*D+0],                                        // nx
				   normals[t*D+1],                                        // ny
				   normals[t*D+2],                                        // nz
				   &buffer[Fidx],&buffer[ncidx]);
	}
#endif //end if for slau2 flux
      }
#endif //end if for singlefluid physics
 
#endif // dimension if

      // Apply the fluxes
      for(int fc = 0; fc < N_F; fc++){
	q[((t*N_F+fc)*2+0)*M_G+g] =-buffer[Fidx+fc] + buffer[ncidx+fc];
	q[((t*N_F+fc)*2+1)*M_G+g] = buffer[Fidx+fc] + buffer[ncidx+fc];
      }
      
#ifdef USE_CPU
  } // end loop on g
    delete[] buffer;
#endif
    } // end loop on t
  } 
  


  arch_global void evaluate_sf_rad(int N_G, int N_E, int order, scalar* DivMax, scalar Beta_S, scalar Mew_S, scalar* CsMax, scalar* LamMax, scalar eps_gen, scalar* ADeps, scalar Cthresh, scalar* elemCmax, int* sensor, scalar* invJac, scalar* Ug, scalar* dUg_phys, scalar* s, scalar* f){//, scalar* xyz){
    /*!
      \Evaluates artificial viscosity flux, adds it to interior inviscid flux
      \param[in] N_G number of gaussian nodes per element
      \param[in] N_E number of elements
      \param[in] order DG solution order
      \param[in] DivMax maximum super-divergence (del dot velocity, with adjusted formula to maximize size) in the entire domain
      \param[in] Beta_S global free paramter for setting strength of AD
      \param[in] Mew_S global free paramter for setting spread speed of AD
      \param[in] CsMax The maximum value of C anywhere in the domain
      \param[in] LamMax maximum directional wavespeeds anywhere in the domain {(u+a),(v+a),(w+a)}
      \param[in] eps_gen small value to go in denominators, prevent division by zero
      \param[in] ADeps each element's (h/p)^2 value to scalar artificial dissipation, different for each direction
      \param[in] elemCmax the maximum C value over the element
      \param[in] sensor integer value indicating whether a shock/contact is in the region
      \param[in] invJac inverse jacobian of the elements
      
      \param[in] Ug solution (evaluated at gaussian nodes)
      \param[in] dUg_phys gradient of solution wrt physical coordinates (evaluated at gaussian nodes)
      \param[in and out] s element source term array
      \param[in and out] f element flux array, incremented by viscous flux
    */
   
    //template: dUg[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f]; a=dimensional component, g=quadrature node, f=field

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    //Run this routine ONLY if the element has some AD relevance
    //elemCmax is a c value on the element, and sensor[e] comes from Marc's sensor array.
    //Note that even if sensor is off, I can still treat the element
    //to hopefully dissipate (actually sink, since I'm using source term) C away.
    if (elemCmax[e] > Cthresh/*eps_gen*/ || sensor[e] > 0)
      {
	//	printf("sfRad executing element %d: elemCmax=%f, sensor=%d, DivMax=(%f,%f), LamMax=(%f,%f) \n",e,elemCmax[e],sensor[e],DivMax[0],DivMax[1],LamMax[0],LamMax[1]);
	int sensorYN = 0;
	//Forcing active whenever sensor is nonzero.
	//Also, if Cthresh is negative, that means I want AD forcing
	//to always be active.
	if (sensor[e] > 0 || Cthresh < 0){sensorYN = 1;}
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

#ifdef RADSINGLEFLUID //=========================================================
    {
    /* Marc's array technique, quoted for reference:
       scalar rho   = Ug[(e*N_F+0)*N_G+g];   
       scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
       scalar Et    = Ug[(e*N_F+2)*N_G+g];
    */
    /*
    //Step 1: Get the artificial dissipation coefficient, Kappa. Needs to be matrix to eventually add streamwise-type dependence
    scalar Kappa[D][D]; //yes, different at each 
    //How do I get this? Requires some global constants that will need to be fed in to this routine
    
    //get gas constant, ratio of specific heats, Prandtl number, and specific heat cp
    scalar Rgas =   constants:: GLOBAL_RGAS; //Fetch the Gas Constant, usually 287.15
    scalar gamma =  constants::GLOBAL_GAMMA;
    scalar Pran =   constants::GLOBAL_PRAN; //Fetch Prandtl Number
    scalar CpGas =  constants::GLOBAL_CPGAS; //specific heat (at constant pressure, I think)
    */
    // scalar eps_gen = pow(10,-12); //to avoid division by zero
    //Get the conserved variables (need velocity gradient):
    scalar rho   = Ug[e*N_F*N_G + 0*N_G + g];
    scalar mox   = Ug[e*N_F*N_G + 1*N_G + g];
    //scalar Et    = Ug[e*N_F*N_G + 2*N_G + g];
    
    //Cs is the only field variable that I need:
    scalar Cs    = Ug[e*N_F*N_G + 3*N_G + g];
    scalar eps_S = ADeps[e];

    //Import gradients od conserved variables + C variable: Subscript means derivative here
    /*
      //If imported gradient is wrt reference coordinates, must multiply by inverse Jacobian
    scalar rho_x = dUg_ref[(e*N_F+0)*N_G+g] * invJac[e*N_G+g];
    scalar mox_x = dUg_ref[(e*N_F+1)*N_G+g] * invJac[e*N_G+g];
    scalar Et_x = dUg_ref[(e*N_F+2)*N_G+g] * invJac[e*N_G+g];
    scalar Cs_x = dUg_ref[(e*N_F+3)*N_G+g] * invJac[e*N_G+g];
    */

    //Sigma (dUinteg_phys) uses different organization than dUg_ref(dUinteg)
    scalar rho_x = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0];
    scalar mox_x = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 1];
    scalar Et_x =  dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 2];
    scalar Cs_x =  dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 3];
    /*
    scalar rho_x = dUg_phys[(e*N_F+0)*N_G+g];
    scalar mox_x = dUg_phys[(e*N_F+1)*N_G+g];
    scalar Et_x = dUg_phys[(e*N_F+2)*N_G+g];
    scalar Cs_x = dUg_phys[(e*N_F+3)*N_G+g];
    */
    //velocity gradient (necessary for source term)
    scalar u_x = (rho*mox_x - mox*rho_x) / (rho*rho);

    //Set strength of AD: DivMax, CsMax are global constants in space: eps_S should be constant in time per element
    scalar factor = Cs * Beta_S * (eps_S*eps_S) * (DivMax[0]/(CsMax[0] + eps_gen));
    //factor = 0.5*factor*(u_x-fabs(u_x));
    // printf("Evaluate_sf_rad: e=%d,g=%d, rho=%f, mox=%f, Cs=%f, rho_x=%f, mox_x=%f, Et_x=%f, Cs_x=%f, u_x=%f, DIF_FACTOR=%f\n", e, g, rho,mox,Cs,rho_x,mox_x,Et_x,Cs_x,u_x,factor);    
    //Set the AD flux
    scalar Qs_rho = factor*rho_x;
    scalar Qs_mox = factor*mox_x;
    scalar Qs_Et = factor*Et_x;
    //scalar Qs_C = LamMax[0] * Mew_S * eps_S * Cs_x;
    //10/24: Local wavespeed approach for C, requires a bit more information:
    scalar gamma =  constants::GLOBAL_GAMMA;
    scalar Et    = Ug[e*N_F*N_G + 2*N_G + g];
    scalar p = (gamma-1)*(Et - 0.5*mox*mox/rho);
    scalar LamLocal = sqrt(gamma*p/rho) + fabs(mox/rho);
    scalar Qs_C = LamLocal * Mew_S * eps_S * Cs_x;
    // Source term also requires alteration for Reisner AD approach (to evolve C)
    //Calculate the G force (Reisner gives two formulae, I'm using the simpler one)
    
    //Modify source term in 4th field variable (C)
    //Reisner original approach: Gforce grows linearly
    //with fabs(u_x)/DivMax.
    /*
    scalar Gforce = fabs(u_x) / (DivMax[0] + eps_gen);
    s[(e*N_F+3)*N_G+g] = LamMax[0] * (Gforce - Cs/eps_S); 
    */
    //scalar Gforce = fmax(pow(fabs(u_x)/(DivMax[0]+eps_gen) , order+1), 1.0); //grow with square instead
    //    scalar Gforce = fabs(u_x) / (DivMax[0] + eps_gen); //original approach
    //Cool idea: set Gforce to be nonzero only when sensor
    //is triggered inside element, otherwise set to zero.
    //Then, the forcing term's solution is a decaying exponential, very nice.
    //Must be capped at 1.0 because DivMax and u_x are calculated using different methods,
    //So DivMax likely under-predicts the divergence.
    scalar Gforce = fmin(fabs(u_x) / (DivMax[0] + eps_gen), 1.0) * (sensorYN + 0.0);
    //    scalar Gforce = fabs(u_x) / (DivMax[0] + eps_gen) * (sensorYN + 0.0);
    //   printf("Evaluate_sf_rad: e=%d,g=%d, fabs(u_x)=%f, DivMax[0]=%f, sensorYN=%d, Gforce=%f\n", e, g, fabs(u_x), DivMax[0], sensorYN, Gforce);
    /*
    if (Gforce < 0.5)
      {
	Gforce = 0.0;
      }
    */
    s[(e*N_F+3)*N_G+g] = LamMax[0] * (Gforce - Cs/eps_S); 
    //    printf("Evaluate_sf_rad: e=%d,g=%d, Cs forcing = %f\n", e, g, s[(e*N_F+3)*N_G+g]);
    //  printf("Evaluate_sf_rad: e=%d,g=%d, Qs(rho,mox,Et,C) = (%f,%f,%f,%f)\n", e, g, Qs_rho,Qs_mox,Qs_Et,Qs_C);
    
    // Flux: the viscous contribution is subtracted from the previously calculated inviscid contribution
#ifdef USE_CPU
    f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - Qs_rho;
    f[((e*N_F+1)*N_G+g)*D+0] = f[((e*N_F+1)*N_G+g)*D+0] - Qs_mox; 
    f[((e*N_F+2)*N_G+g)*D+0] = f[((e*N_F+2)*N_G+g)*D+0] - Qs_Et;
    f[((e*N_F+3)*N_G+g)*D+0] = f[((e*N_F+3)*N_G+g)*D+0] - Qs_C;
    
#endif
#ifdef USE_GPU
    scalar addend_rho = -Qs_C;
    scalar addend_mox = -Qs_mox;
    scalar addend_en = -Qs_Et;
    scalar addend_C = -Qs_C;
    AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_rho);
    AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+0], addend_mox);
    AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+0], addend_en);
    AtomicAdd_vp(&f[((e*N_F+3)*N_G+g)*D+0], addend_C);
#endif
      
    }
#endif // end ifs on physics
    
#elif TWOD //============BEGIN 2D SECTION

#ifdef RADSINGLEFLUID
    {

      //For the gradients: I need rho and momentum to calculate
      //velocity gradients, but energy can be ignored.
      //However, gradients of all CONSERVED variables are needed.
      /*A thought: Eventually, maybe change the artificial
	dissipation applied to different fields. The first thing
	I can think of is that if pressure is much higher than
	density, then there should be extra dissipation,
	by some factor like p/rho, in the energy flux
      */
      //printf("ENTERED 2D evaluate_sf_rad procedure\n");
      //Conserved field variables:
      scalar rho   = Ug[(e*N_F+0)*N_G+g];
      scalar mox   = Ug[(e*N_F+1)*N_G+g];
      scalar moy   = Ug[(e*N_F+2)*N_G+g];
      //No need for energy
      scalar Cx    = Ug[(e*N_F+4)*N_G+g]; //(field 3 is energy)
      scalar Cy    = Ug[(e*N_F+5)*N_G+g]; //(field 3 is energy)

      //Gradients of conserved variables: subscript means derivative
      /*
      scalar rho_x = 0.0; scalar rho_y = 0.0;
      scalar mox_x = 0.0; scalar mox_y = 0.0;
      scalar moy_x = 0.0; scalar moy_y = 0.0;
      scalar Eg_x  = 0.0; scalar Eg_y = 0.0;
      scalar Cx_x  = 0.0; scalar Cx_y = 0.0;
      scalar Cy_x  = 0.0; scalar Cy_y = 0.0;
      
      //If imported gradient is wrt reference coordinates, must mulitply by inverse Jacobian
      //10/23/2017: I've decided that for invJac, the argument should be (alpha*D + a)
      for (int alpha = 0; alpha < D; alpha++)
	{
	  //both xi and eta derivates factor in to x and y derivatves in general.
	  //dUg_ref index alpha hits reference coordinates while 0,1 in invJac hit phys coords.
	  scalar dx = invJac[e*N_G*D*D + g*D*D + alpha*D + 0];
	  scalar dy = invJac[e*N_G*D*D + g*D*D + alpha*D + 1];

	  rho_x += dUg_ref[((e*N_F+0)*N_G+g)*D + alpha] * dx; 
	  rho_y += dUg_ref[((e*N_F+0)*N_G+g)*D + alpha] * dy;

	  mox_x += dUg_ref[((e*N_F+1)*N_G+g)*D + alpha] * dx;
	  mox_y += dUg_ref[((e*N_F+1)*N_G+g)*D + alpha] * dy;

	  moy_x += dUg_ref[((e*N_F+2)*N_G+g)*D + alpha] * dx;
	  moy_y += dUg_ref[((e*N_F+2)*N_G+g)*D + alpha] * dy;

	  Eg_x += dUg_ref[((e*N_F+3)*N_G+g)*D + alpha] * dx;
	  Eg_y += dUg_ref[((e*N_F+3)*N_G+g)*D + alpha] * dy;

	  Cx_x += dUg_ref[((e*N_F+4)*N_G+g)*D + alpha] * dx;
	  Cx_y += dUg_ref[((e*N_F+4)*N_G+g)*D + alpha] * dy;

	  Cy_x += dUg_ref[((e*N_F+5)*N_G+g)*D + alpha] * dx;
	  Cy_y += dUg_ref[((e*N_F+5)*N_G+g)*D + alpha] * dy;
	}
      */
      //Gradients of conserved variables:
      scalar rho_x = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 0]; scalar rho_y = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 0];
      scalar mox_x = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 1]; scalar mox_y = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 1];
      scalar moy_x = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 2]; scalar moy_y = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 2];
      scalar Eg_x  = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 3]; scalar Eg_y  = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 3];
      scalar Cx_x  = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 4]; scalar Cx_y  = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 4];
      scalar Cy_x  = dUg_phys[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 5]; scalar Cy_y  = dUg_phys[e*D*N_G*N_F + 1*N_G*N_F + g*N_F + 5];

      //Get velocity gradients so that we can properly scale C forcing
      scalar Inv_RhoSq = 1.0 / (rho*rho);
      scalar u_x = (rho*mox_x - mox*rho_x) * Inv_RhoSq;
      scalar v_y = (rho*moy_y - moy*rho_y) * Inv_RhoSq;
      
      //    printf("Evaluate_sf_rad: e=%d,g=%d, rho=%f, mox=%f, Cx=%f, rho_x=%f, mox_x=%f, Et_x=%f, Cx_x=%f, u_x=%f\n", e, g, rho,mox,Cx,rho_x,mox_x,Eg_x,Cx_x,u_x);    
      //Get the hp epsilons from global storage:
      scalar eps_S_x = ADeps[e*D + 0];
      scalar eps_S_y = ADeps[e*D + 1];
      
      //Get the Kappa diffusion matrix (it's a diagonal DxD matrix)
      //Here, It would be good to introduce a ratio between KappaX and KappaY;
      //if all of the C forcing is due to action in the x direction,
      //then KappaY should tend to zero.
      scalar KappaX = Cx * Beta_S * (eps_S_x*eps_S_x) * DivMax[0] / (CsMax[0] + eps_gen);
      scalar KappaY = Cy * Beta_S * (eps_S_y*eps_S_y) * DivMax[1] / (CsMax[1] + eps_gen);
      //Alternative idea: denominator is sum of max Cx and max Cy, to avoid nasty cross-diffusion
      
      //Get the G forcing:
      //Ideally, G forcing in each direction would use more detailed sensor information.
      //I leave this task to a future student.
      scalar DivMag = sqrt(DivMax[0]*DivMax[0] + DivMax[1]*DivMax[1]);
      scalar GforceX = fmin(fabs(u_x) / (DivMag + eps_gen),1.0) * (sensorYN + 0.0); 
      scalar GforceY = fmin(fabs(v_y) / (DivMag + eps_gen),1.0) * (sensorYN + 0.0); 
      //scalar LAM = fmax(LamMax[0], LamMax[1]);
      //scalar LAM = LamMax[0] + LamMax[1];
      //Adjust source term to force the C equations:
      //    printf("\t\tLamMax[0]=%f, LamMax[1]=%f, GforceX=%f, GforceY=%f,eps_S_x=%f,eps_S_y=%f\n",LamMax[0], LamMax[1], GforceX, GforceY, eps_S_x, eps_S_y);
      s[(e*N_F+4)*N_G+g] = LamMax[0] * (GforceX - Cx/eps_S_x); //Cx forcing
      s[(e*N_F+5)*N_G+g] = LamMax[1] * (GforceY - Cy/eps_S_y); //Cy forcing
      //s[(e*N_F+4)*N_G+g] = LAM * ( 0.5*(GforceX+GforceY) - C/fmin(eps_S_x,eps_S_y));
      //printf("\t\tCx source = %f, Cy source = %f\n", s[(e*N_F+4)*N_G+g],s[(e*N_F+5)*N_G+g]);
      
      //Use gradients to calculate AD viscous flux
      scalar Qs_rho_x = KappaX * rho_x; scalar Qs_rho_y = KappaY * rho_y;
      scalar Qs_mox_x = KappaX * mox_x; scalar Qs_mox_y = KappaY * mox_y;
      scalar Qs_moy_x = KappaX * moy_x; scalar Qs_moy_y = KappaY * moy_y;
      scalar Qs_Eg_x  = KappaX * Eg_x;  scalar Qs_Eg_y  = KappaY * Eg_y;
      //For C flux: Cx only diffuses in x direction, Cy only diffuses in y
      //scalar Qs_Cx_x = LamMax[0] * eps_S_x * Cx_x; scalar Qs_Cx_y = 0.0;
      //scalar Qs_Cy_x = 0.0;                        scalar Qs_Cy_y = LamMax[1] * eps_S_y * Cy_y;
      //scalar Qs_c_x = LamMax[0] * eps_S_x * C_x; 
      //scalar Qs_c_y = LamMax[1] * eps_S_y * C_y;
      //scalar Qs_c_x = LAM * eps_S_x * C_x; 
      //scalar Qs_c_y = LAM * eps_S_y * C_y;
      //10/24: Local wavespeed approach for C, requires a bit more information:
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Et    = Ug[e*N_F*N_G + 3*N_G + g];
      scalar p = (gamma-1)*(Et - 0.5/rho*(mox*mox + moy*moy));
      scalar LamXLocal = sqrt(gamma*p/rho) + fabs(mox/rho);
      scalar LamYLocal = sqrt(gamma*p/rho) + fabs(moy/rho);
      scalar Qs_Cx_x = LamXLocal * Mew_S * eps_S_x * Cx_x; scalar Qs_Cx_y = 0.0;
      scalar Qs_Cy_y = LamYLocal * Mew_S * eps_S_y * Cy_y; scalar Qs_Cy_x = 0.0;

      //for CPU case, I can subtract the viscous contribution from flux.
      //Need to use atomic add in the GPU case
#ifdef USE_CPU
      //Flux in x direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+0] = f[((e*N_F+0)*N_G+g)*D+0] - Qs_rho_x; 
      f[((e*N_F+1)*N_G+g)*D+0] = f[((e*N_F+1)*N_G+g)*D+0] - Qs_mox_x;
      f[((e*N_F+2)*N_G+g)*D+0] = f[((e*N_F+2)*N_G+g)*D+0] - Qs_moy_x;
      f[((e*N_F+3)*N_G+g)*D+0] = f[((e*N_F+3)*N_G+g)*D+0] - Qs_Eg_x;
      f[((e*N_F+4)*N_G+g)*D+0] = f[((e*N_F+4)*N_G+g)*D+0] - Qs_Cx_x;
      f[((e*N_F+5)*N_G+g)*D+0] = f[((e*N_F+5)*N_G+g)*D+0] - Qs_Cy_x;
      
      // Flux in y direction: Subtract viscous contribution from inviscid flux
      f[((e*N_F+0)*N_G+g)*D+1] = f[((e*N_F+0)*N_G+g)*D+1] - Qs_rho_y;  
      f[((e*N_F+1)*N_G+g)*D+1] = f[((e*N_F+1)*N_G+g)*D+1] - Qs_mox_y;  
      f[((e*N_F+2)*N_G+g)*D+1] = f[((e*N_F+2)*N_G+g)*D+1] - Qs_moy_y; 
      f[((e*N_F+3)*N_G+g)*D+1] = f[((e*N_F+3)*N_G+g)*D+1] - Qs_Eg_y;
      f[((e*N_F+4)*N_G+g)*D+1] = f[((e*N_F+4)*N_G+g)*D+1] - Qs_Cx_y;
      f[((e*N_F+5)*N_G+g)*D+1] = f[((e*N_F+5)*N_G+g)*D+1] - Qs_Cy_y;
#endif
#ifdef USE_GPU
      scalar addend_x_rho = -Qs_rho_X;
      scalar addend_x_mox = -Qs_mox_X;
      scalar addend_x_moy = -Qs_moy_X;
      scalar addend_x_en = -Qs_en_X;
      scalar addend_x_cx = -Qs_Cx_X;
      scalar addend_x_cy = -Qs_Cy_X;

      scalar addend_y_rho = -Qs_rho_Y;
      scalar addend_y_mox = -Qs_mox_Y;
      scalar addend_y_moy = -Qs_moy_Y;
      scalar addend_y_en = -Qs_en_Y;
      scalar addend_y_cx = -Qs_Cx_Y;
      scalar addend_y_cy = -Qs_Cy_Y;


      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+0], addend_x_rho);
      AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+0], addend_x_mox);
      AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+0], addend_x_moy);
      AtomicAdd_vp(&f[((e*N_F+3)*N_G+g)*D+0], addend_x_en);
      AtomicAdd_vp(&f[((e*N_F+4)*N_G+g)*D+0], addend_x_cx);
      AtomicAdd_vp(&f[((e*N_F+5)*N_G+g)*D+0], addend_x_cy);
      
      AtomicAdd_vp(&f[((e*N_F+0)*N_G+g)*D+1], addend_y_rho);
      AtomicAdd_vp(&f[((e*N_F+1)*N_G+g)*D+1], addend_y_mox);
      AtomicAdd_vp(&f[((e*N_F+2)*N_G+g)*D+1], addend_y_moy);
      AtomicAdd_vp(&f[((e*N_F+3)*N_G+g)*D+1], addend_y_en);
      AtomicAdd_vp(&f[((e*N_F+4)*N_G+g)*D+1], addend_y_cx);
      AtomicAdd_vp(&f[((e*N_F+5)*N_G+g)*D+1], addend_y_cy);
      
#endif
    } //close bracket for 2D RADSINGLEFLUID
#endif //end if on RADSINGLEFLUID physics
#endif // end ifs on dimensions  
#ifdef USE_CPU
  }
#endif
      }// end  giant "if" loop
  } //end e loop
  } //End evalueate_sf_vis subroutine

  arch_global void evaluate_q_rad(int M_G, int M_T, scalar* q, scalar* UgF, scalar* gradCommon, scalar* DivMax, scalar Beta_S, scalar Mew_S, scalar* CsMax, scalar* LamMax, scalar eps_gen, scalar* ADepsF, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* normals){
  /*!
    \brief Evaluate interface fluxes for artificial dissipation
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[in and out] q interface flux array
    \param[in] UgF interface solution (evaluated at gaussian nodes)
    \param[in] gradCommon interface gradient wrt physical coordinates (naive gradient)
    \param[in] DivMax max divergence in domain
    \param[in] Beta_S global free paramter for setting strength of AD
    \param[in] Beta_S global free paramter for setting spread speed of AD
    \param[in] CsMax The maximum value of C anywhere in the domain
    \param[in] LamMax maximum wavespeed anywhere in the domain (u+a)
    \param[in] eps_gen small value to go in denominators, prevent division by zero
    \param[in] ADepsF Face average of each element's (h/p)^2 
    \param[in] elemCmax the maximum C value over each element
    \param[in] sensor integer value indicating whether a shock/contact is near each element
    \param[in] BR2_Map tells each interface which elements border it.
    \param[in] normals interface normals
  */
 
#ifdef USE_CPU
  int blk = 0;
  for(int t = 0; t < M_T; t++){
    int omA = BR2_Map[t*4 + 2*0 + 0];
    int omB = BR2_Map[t*4 + 2*1 + 0];
    if (elemCmax[omA] > Cthresh/*eps_gen*/ || sensor[omA] > 0 || elemCmax[omB] > Cthresh/*eps_gen*/ || sensor[omB] > 0)
      {
	//	printf("Executing evalueate_q_rad for t=%d: omA=%d, omB=%d\n",t,omA,omB);
	//only run this routine if an element bordering the interface
	//is in need of AD
    scalar* buffer = new scalar[M_G*2*N_F]; // buffer holds F and ncterm (M_G x 2*N_F: [F,ncterm])
    for(int g = 0; g < M_G; g++){
      
#elif USE_GPU
  int blk = threadIdx.z; // buffer needs room for all the elements in the block
  int t = blockIdx.x*blkT+blk;
  if ( t < M_T){
    int g = threadIdx.x;
    extern __shared__ scalar buffer[];    // buffer holds F and ncterm (blkT x M_G x 2*N_F: [F,ncterm])
#endif
      
      // Initialize the buffer
      int Fidx = (blk*M_G+g)*2*N_F;       // index for start of F
      int ncidx = (blk*M_G+g)*2*N_F+N_F;  // index for start of ncterm
      for(int k = 0; k < 2*N_F; k++) buffer[Fidx+k]  = 0;

     
#ifdef ONED
#ifdef RADSINGLEFLUID //=========================================================
      {
      //Very similar to element-wise population of AD flux;
      //each interface has been populated with gradient wrt phys coordinates,
      //I just need to plug in to flux formula.
      
      //Difference from UhCommon for viscous physics:
      //here, I take Ug directly from the inviscid physics
      //appraoch (be it simple or ibc), then take average
      //for interface state
      
      //Get the conserved Variables:
      scalar rho = 0.5*(UgF[((t*N_F+0)*2+0)*M_G+g] + UgF[((t*N_F+0)*2+1)*M_G+g]);
      scalar mox = 0.5*(UgF[((t*N_F+1)*2+0)*M_G+g] + UgF[((t*N_F+1)*2+1)*M_G+g]);
      scalar Cs   = 0.5*(UgF[((t*N_F+3)*2+0)*M_G+g] + UgF[((t*N_F+3)*2+1)*M_G+g]);
      
      //local (h/p)^2 value:
      scalar eps_S = ADepsF[t];
      
      //Import gradients: Subscript means derivative here
      scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
      scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
      scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
      scalar Cs_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 0];
      
      //u_x is veclocity gradient
      scalar u_x = (rho*mox_x - mox*rho_x)/(rho*rho);
      //      printf("Evaluate_q_rad: t=%d,g=%d, rho=%f, mox=%f, rho_x=%f, mox_x=%f, u_x=%f\n", t, g, rho,mox,rho_x,mox_x,u_x);
      //Set strength of AD: DivMax, CsMax are global constants in space: eps_S should be constant in time per element
      scalar factor = Cs * Beta_S * (eps_S*eps_S) * (DivMax[0]/(CsMax[0]+eps_gen));
      //if divergence is positive, set factor=0, no need for AD
      //factor = 0.5*factor*(u_x-fabs(u_x));
      //     printf("Evaluate_q_rad: t=%d,g=%d, rho=%f, mox=%f, Cs=%f, rho_x=%f, mox_x=%f, Et_x=%f, Cs_x=%f, u_x=%f, DIF_FACTOR=%f\n", t, g, rho,mox,Cs,rho_x,mox_x,Et_x,Cs_x,u_x,factor);    
      //interface normal:
      //scalar nx = -normals[t*D+0];
      scalar nx = -normals[t*D+0];

      //Set the AD flux
      scalar Qs_rho = nx * factor*rho_x;
      scalar Qs_mox = nx * factor*mox_x;
      scalar Qs_en = nx * factor*Et_x;
      //scalar Qs_C = nx * LamMax[0] * Mew_S * eps_S * Cs_x;
      //10/24/2017: Use local wavespeed for C diffusivity coeffien:
      scalar gamma =  constants::GLOBAL_GAMMA;
      scalar Et = 0.5*(UgF[((t*N_F+2)*2+0)*M_G+g] + UgF[((t*N_F+2)*2+1)*M_G+g]);
      scalar p = (gamma-1)*(Et - 0.5*mox*mox/rho);
      scalar LamLocal = sqrt(gamma*p/rho) + fabs(mox/rho);
      scalar Qs_C = nx * LamLocal * Mew_S * eps_S * Cs_x;
      //printf("Evaluate_q_rad: t=%d,g=%d, Qs(rho,mox,Et,C) = (%f,%f,%f,%f)\n", t, g, Qs_rho,Qs_mox,Qs_en,Qs_C);
      
      //Adjust the Q vector based on your findings:
      //The sign is based on the "Q dot n" term in surface residual
#ifdef USE_CPU
      q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qs_rho;
      q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qs_rho;
      
      q[((t*N_F+1)*2+0)*M_G+g] = q[((t*N_F+1)*2+0)*M_G+g] - Qs_mox;
      q[((t*N_F+1)*2+1)*M_G+g] = q[((t*N_F+1)*2+1)*M_G+g] + Qs_mox;
      
      q[((t*N_F+2)*2+0)*M_G+g] = q[((t*N_F+2)*2+0)*M_G+g] - Qs_en;
      q[((t*N_F+2)*2+1)*M_G+g] = q[((t*N_F+2)*2+1)*M_G+g] + Qs_en;
      
      q[((t*N_F+3)*2+0)*M_G+g] = q[((t*N_F+3)*2+0)*M_G+g] - Qs_C;
      q[((t*N_F+3)*2+1)*M_G+g] = q[((t*N_F+3)*2+1)*M_G+g] + Qs_C;
#endif
#ifdef USE_GPU
      scalar addend_A_rho = -Qs_rho;
      scalar addend_B_rho = Qs_rho;
      scalar addend_A_mox = -Qs_mox;
      scalar addend_B_mox = Qs_mox;
      scalar addend_A_en = -Qs_en;
      scalar addend_B_en = Qs_en;
      scalar addend_A_C = -Qs_C;
      scalar addend_B_C = Qs_C;
      
      AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A_rho);
      AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B_rho);
      AtomicAdd_vp(&q[((t*N_F+1)*2+0)*M_G+g], addend_A_mox);
      AtomicAdd_vp(&q[((t*N_F+1)*2+1)*M_G+g], addend_B_mox);
      AtomicAdd_vp(&q[((t*N_F+2)*2+0)*M_G+g], addend_A_en);
      AtomicAdd_vp(&q[((t*N_F+2)*2+1)*M_G+g], addend_B_en);
      AtomicAdd_vp(&q[((t*N_F+3)*2+0)*M_G+g], addend_A_C);
      AtomicAdd_vp(&q[((t*N_F+3)*2+1)*M_G+g], addend_B_C);
#endif
      } //end RADSINGLEFLUID bracket

#endif // physics if
#elif TWOD //=========BEGIN 2D REISNER AD FLUX ADJUSTMENT===
#ifdef RADSINGLEFLUID
      {
	//printf("Entered 2D q_rad function\n");
	//Difference from UhCommon for viscous physics:
	//here, I take Ug directly from the inviscid physics
	//appraoch (be it simple or ibc), then take average
	//for interface state
	
	//Get the C Variable:
	scalar Cx = 0.5*(UgF[((t*N_F+4)*2+0)*M_G+g] + UgF[((t*N_F+4)*2+1)*M_G+g]);
	scalar Cy = 0.5*(UgF[((t*N_F+5)*2+0)*M_G+g] + UgF[((t*N_F+5)*2+1)*M_G+g]);
	
	//Import gradients: Subscript means derivative here.
	//Unlike sf_rad function, these gradients are already
	//with respect to phys coordinates
	/*
	scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
	scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
	scalar moy_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
	scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 0];
	scalar Cx_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 0];
	scalar Cy_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 5*D + 0];
	*/
	scalar rho_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 0];
	scalar rho_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 0*D + 1];
	
	scalar mox_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 0];
	scalar mox_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 1*D + 1];
	
	scalar moy_x = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 0];
	scalar moy_y = gradCommon[t*M_G*N_F*D + g*N_F*D + 2*D + 1];
	
	scalar Et_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 0];
	scalar Et_y  = gradCommon[t*M_G*N_F*D + g*N_F*D + 3*D + 1];

	scalar Cx_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 0];
	scalar Cx_y  = gradCommon[t*M_G*N_F*D + g*N_F*D + 4*D + 1];

	scalar Cy_x  = gradCommon[t*M_G*N_F*D + g*N_F*D + 5*D + 0];
	scalar Cy_y  = gradCommon[t*M_G*N_F*D + g*N_F*D + 5*D + 1];

	//Get the face's hp epsilons from global storage:
	scalar eps_S_x = ADepsF[t*D + 0];
	scalar eps_S_y = ADepsF[t*D + 1];
	//	printf("q_rad: eps_S_x=%f, eps_S_y=%f\n", eps_S_x, eps_S_y);
	
	//Get the Kappa diffusion matrix (it's a diagonal DxD matrix)
	//THIS MUCH MATCH WITH HOW WE TREAT THE ELEMENT FLUX
	scalar KappaX = Cx * Beta_S * (eps_S_x*eps_S_x) * DivMax[0] / (CsMax[0] + eps_gen);
	scalar KappaY = Cy * Beta_S * (eps_S_y*eps_S_y) * DivMax[1] / (CsMax[1] + eps_gen);
	//scalar factor = Cs * Beta_S * (eps_S*eps_S) * (DivMax[0]/(CsMax[0]+eps_gen));
	//No source term to worry about on interfaces
	//Grab normals to prepare for interface thruflux calculation
	scalar nx = -normals[t*D+0];
	scalar ny = -normals[t*D+1];
	//scalar LAM = fmax(LamMax[0], LamMax[1]);
	//scalar LAM = LamMax[0] + LamMax[1];
	//Set the AD flux (nx * xFlux + ny*yFlux)
	scalar Qs_rho = nx * (KappaX*rho_x)               + ny * (KappaY*rho_y);
	scalar Qs_mox = nx * (KappaX*mox_x)               + ny * (KappaY*mox_y);
	scalar Qs_moy = nx * (KappaX*moy_x)               + ny * (KappaY*moy_y);
	scalar Qs_en  = nx * (KappaX*Et_x)                + ny * (KappaY*Et_y);
	//scalar Qs_C   = nx * (LamMax[0] * eps_S_x * C_x)   + ny * (LamMax[1] * eps_S_y * C_y);
	//scalar Qs_Cx   = nx * (LamMax[0] * eps_S_x * Cx_x) + 0.0;
	//scalar Qs_Cy   = 0.0                               + ny * (LamMax[1] * eps_S_y * Cy_y);
	//For C equation diffusivity, use local flow properties:
	scalar gamma =  constants::GLOBAL_GAMMA;

	scalar rho = 0.5*(UgF[((t*N_F+0)*2+0)*M_G+g] + UgF[((t*N_F+0)*2+1)*M_G+g]);
	scalar u   = 0.5*(UgF[((t*N_F+1)*2+0)*M_G+g] + UgF[((t*N_F+1)*2+1)*M_G+g]) / rho;
	scalar v   = 0.5*(UgF[((t*N_F+2)*2+0)*M_G+g] + UgF[((t*N_F+2)*2+1)*M_G+g]) / rho;
	scalar Et  = 0.5*(UgF[((t*N_F+3)*2+0)*M_G+g] + UgF[((t*N_F+3)*2+1)*M_G+g]);
	scalar p = (gamma-1)*(Et - 0.5*rho*(u*u + v*v));
	scalar LamXLocal = sqrt(gamma*p/rho) + fabs(u);
	scalar LamYLocal = sqrt(gamma*p/rho) + fabs(v);
	scalar Qs_Cx   = nx * (LamXLocal * Mew_S * eps_S_x * Cx_x);
	scalar Qs_Cy   = ny * (LamYLocal * Mew_S * eps_S_y * Cy_y);

	//Adjust the Q vector based on your findings:
	//The sign is based on the "Q dot n" term in surface residual
#ifdef USE_CPU
	q[((t*N_F+0)*2+0)*M_G+g] = q[((t*N_F+0)*2+0)*M_G+g] - Qs_rho;
	q[((t*N_F+0)*2+1)*M_G+g] = q[((t*N_F+0)*2+1)*M_G+g] + Qs_rho;
	
	q[((t*N_F+1)*2+0)*M_G+g] = q[((t*N_F+1)*2+0)*M_G+g] - Qs_mox;
	q[((t*N_F+1)*2+1)*M_G+g] = q[((t*N_F+1)*2+1)*M_G+g] + Qs_mox;

	q[((t*N_F+2)*2+0)*M_G+g] = q[((t*N_F+2)*2+0)*M_G+g] - Qs_moy;
	q[((t*N_F+2)*2+1)*M_G+g] = q[((t*N_F+2)*2+1)*M_G+g] + Qs_moy;
	
	q[((t*N_F+3)*2+0)*M_G+g] = q[((t*N_F+3)*2+0)*M_G+g] - Qs_en;
	q[((t*N_F+3)*2+1)*M_G+g] = q[((t*N_F+3)*2+1)*M_G+g] + Qs_en;
	
	q[((t*N_F+4)*2+0)*M_G+g] = q[((t*N_F+4)*2+0)*M_G+g] - Qs_Cx;
	q[((t*N_F+4)*2+1)*M_G+g] = q[((t*N_F+4)*2+1)*M_G+g] + Qs_Cx;

	q[((t*N_F+5)*2+0)*M_G+g] = q[((t*N_F+5)*2+0)*M_G+g] - Qs_Cy;
	q[((t*N_F+5)*2+1)*M_G+g] = q[((t*N_F+5)*2+1)*M_G+g] + Qs_Cy;
#endif
#ifdef USE_GPU
      scalar addend_A_rho = -Qs_rho;
      scalar addend_B_rho = Qs_rho;
      scalar addend_A_mox = -Qs_mox;
      scalar addend_B_mox = Qs_mox;
      scalar addend_A_moy = -Qs_moy;
      scalar addend_B_moy = Qs_moy;
      scalar addend_A_en = -Qs_en;
      scalar addend_B_en = Qs_en;
      scalar addend_A_Cx = -Qs_Cx;
      scalar addend_B_Cx = Qs_Cx;
      scalar addend_A_Cy = -Qs_Cy;
      scalar addend_B_Cy = Qs_Cy;
      
      AtomicAdd_vp(&q[((t*N_F+0)*2+0)*M_G+g], addend_A_rho);
      AtomicAdd_vp(&q[((t*N_F+0)*2+1)*M_G+g], addend_B_rho);
      AtomicAdd_vp(&q[((t*N_F+1)*2+0)*M_G+g], addend_A_mox);
      AtomicAdd_vp(&q[((t*N_F+1)*2+1)*M_G+g], addend_B_mox);
      AtomicAdd_vp(&q[((t*N_F+2)*2+0)*M_G+g], addend_A_moy);
      AtomicAdd_vp(&q[((t*N_F+2)*2+1)*M_G+g], addend_B_moy);
      AtomicAdd_vp(&q[((t*N_F+3)*2+0)*M_G+g], addend_A_en);
      AtomicAdd_vp(&q[((t*N_F+3)*2+1)*M_G+g], addend_B_en);
      AtomicAdd_vp(&q[((t*N_F+4)*2+0)*M_G+g], addend_A_Cx);
      AtomicAdd_vp(&q[((t*N_F+4)*2+1)*M_G+g], addend_B_Cx);
      AtomicAdd_vp(&q[((t*N_F+5)*2+0)*M_G+g], addend_A_Cy);
      AtomicAdd_vp(&q[((t*N_F+5)*2+1)*M_G+g], addend_B_Cy);
#endif
      }//End RADSINGLEFLUID bracket
#endif //endif for RADSINGLEFLUID physics

#endif // dimension if
      /*
	This here is Marc's procedure for handling interface
	fluxes that I haven't been able to figure out yet.
      // Apply the fluxes
      for(int fc = 0; fc < N_F; fc++){
	q[((t*N_F+fc)*2+0)*M_G+g] =-buffer[Fidx+fc] + buffer[ncidx+fc];
	q[((t*N_F+fc)*2+1)*M_G+g] = buffer[Fidx+fc] + buffer[ncidx+fc];
      }
      */
#ifdef USE_CPU
  } // end loop on g
    delete[] buffer;
#endif
      } //end  giant "if" structure to engage AD equations
  } // end loop on t
  } //End evaluate_q_rad routine
  

  arch_global void GrabCparameters(int N_G, int N_E, scalar* Ug, scalar* dUg_ref, scalar* invJac, scalar* DivMax, scalar* LamMax, scalar* CsMax, scalar* elemCmax)
  {
    /*!
      \brief Calculates some global parameters for Reisner artificial viscosity approach
      \param[in] N_G gaussian quadrature nodes per element
      \param[in] N_E number of elements in the domain
      \param[in] Ug solution at gaussian quadrature nodes
      \param[in] dUg_ref gradient wrt reference coordinates at quadrature nodes
      \param[in] invJac inverse of element Jacobian at each quadrature node
      \param[out] DivMax max components of divergence calculation
      \param[out] LamMax max wavespeed in each direction in entire domain
      \param[out] CsMax max C values in entire domain, each direction
      \param[out] elemCmax max C value per eleah element, no direction, just 1 per element
    */
    //Before anything else, initialize the outputs to zero
    
    for (int a = 0; a < D; a++)
      {
	DivMax[a] = 0.0;
	CsMax[a] = 0.0;
	LamMax[a] = 0.0;
      }
    for (int e = 0; e < N_E; e++)
      {
	elemCmax[e] = 0.0;
      }
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
      
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int g = threadIdx.x;
#endif
      
    //1D problem=======================================================
#ifdef ONED
#ifdef RADSINGLEFLUID
    {
    //get some field variables
    scalar rho   = Ug[(e*N_F+0)*N_G+g];   
    scalar mox   = Ug[(e*N_F+1)*N_G+g];   
    scalar u     = mox / rho;  // (rho u / rho) = u
    scalar Et    = Ug[(e*N_F+2)*N_G+g];
    scalar C     = Ug[(e*N_F+3)*N_G+g];
    scalar gamma = constants::GLOBAL_GAMMA;
    scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
    
    //get some gradients
    //Import gradients od conserved variables + C variable: Subscript means derivative here
    //Proper dUg syntax from Marc's evaluate_sf routine:
    //   s[(e*N_F+4)*N_G+g] = -u*dUg[(e*N_F+4)*N_G+g]*invJac[e*N_G+g];// = -u*dphincdx;
    scalar rho_x = dUg_ref[(e*N_F+0)*N_G+g] * invJac[e*N_G+g];
    scalar mox_x = dUg_ref[(e*N_F+1)*N_G+g] * invJac[e*N_G+g];
    //scalar Et_x =  dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 2] * invJac[e*N_G+g];
    //scalar C_x =   dUg[e*D*N_G*N_F + 0*N_G*N_F + g*N_F + 2] * invJac[e*N_G+g];
      
    scalar u_x = (rho*mox_x - mox*rho_x) / (rho*rho);
    //printf("GrabCParameters: e=%d,g=%d, rho=%f, mox=%f, rho_x=%f, mox_x=%f, u_x=%f\n", e, g, rho,mox,rho_x,mox_x,u_x);
    //Calculate the local divergenve
    scalar Div = u_x + 0.0 + 0.0;
    //Calculate the local wavespeed
    scalar upa = fabs(u) + sqrt(gamma*p/rho);
    //Update global maxima
    DivMax[0] = fmax(fabs(Div), DivMax[0]);
    CsMax[0] = fmax(C, CsMax[0]);
    LamMax[0] = fmax(upa, LamMax[0]);
    elemCmax[e] = fmax(elemCmax[e], C);
    //printf("element %d: Cmax=%f\n",e, elemCmax[e]);
    
    //That's all, folks. dg_solver should now be holding proper CsMax,DivMax,LamMax values
  } //end of RAD Singlefluid case, 1D
#endif //end if for RAD Singlefluid physics
#endif //end if for 1D

#ifdef TWOD
#ifdef RADSINGLEFLUID
    {
      //get some field variables
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar mox   = Ug[(e*N_F+1)*N_G+g];   
      scalar moy   = Ug[(e*N_F+2)*N_G+g];   
      scalar u     = mox / rho;  // (rho u / rho) = u
      scalar v     = moy / rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar Cx    = Ug[(e*N_F+4)*N_G+g];
      scalar Cy    = Ug[(e*N_F+5)*N_G+g];
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u - 0.5*rho*v*v);

      //Gradients of conserved variables: subscript means derivative
      //I only need density and momentum here
      scalar rho_x = 0.0; scalar rho_y = 0.0;
      scalar mox_x = 0.0; scalar mox_y = 0.0;
      scalar moy_x = 0.0; scalar moy_y = 0.0;

      //10/23/2017: I've changed my opinion on how to use invJac
      for (int alpha = 0; alpha < D; alpha++)
	{
	  //both xi and eta derivates factor in to x and y derivatves in general.
	  //dUg_ref index alpha hits reference coordinates while 0,1 in invJac hit phys coords.
	  scalar dx = invJac[e*N_G*D*D + g*D*D + alpha*D + 0];
	  scalar dy = invJac[e*N_G*D*D + g*D*D + alpha*D + 1];
	  rho_x += dUg_ref[((e*N_F+0)*N_G+g)*D + alpha] * dx;
	  rho_y += dUg_ref[((e*N_F+0)*N_G+g)*D + alpha] * dy;

	  mox_x += dUg_ref[((e*N_F+1)*N_G+g)*D + alpha] * dx;
	  mox_y += dUg_ref[((e*N_F+1)*N_G+g)*D + alpha] * dy;

	  moy_x += dUg_ref[((e*N_F+2)*N_G+g)*D + alpha] * dx;
	  moy_y += dUg_ref[((e*N_F+2)*N_G+g)*D + alpha] * dy;
	}
      //Get velocity gradients for divergence calculation
      scalar Inv_RhoSq = 1.0/(rho*rho);
      scalar u_x = (rho*mox_x - mox*rho_x) * Inv_RhoSq;
      scalar v_y = (rho*moy_y - moy*rho_y) * Inv_RhoSq;

      //Calculate the local super-divergence (this is not the actual divergence)
      //scalar Div = fabs(u_x) + fabs(v_y) + 0.0;
      //Calculate the local wavespeeds
      scalar upa = fabs(u) + sqrt(gamma*p/rho);
      scalar vpa = fabs(v) + sqrt(gamma*p/rho);
      //Update global maxima
      DivMax[0] = fmax(DivMax[0], fabs(u_x));
      DivMax[1] = fmax(DivMax[1], fabs(v_y));
      
      CsMax[0] = fmax(Cx, CsMax[0]);
      CsMax[1] = fmax(Cy, CsMax[1]);

      LamMax[0] = fmax(upa, LamMax[0]);
      LamMax[1] = fmax(vpa, LamMax[1]);
      elemCmax[e] = fmax(elemCmax[e], Cx);
      elemCmax[e] = fmax(elemCmax[e], Cy);
      //printf("element %d: Cmax=%f\n",e, elemCmax[e]);
      //We have the global maxima, this subroutine's job is done :)
    }
#endif //end if for RADSINGLEFLUID physics
#endif //end if for 2D

#ifdef USE_CPU
    } //end g loop for CGU case
#endif
  } //end element loop

  } //END SUBROUTINE

  arch_global void GrabCmax_Ghosts(int N_G, int N_E, int Ne_AUG, scalar* Ug, scalar* elemCmax)
  {
    /*!
      \brief Calculates max C values in ghost elements
      \param[in] N_G gaussian quadrature nodes per element
      \param[in] N_E number of flesh elements in the partition
      \param[in] Ne_AUG number of flesh+ghost elements in the partition
      \param[in] Ug solution at gaussian quadrature nodes
      \param[in] dUg_ref gradient wrt reference coordinates at quadrature nodes
      \param[in] invJac inverse of element Jacobian at each quadrature node
      \param[out] elemCmax max C value per each element, no direction, just 1 per element
    */
    //This routine runs after GrabCparameters; elemCmax is mostly full,
    //I just need to grab the value in each ghost element
   
#ifdef USE_CPU
    for(int e = N_E; e < Ne_AUG; e++){ //like I said, just the ghost elements
      elemCmax[e] = 0.0; //initialize element's C variable to zero
      for(int g = 0; g < N_G; g++){
	
#elif USE_GPU
	ERROR
#endif
	  
	  //1D problem=======================================================
#ifdef ONED
#ifdef RADSINGLEFLUID
    {
    scalar C     = Ug[(e*N_F+3)*N_G+g];
    elemCmax[e] = fmax(elemCmax[e], C);
    //printf("element %d: Cmax=%f\n",e, elemCmax[e]);
    
  } //end of RAD Singlefluid case, 1D
#endif //end if for RAD Singlefluid physics
#endif //end if for 1D

#ifdef TWOD
#ifdef RADSINGLEFLUID
    {
      scalar Cx    = Ug[(e*N_F+4)*N_G+g];
      scalar Cy    = Ug[(e*N_F+5)*N_G+g];

      elemCmax[e] = fmax(elemCmax[e], Cx);
      elemCmax[e] = fmax(elemCmax[e], Cy);
      //printf("element %d: Cmax=%f\n",e, elemCmax[e]);
    }
#endif //end if for RADSINGLEFLUID physics
#endif //end if for 2D

#ifdef USE_CPU
    } //end g loop for CPU case
#endif
  } //end element loop

  } //END SUBROUTINE


//==========================================================================
arch_global void kinetic_energy1D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* K){
  /*!
    \brief Calculate the 1D kinetic energy.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] rho density
    \param[in] rhou momentum
    \param[out] K kinetic energy
    \section Description
    Can be tested with the following:
    scalar* a = new scalar[18];
    scalar* b = new scalar[18];
    scalar* c = new scalar[18];     
    for(int i=0; i<18; i++){a[i] = i+1;}
    for(int i=0; i<18; i++){b[i] = i+5;}
    scalar* d_a;
    scalar* d_b;
    scalar* d_c;
    cudaMalloc((void**) &d_a,18*sizeof(scalar));
    cudaMalloc((void**) &d_b,18*sizeof(scalar));
    cudaMalloc((void**) &d_c,18*sizeof(scalar));
    cudaMemcpy(d_a, a, 18*sizeof(scalar), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, 18*sizeof(scalar), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c, c, 18*sizeof(scalar), cudaMemcpyHostToDevice);
    //Lkinetic_energy1D(1,18,a,b,c);
    Lkinetic_energy1D(1,18,d_a,d_b,d_c);
    cudaMemcpy(c, d_c, 18*sizeof(scalar), cudaMemcpyDeviceToHost);
    for(int i=0; i<18; i++){printf("%i: %f =? %f\n",i,c[i],0.5*b[i]*b[i]/a[i]);}
    delete[] a;
    delete[] b;
    delete[] c;
    exit(0);
  */  
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif
    int idx = e*N_s+i;
    K[idx] = 0.5*rhou[idx]*rhou[idx]/rho[idx];
#ifdef USE_CPU
  }
#endif
  }
}

//==========================================================================
arch_global void kinetic_energy2D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* rhov, scalar* K){
  /*!
    \brief Calculate the 2D kinetic energy.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] rho density
    \param[in] rhou momentum
    \param[in] rhov momentum
    \param[out] K kinetic energy
  */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif
    int idx = e*N_s+i;
    K[idx] = 0.5*(rhou[idx]*rhou[idx]+rhov[idx]*rhov[idx])/rho[idx];
#ifdef USE_CPU
  }
#endif
  }
}

//==========================================================================
arch_global void pressure(int N_s, int N_E, scalar* U, scalar* p){
  /*!
    \brief Calculate the pressure.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution
    \param[out] p pressure
  */
  /*
    PEJ 05/29/2017: Modifying this subroutine for scalar advection-diffusion.
    In this case, pressure is irrelevant, so this function needs to
    do nothing. To this end, I'm putting the p calculation itself INSIDE
    all appropriate compiler-ifs.
   */
#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int i = 0; i < N_s; i++){
#elif USE_GPU
  int e = blockIdx.x*blkE+threadIdx.z;
  if (e < N_E){
    int i = threadIdx.x;
#endif

#ifdef ONED
      // Get the pressure field
      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      scalar E    = U[(e*N_F+2)*N_s+i];
#ifdef SINGLEFLUID
      scalar gamma = constants::GLOBAL_GAMMA;
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);
#elif MULTIFLUID
#ifdef GAMCONS
      scalar gamma=1.0+rho/U[(e*N_F+3)*N_s+i];
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);
#elif GAMNCON
      scalar gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);
#endif
#elif PASSIVE
      scalar gamma = constants::GLOBAL_GAMMA;
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);
#elif STIFFENED
      scalar gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
      scalar beta = U[(e*N_F+4)*N_s+i];
      E = E - beta; // correct the energy for stiffened EOS
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho); 
#endif //end if for stiffened physics
   

#elif TWOD   
   // Get the pressure field
      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      scalar rhov = U[(e*N_F+2)*N_s+i];
      scalar E    = U[(e*N_F+3)*N_s+i];
#ifdef SINGLEFLUID
      scalar gamma = constants::GLOBAL_GAMMA;
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#elif MULTIFLUID
#ifdef GAMCONS
      scalar gamma=1.0+rho/U[(e*N_F+4)*N_s+i];
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#elif GAMNCON
      scalar gamma=1.0+1.0/U[(e*N_F+4)*N_s+i];
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#endif //end gamma model if
#elif PASSIVE
      scalar gamma = constants::GLOBAL_GAMMA;
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#elif STIFFENED
      scalar gamma=1.0+1.0/U[(e*N_F+4)*N_s+i];
      scalar beta = U[(e*N_F+5)*N_s+i];
      E = E - beta; // correct the energy for stiffened EOS
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#endif //end if for stiffend physics

#elif THREED
      //get the pressure field
      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      scalar rhov = U[(e*N_F+2)*N_s+i];
      scalar rhow = U[(e*N_F+3)*N_s+i];
      scalar E    = U[(e*N_F+4)*N_s+i];
#ifdef SINGLEFLUID
      scalar gamma = constants::GLOBAL_GAMMA;
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow) / rho);
#endif //end if on physics

#endif //end if on dimensions
      
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
void Levaluate_sf(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){//, scalar* xyz){
  /*!
    \brief Host C function to lauch evaluate_sf kernel.
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_E number of elements
    \param[out] s source array
    \param[out] f element flux array
    \param[in] Ug solution (evaluated at gaussian nodes)
    \param[in] dUg derivative of solution (evaluated at gaussian nodes)
    \param[in] invJac inverse Jacobian of the elements
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

  evaluate_sf arch_args (N_G, N_E, s, f, Ug, dUg, invJac, constants::GLOBAL_GX, constants::GLOBAL_GY);//, xyz);
}


extern "C" 
void Levaluate_q(int M_G, int M_T, scalar* q, scalar* UgF, scalar* normals){//, scalar* xyzf){
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

  evaluate_q arch_args_array(blkT*M_G*2*N_F*sizeof(scalar)) (M_G, M_T, q, UgF, normals);//, xyzf);
}

extern "C"
  void Levaluate_sf_rad(int N_G, int N_E, int order, scalar* DivMax, scalar Beta_S, scalar Mew_S, scalar* CsMax, scalar* LamMax, scalar eps_gen, scalar* ADeps, scalar Cthresh, scalar* elemCmax, int* sensor, scalar* invJac, scalar* Ug, scalar* dUg_phys, scalar* s, scalar* f)
{
   /*!
      \brief host C function to launch evaluate_sf_rad
      \param[in] N_G number of gaussian nodes per element
      \param[in] N_E number of elements
      \param[in] order DG solution order
      \param[in] DivMax maximum divergence (del dot velocity) in the entire domain
      \param[in] Beta_S global free paramter for setting strength of AD
      \param[in] CsMax The maximum value of C anywhere in the domain
      \param[in] LamMax maximum wavespeed anywhere in the domain (u+a)
      \param[in] eps_gen small value to go in denominators, prevent division by zero
      \param[in] ADeps each element's (h/p)^2 value to scalar artificial dissipation
      \param[in] invJac inverse jacobian of the elements
      
      \param[in] Ug solution (evaluated at gaussian nodes)
      \param[in] dUg_phys gradient of solution wrt physical coordinates (evaluated at gaussian nodes)
      \param[in and out] s element source term array
      \param[in and out] f element flux array, incremented by viscous flux
    */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  evaluate_sf_rad arch_args (N_G, N_E, order, DivMax, Beta_S, Mew_S, CsMax, LamMax, eps_gen, ADeps, Cthresh, elemCmax, sensor, invJac, Ug, dUg_phys, s, f);
}

extern "C" 
  void Levaluate_q_rad(int M_G, int M_T, scalar* q, scalar* UgF, scalar* gradCommon, scalar* DivMax, scalar Beta_S, scalar Mew_S, scalar* CsMax, scalar* LamMax, scalar eps_gen, scalar* ADepsF, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* normals){//, scalar* xyzf){
  /*!
    \brief Host C function to lauch evaluate_q_rad kernel.
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[in and out] q interface flux array
    \param[in] UgF interface solution (evaluated at gaussian nodes)
    \param[in] gradCommon interface gradient wrt physical coordinates
    \param[in] DivMax max divergence in domain
    \param[in] Beta_S global free paramter for setting strength of AD
    \param[in] CsMax The maximum value of C anywhere in the domain
    \param[in] LamMax maximum wavespeed anywhere in the domain (u+a)
    \param[in] eps_gen small value to go in denominators, prevent division by zero
    \param[in] ADepsF Face average of each element's (h/p)^2 
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

  evaluate_q_rad arch_args_array(blkT*M_G*2*N_F*sizeof(scalar)) (M_G, M_T, q, UgF, gradCommon, DivMax, Beta_S, Mew_S, CsMax, LamMax, eps_gen, ADepsF, Cthresh, elemCmax, sensor, BR2_Map, normals);//, xyzf);
}

extern "C"
void LGrabCparameters(int N_G, int N_E, scalar* Ug, scalar* dUg_ref, scalar* invJac, scalar* DivMax, scalar* LamMax, scalar* CsMax, scalar* elemCmax)
{
  /*!
     \brief Host C function to GrabCparameters
     \param[in] N_G gaussian quadrature nodes per element
     \param[in] N_E number of elements in the domain
     \param[in] Ug solution at gaussian quadrature nodes
     \param[in] dUg_ref gradient wrt reference coordinates at quadrature nodes
     \param[in] invJac inverse of element Jacobian at each quadrature node
     \param[out] DivMax max divergence in entire domain
     \param[out] LamMax max wavespeed in entire domain
     \param[out] CsMax max Cs value in entire domain
  */
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(N_G,1,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  GrabCparameters arch_args (N_G, N_E, Ug, dUg_ref, invJac, DivMax, LamMax, CsMax, elemCmax);
}

extern "C"
  void LGrabCmax_Ghosts(int N_G, int N_E, int Ne_AUG, scalar* Ug, scalar* elemCmax)
{
#ifdef USE_GPU
  printf("GrabCParameters_CsGhosts not ready for GPU\n");
  exit(1);
#endif
#ifdef USE_CPU
  GrabCmax_Ghosts arch_args(N_G, N_E, Ne_AUG, Ug, elemCmax);
#endif
}

extern "C"
void Lkinetic_energy1D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* K){
  /*!
    \brief Host C function to lauch kinetic_energy1D kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] rho density
    \param[in] rhou momentum
    \param[out] K kinetic energy
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

  kinetic_energy1D arch_args (N_s, N_E, rho, rhou, K);
}

extern "C"
void Lkinetic_energy2D(int N_s, int N_E, scalar* rho, scalar* rhou, scalar* rhov, scalar* K){
  /*!
    \brief Host C function to lauch kinetic_energy1D kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] rho density
    \param[in] rhou momentum
    \param[in] rhov momentum
    \param[out] K kinetic energy
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

  kinetic_energy2D arch_args (N_s, N_E, rho, rhou, rhov, K);
}

extern "C"
void Lpressure(int N_s, int N_E, scalar* U, scalar* p){
  /*!
    \brief Host C function to lauch pressure kernel.
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] U solution
    \param[out] p pressure
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

  pressure arch_args (N_s, N_E, U, p);
}
