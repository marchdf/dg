/*!
  \file physics.cu
  \brief Kernels used for the physics
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "physics.h"
#include "basic_fluxes.h"
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
#ifdef PASSIVE //===========================================================
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

#elif SINGLEFLUID //=========================================================
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
      
#elif MULTIFLUID //=========================================================
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

#elif STIFFENED //=========================================================
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

#endif // end ifs on physics


      //===================================================================
      //
      // 2D problem
      //
      //===================================================================
#elif TWOD
#ifdef PASSIVE //==========================================================
      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar phic  = Ug[(e*N_F+4)*N_G+g]/rho; 
      scalar vdotgradphinc = 0; // = -u*dphincdx-v*dphincdy
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradphinc += 
	  -u*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
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

#elif SINGLEFLUID //========================================================
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
      
#elif MULTIFLUID //========================================================
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

#elif STIFFENED //========================================================
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

#endif // end ifs on physics

#endif // end ifs on dimensions
     
#ifdef USE_CPU
    }
#endif
  }
}


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

     
      // Send the data to the Riemann solvers
#ifdef ONED

#ifdef PASSIVE //=========================================================
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

#elif SINGLEFLUID //=========================================================
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
#endif // flux if

#elif MULTIFLUID //=========================================================

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

#elif STIFFENED //=========================================================

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

#endif // physics if

#elif TWOD

#ifdef PASSIVE

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

#elif SINGLEFLUID //=========================================================

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
#endif // flux if

#elif MULTIFLUID //=========================================================

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

#elif STIFFENED //=========================================================

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

#endif // physics if

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
#elif MULTIFLUID
#ifdef GAMCONS
      scalar gamma=1.0+rho/U[(e*N_F+3)*N_s+i];
#elif GAMNCON
      scalar gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
#endif
#elif PASSIVE
      scalar gamma = constants::GLOBAL_GAMMA;
#elif STIFFENED
      scalar gamma=1.0+1.0/U[(e*N_F+3)*N_s+i];
      scalar beta = U[(e*N_F+4)*N_s+i];
      E = E - beta; // correct the energy for stiffened EOS
#endif
      p[e*N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);

#elif TWOD
      // Get the pressure field
      scalar rho  = U[(e*N_F+0)*N_s+i];
      scalar rhou = U[(e*N_F+1)*N_s+i];
      scalar rhov = U[(e*N_F+2)*N_s+i];
      scalar E    = U[(e*N_F+3)*N_s+i];
#ifdef SINGLEFLUID
      scalar gamma = constants::GLOBAL_GAMMA;
#elif MULTIFLUID
#ifdef GAMCONS
      scalar gamma=1.0+rho/U[(e*N_F+4)*N_s+i];
#elif GAMNCON
      scalar gamma=1.0+1.0/U[(e*N_F+4)*N_s+i];
#endif
#elif PASSIVE
      scalar gamma = constants::GLOBAL_GAMMA;
#elif STIFFENED
      scalar gamma=1.0+1.0/U[(e*N_F+4)*N_s+i];
      scalar beta = U[(e*N_F+5)*N_s+i];
      E = E - beta; // correct the energy for stiffened EOS
#endif
      p[e*N_s+i] = (gamma-1)*(E - 0.5*(rhou*rhou+rhov*rhov)/rho);
#endif
      
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


// //
// // Possibly broken
// // 
// //==========================================================================
// arch_global void evaluate_sf_shallow(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

// #ifdef USE_CPU
//   for(int e = 0; e < N_E; e++){
//     for(int g = 0; g < N_G; g++){
// #elif USE_GPU
//       int e = blockIdx.x;
//       int g = threadIdx.x;
// #endif
  
//       scalar eta  =  Ug[(e*N_F+0)*N_G+g];
  
//       s[(e*N_F+0)*N_G+g] = 0;
//       s[(e*N_F+1)*N_G+g] = 0;
//       s[(e*N_F+2)*N_G+g] = 0;
  
//       // Flux derive par rapport a x
//       f[((e*N_F+0)*N_G+g)*D+0] = H0*Ug[(e*N_F+1)*N_G+g]; // u_x
//       f[((e*N_F+1)*N_G+g)*D+0] = G0*eta; // eta 
//       f[((e*N_F+2)*N_G+g)*D+0] = 0; 
  
//       // Flux derive par rapport a y
//       f[((e*N_F+0)*N_G+g)*D+1] = H0*Ug[(e*N_F+2)*N_G+g]; // u_y
//       f[((e*N_F+1)*N_G+g)*D+1] = 0;
//       f[((e*N_F+2)*N_G+g)*D+1] = G0*eta; // eta

// #ifdef USE_CPU
//     }
//   }
// #endif
// }

// arch_device scalar flux1_mhd(scalar rho, scalar u){return rho*u;}                                        // for f0X, f0X, f0X
// arch_device scalar flux2_mhd(scalar rho, scalar u, scalar Bx, scalar pbar){return rho*u*u-Bx*Bx+pbar;}   // for f1X, f2Y, f3Z
// arch_device scalar flux3_mhd(scalar rho, scalar u, scalar v, scalar Bx, scalar By){return rho*u*v-Bx*By;}// for f1Y, f1Z, f2X, f2Z, f3X, f3Y
// arch_device scalar flux4_mhd(scalar u, scalar By, scalar v, scalar Bx){return u*By-v*Bx;}                // for f4Y, f4Z, f5X, f5Z, f6X, f6Y
// arch_device scalar flux5_mhd(scalar EtplusPbar, scalar u, scalar vdotB, scalar Bx) {return EtplusPbar*u - vdotB*Bx;} // for f7X, f7Y, f7Z

// //==========================================================================
// arch_global void evaluate_sf_mhd(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

// #ifdef USE_CPU
//   for(int e = 0; e < N_E; e++){
//     for(int g = 0; g < N_G; g++){
// #elif USE_GPU
//       int e = blockIdx.x;
//       int g = threadIdx.x;
// #endif

//       scalar rho = Ug[(e*N_F+0)*N_G+g];   
//       scalar u   = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
//       scalar v   = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
//       scalar Bx  = Ug[(e*N_F+3)*N_G+g];
//       scalar By  = Ug[(e*N_F+4)*N_G+g];
//       scalar Et  = Ug[(e*N_F+5)*N_G+g];
//       scalar vdotv  = u*u+v*v;       
//       scalar BdotB  = Bx*Bx+By*By;
//       scalar vdotB  = u*Bx+v*By;
//       scalar divB   = 0;
//       for(int alpha = 0; alpha < D; alpha++){
// 	divB += dUg[((e*N_F+3)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] + dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];
//       }
//       scalar p = (gamma-1)*(Et - 0.5*(rho*vdotv+BdotB));
//       scalar Pbar = p+0.5*BdotB;
//       scalar EtPbar = Et+Pbar;

//       // Powel source
//       s[(e*N_F+0)*N_G+g] = 0;
//       s[(e*N_F+1)*N_G+g] = 0;//-divB*Bx;
//       s[(e*N_F+2)*N_G+g] = 0;//-divB*By;
//       s[(e*N_F+3)*N_G+g] = 0;//-divB*u;
//       s[(e*N_F+4)*N_G+g] = 0;//-divB*v;
//       s[(e*N_F+5)*N_G+g] = 0;//-divB*vdotB;
            
//       // Flux derive par rapport a x
//       f[((e*N_F+0)*N_G+g)*D+0] = flux1_mhd(rho,u);                //rho*u; 
//       f[((e*N_F+1)*N_G+g)*D+0] = flux2_mhd(rho,u,Bx,Pbar);        //rho*u*u-Bx*Bx+Pbar; 
//       f[((e*N_F+2)*N_G+g)*D+0] = flux3_mhd(rho,u,v,Bx,By);        //rho*u*v-Bx*By; 
//       f[((e*N_F+3)*N_G+g)*D+0] = 0;                                   //0;
//       f[((e*N_F+4)*N_G+g)*D+0] = flux4_mhd(u,By,v,Bx);            //u*By-v*Bx;
//       f[((e*N_F+5)*N_G+g)*D+0] = flux5_mhd(EtPbar,u,vdotB,Bx);    //EtplusPbar*u-vdotB*Bx;
      
//       // Flux derive par rapport a y
//       f[((e*N_F+0)*N_G+g)*D+1] = flux1_mhd(rho,v);                //rho*v;
//       f[((e*N_F+1)*N_G+g)*D+1] = flux3_mhd(rho,v,u,By,Bx);        //rho*v*u-By*Bx;
//       f[((e*N_F+2)*N_G+g)*D+1] = flux2_mhd(rho,v,By,Pbar);        //rho*v*v-By*By+Pbar;
//       f[((e*N_F+3)*N_G+g)*D+1] = flux4_mhd(v,Bx,u,By);            //v*Bx-u*By;
//       f[((e*N_F+4)*N_G+g)*D+1] = 0;                                   //0;
//       f[((e*N_F+5)*N_G+g)*D+1] = flux5_mhd(EtPbar,v,vdotB,By);    //EtplusPbar*v-vdotB*By;

// #ifdef USE_CPU
//     }
//   }
// #endif
// }

// //==========================================================================
// arch_global void evaluate_q_shallow(int M_G, int M_T, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

// #ifdef USE_CPU
//   for(int t = 0; t < M_T; t++){
//     for(int g = 0; g < M_G; g++){
// #elif USE_GPU
//       int t = blockIdx.x;
//       int g = threadIdx.x;
// #endif

//       scalar nx = normals[t*2+0];
//       scalar ny = normals[t*2+1];
//       scalar etaL= UgF[((t*N_F+0)*2+0)*M_G+g];
//       scalar etaR= UgF[((t*N_F+0)*2+1)*M_G+g];
//       scalar uLn = UgF[((t*N_F+1)*2+0)*M_G+g] * nx + UgF[((t*N_F+2)*2+0)*M_G+g] * ny;
//       scalar uRn = UgF[((t*N_F+1)*2+1)*M_G+g] * nx + UgF[((t*N_F+2)*2+1)*M_G+g] * ny;
      
//       scalar h0 = H0;
//       scalar g0 = G0;
      
//       // first equation
//       scalar qL = -0.5*h0*(uLn + uRn + sqrt(g0/h0)*(etaL-etaR)); // Left
//       q[((t*N_F+0)*2+0)*M_G+g] = qL;
//       q[((t*N_F+0)*2+1)*M_G+g] = -qL;
//       // second
//       qL = -0.5*g0*nx*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
//       q[((t*N_F+1)*2+0)*M_G+g] = qL;
//       q[((t*N_F+1)*2+1)*M_G+g] = -qL;
//       // third
//       qL = -0.5*g0*ny*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
//       q[((t*N_F+2)*2+0)*M_G+g] = qL;
//       q[((t*N_F+2)*2+1)*M_G+g] = -qL;

// #ifdef USE_CPU
//     }
//   }
// #endif
// }


// //==========================================================================
// arch_global void evaluate_q_mhd(int M_G, int M_T, scalar* q, scalar* UgF, scalar gamma, scalar* normals){
  
// #ifdef USE_CPU
//   for(int t = 0; t < M_T; t++){
//     scalar* vap = new scalar[M_G*2*8];
//     for(int g = 0; g < M_G; g++){
// #elif USE_GPU
//       int t = blockIdx.x;
//       int g = threadIdx.x;
//       extern __shared__ scalar vap[];
// #endif

//       scalar nx = normals[t*2+0];
//       scalar ny = normals[t*2+1];
//       scalar rhoL= UgF[((t*N_F+0)*2+0)*M_G+g];
//       scalar rhoR= UgF[((t*N_F+0)*2+1)*M_G+g];
//       scalar uL  = UgF[((t*N_F+1)*2+0)*M_G+g]/rhoL;
//       scalar uR  = UgF[((t*N_F+1)*2+1)*M_G+g]/rhoR;
//       scalar vL  = UgF[((t*N_F+2)*2+0)*M_G+g]/rhoL;
//       scalar vR  = UgF[((t*N_F+2)*2+1)*M_G+g]/rhoR;
//       scalar BxL = UgF[((t*N_F+3)*2+0)*M_G+g];
//       scalar BxR = UgF[((t*N_F+3)*2+1)*M_G+g];
//       scalar ByL = UgF[((t*N_F+4)*2+0)*M_G+g];
//       scalar ByR = UgF[((t*N_F+4)*2+1)*M_G+g];
//       scalar EtL = UgF[((t*N_F+5)*2+0)*M_G+g];
//       scalar EtR = UgF[((t*N_F+5)*2+1)*M_G+g];
//       scalar vdotvL = uL*uL+vL*vL;
//       scalar vdotvR = uR*uR+vR*vR;
//       scalar BdotBL = BxL*BxL+ByL*ByL;
//       scalar BdotBR = BxR*BxR+ByR*ByR;
//       scalar vdotBL = uL*BxL+vL*ByL;
//       scalar vdotBR = uR*BxR+vR*ByR;
//       scalar vdotnL = uL*nx+vL*ny;
//       scalar vdotnR = uR*nx+vR*ny;
//       scalar BdotnL = BxL*nx+ByL*ny;
//       scalar BdotnR = BxR*nx+ByR*ny;
//       //scalar aveBdotn = 0.5*(BdotnL + BdotnR);      
//       scalar pL = (gamma-1)*(EtL - 0.5*(rhoL*vdotvL+BdotBL));
//       scalar pR = (gamma-1)*(EtR - 0.5*(rhoR*vdotvR+BdotBR));
//       scalar PbarL = pL+0.5*BdotBL;
//       scalar PbarR = pR+0.5*BdotBR;
//       scalar EtPbarL = EtL+PbarL;
//       scalar EtPbarR = EtR+PbarR;

//       // Evaluate the right and left eigenvalues
//       int sizevap = 8;
//       scalar alfenL  = BdotnL/sqrt(rhoL);
//       scalar alfenR  = BdotnR/sqrt(rhoR);
//       scalar a2L = (gamma*pL+BdotBL)/rhoL;
//       scalar a2R = (gamma*pR+BdotBR)/rhoR;
//       scalar cfL = sqrt(0.5*(a2L + sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
//       scalar cfR = sqrt(0.5*(a2R + sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));
//       scalar csL = sqrt(0.5*(a2L - sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
//       scalar csR = sqrt(0.5*(a2R - sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));

//       vap[(g*2+0)*sizevap+0] = fabs(vdotnL);
//       vap[(g*2+0)*sizevap+1] = fabs(vdotnL) + alfenL;
//       vap[(g*2+0)*sizevap+2] = fabs(vdotnL) - alfenL;
//       vap[(g*2+0)*sizevap+3] = fabs(vdotnL) + cfL;
//       vap[(g*2+0)*sizevap+4] = fabs(vdotnL) - cfL;
//       vap[(g*2+0)*sizevap+5] = fabs(vdotnL) + csL;
//       vap[(g*2+0)*sizevap+6] = fabs(vdotnL) - csL;
//       vap[(g*2+0)*sizevap+7] = fabs(vdotnL);
      
//       vap[(g*2+1)*sizevap+0] = fabs(vdotnR);
//       vap[(g*2+1)*sizevap+1] = fabs(vdotnR) + alfenR;
//       vap[(g*2+1)*sizevap+2] = fabs(vdotnR) - alfenR;
//       vap[(g*2+1)*sizevap+3] = fabs(vdotnR) + cfR;
//       vap[(g*2+1)*sizevap+4] = fabs(vdotnR) - cfR;
//       vap[(g*2+1)*sizevap+5] = fabs(vdotnR) + csR;
//       vap[(g*2+1)*sizevap+6] = fabs(vdotnR) - csR;
//       vap[(g*2+1)*sizevap+7] = fabs(vdotnR);

//       scalar maxvap = 0;
//       for (int k = 0; k < 2*sizevap; k++){
// 	if (maxvap<vap[g*16+k]) maxvap = vap[g*16+k];
//       }

//       // Upwinding on the source term
//       // scalar upBdotnL = aveBdotn;
//       // scalar upBdotnR = aveBdotn;
//       scalar upBdotnL = 0.0;
//       scalar upBdotnR = 0.0;
//       if      (vdotnL >= 0) upBdotnL = BdotnL;
//       else if (vdotnL <  0) upBdotnL = BdotnR;
//       if      (vdotnR >= 0) upBdotnR = BdotnL;
//       else if (vdotnR <  0) upBdotnR = BdotnR;  
      
//       //
//       // Evaluate the fluxes on the right and left
//       //

//       //first: fx = rho*u; fy = rho*v; fz = rho*w; 
//       scalar qL = -0.5*((flux1_mhd(rhoL,uL) + flux1_mhd(rhoR,uR))*nx +
// 			(flux1_mhd(rhoL,vL) + flux1_mhd(rhoR,vR))*ny 
// 			-maxvap*(rhoR-rhoL));
//       q[((t*N_F+0)*2+0)*M_G+g] = qL;
//       q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      
//       //second: fx = rho*u*u+Bx*Bx+Pbar; fy = rho*v*u-By*Bx; fz = rho*w*u-Bz*Bx;
//       qL = -0.5*((flux2_mhd(rhoL,uL,BxL,PbarL)  + flux2_mhd(rhoR,uR,BxR,PbarR) )*nx +
// 		 (flux3_mhd(rhoL,vL,uL,ByL,BxL) + flux3_mhd(rhoR,vR,uR,ByR,BxR))*ny
// 		 -maxvap*(rhoR*uR-rhoL*uL));
//       q[((t*N_F+1)*2+0)*M_G+g] = qL  - BxL*upBdotnL;
//       q[((t*N_F+1)*2+1)*M_G+g] = -qL + BxR*upBdotnR;

//       //third: fx = rho*u*v-Bx*By; fy = rho*v*v-By*By+Pbar; fz = rho*w*v-Bz*By;
//       qL = -0.5*((flux3_mhd(rhoL,uL,vL,BxL,ByL) + flux3_mhd(rhoR,uR,vR,BxR,ByR))*nx +
// 		 (flux2_mhd(rhoL,vL,ByL,PbarL)  + flux2_mhd(rhoR,vR,ByR,PbarR) )*ny
// 		 -maxvap*(rhoR*vR-rhoL*vL));
//       q[((t*N_F+2)*2+0)*M_G+g] = qL  - ByL*upBdotnL;
//       q[((t*N_F+2)*2+1)*M_G+g] = -qL + ByR*upBdotnR;
      
//       //fourth: fx = 0;  fy = v*Bx-u*By; fz = w*Bx-u*Bz;
//       qL = -0.5*((0                            + 0                           )*nx +
// 		 (flux4_mhd(vL,BxL,uL,ByL) + flux4_mhd(vR,BxR,uR,ByR))*ny
// 		 -maxvap*(BxR-BxL));
//       q[((t*N_F+3)*2+0)*M_G+g] = qL  - uL*upBdotnL;
//       q[((t*N_F+3)*2+1)*M_G+g] = -qL + uR*upBdotnR;

//       //fifth: fx = u*By-v*Bx;  fy = 0; fz = w*By-v*Bz;
//       qL = -0.5*((flux4_mhd(uL,ByL,vL,BxL) + flux4_mhd(uR,ByR,vR,BxR))*nx +
// 		 (0                            + 0                           )*ny 
// 		 -maxvap*(ByR-ByL));
//       q[((t*N_F+4)*2+0)*M_G+g] = qL  - vL*upBdotnL;
//       q[((t*N_F+4)*2+1)*M_G+g] = -qL + vR*upBdotnR;

//       //sixth: fx = EtplusPbar*u-vdotB*Bx; fy = EtplusPbar*v-vdotB*By; fz = EtplusPbar*w-vdotB*Bz;
//       qL = -0.5*((flux5_mhd(EtPbarL,uL,vdotBL,BxL) + flux5_mhd(EtPbarR,uR,vdotBR,BxR))*nx +
// 		 (flux5_mhd(EtPbarL,vL,vdotBL,ByL) + flux5_mhd(EtPbarR,vR,vdotBR,ByR))*ny 
// 		 -maxvap*(EtR-EtL));
//       q[((t*N_F+5)*2+0)*M_G+g] = qL  - vdotBL*upBdotnL;
//       q[((t*N_F+5)*2+1)*M_G+g] = -qL + vdotBR*upBdotnR;

// #ifdef USE_CPU
//     }
//     delete[] vap;
//   }
// #endif
// }

// extern "C" 
// void Levaluate_sf_shallow(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

// #ifdef USE_GPU
//   dim3 dimBlock(N_G,1,1);
//   dim3 dimGrid(N_E,1);
// #endif

//   evaluate_sf_shallow arch_args (N_G, N_E, s, f, Ug, H0, G0);
// }

// extern "C" 
// void Levaluate_sf_mhd(int N_G, int N_E, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

// #ifdef USE_GPU
//   dim3 dimBlock(N_G,1,1);
//   dim3 dimGrid(N_E,1);
// #endif

//   evaluate_sf_mhd arch_args (N_G, N_E, s, f, Ug, dUg, invJac, gamma);
// }

// extern "C" 
// void Levaluate_q_shallow(int M_G, int M_T, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

// #ifdef USE_GPU
//   dim3 dimBlock(M_G,1,1);
//   dim3 dimGrid(M_T,1);
// #endif

//   evaluate_q_shallow arch_args (M_G, M_T, q, UgF, H0, G0, normals);
// }

// extern "C" 
// void Levaluate_q_mhd(int M_G, int M_T, scalar* q, scalar* UgF, scalar gamma, scalar* normals){

// #ifdef USE_GPU
//   dim3 dimBlock(M_G,1,1);
//   dim3 dimGrid(M_T,1);
// #endif

//   evaluate_q_mhd arch_args_array(M_G*2*8*sizeof(scalar)) (M_G, M_T, q, UgF, gamma, normals);
// }
