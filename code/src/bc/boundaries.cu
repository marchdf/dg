/*!
  \file boundaries.cu
  \brief Kernels to implement special boundary conditions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup boundaries
*/
#include "boundaries.h"
#include <cstdlib>
#include <stdio.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void rflctive(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement reflective boundary conditions
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of reflective BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED
      // velocity flip
      UF[((t*N_F+1)*2+1)*M_s+j] = -UF[((t*N_F+1)*2+1)*M_s+j];
#elif TWOD
      // normal
      scalar nx = normals[t*D+0]; 
      scalar ny = normals[t*D+1]; 
      scalar nx2 = nx*nx;
      scalar ny2 = ny*ny;
      scalar m2nxny = -2*nx*ny;
      scalar invnx2ny2 = 1.0/(nx2+ny2);
      
      // Normal and tangential velocities
      scalar vxL = UF[((t*N_F+1)*2+0)*M_s+j];
      scalar vyL = UF[((t*N_F+2)*2+0)*M_s+j];
      scalar vxR = invnx2ny2 * ((-nx2+ny2)*vxL +    m2nxny*vyL);
      scalar vyR = invnx2ny2 * (   m2nxny *vxL + (nx2-ny2)*vyL);
      
      // velocities 
      UF[((t*N_F+1)*2+1)*M_s+j] = vxR;
      UF[((t*N_F+2)*2+1)*M_s+j] = vyR;

#endif //on dimensions
      
#ifdef USE_CPU
    }
  }
#endif
}
//========================================================================================
arch_global void noslip(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement noslip boundary conditions
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of reflective BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED
      // velocity flip (same as reflective boundary)
      UF[((t*N_F+1)*2+1)*M_s+j] = -UF[((t*N_F+1)*2+1)*M_s+j];
#elif TWOD
      //For the 2D case, the no-slip boundary flips both normal and 
      //tangential velocity
      //In the future, this function should be modified to
      //support a tangential velocity, like driven cavity.
      //For that, you will need to make use of normals,
      //similar to rflctive approach.

      //It's not important here, but note in future applications
      //that this procedure flips the momentum vectors
      //directly, not the velocity

      scalar vxL = UF[((t*N_F+1)*2+0)*M_s+j];
      scalar vyL = UF[((t*N_F+2)*2+0)*M_s+j];
      
      //velocities: For the ghost state, reverse both velocity components 
      UF[((t*N_F+1)*2+1)*M_s+j] = -vxL;
      UF[((t*N_F+2)*2+1)*M_s+j] = -vyL;

#endif //on dimensions
      
#ifdef USE_CPU
    }
    }
#endif
}
//==========================================================================
//Note the nograd routine works directly on quadrature points, while other routines
//here deal with DOF
arch_global void nograd(int M_G, int M_B, int* boundaryMap, scalar* normals, int start, scalar* dUgF){
  /*!
    \brief Implement zero-gradient boundary conditions, specifically for viscous flux calculation
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of zero-gradient BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
  //printf("Entered the nograd boundary treatment\n");
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    //printf("\tk = %d, t = %d:\n", k,t);
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      //int j = threadIdx.x; //need to adjust this argument for dealing with g instead. Ask Marc about it
      int g = threadIdx.x; //need to adjust this argument for dealing with g instead. Ask Marc about it
#endif

#ifdef ONED
      //printf("\tEntered 1D loop, t = %d\n",t);
      //set the gradient to zero on BOTH sides of interface
      //Df_Dalpha_S(g) = dUgF(t, f, S, alpha, g)
      //Left side of interface:
      dUgF[t*N_F*2*D*M_G + 0*2*D*M_G + 0*D*M_G + 0*M_G + g] = 0.0;
      dUgF[t*N_F*2*D*M_G + 1*2*D*M_G + 0*D*M_G + 0*M_G + g] = 0.0;
      dUgF[t*N_F*2*D*M_G + 2*2*D*M_G + 0*D*M_G + 0*M_G + g] = 0.0;
      //right side of interface:
      dUgF[t*N_F*2*D*M_G + 0*2*D*M_G + 1*D*M_G + 0*M_G + g] = 0.0;
      dUgF[t*N_F*2*D*M_G + 1*2*D*M_G + 1*D*M_G + 0*M_G + g] = 0.0;
      dUgF[t*N_F*2*D*M_G + 2*2*D*M_G + 1*D*M_G + 0*M_G + g] = 0.0;
      
#elif TWOD
      //set the gradient to zero on BOTH sides of interface
      //Df_Dalpha_S(g) = dUgF(t, f, S, alpha, g)

      //This part is tricky: preserve the tangential gradients, set
      //normal gradients to zero.
      //Just like the normal, dUgF is already organized wrt physical coordintates

      // normal
      //printf("\tEntered 2D loop, t = %d\n",t);
      scalar nx = normals[t*D+0]; 
      scalar ny = normals[t*D+1]; 
      //printf("\tnx, ny = {%f,%f}\n",nx,ny);
      /*
      scalar nx2 = nx*nx;
      scalar ny2 = ny*ny;
      */
      /*
      scalar m2nxny = -2*nx*ny;
      scalar invnx2ny2 = 1.0/(nx2+ny2);
      */
      scalar dUdx_in[N_F];
      scalar dUdy_in[N_F];
      scalar dUdx_out[N_F];
      scalar dUdy_out[N_F];
      for (int f = 0; f < N_F; f++)
	{
          //take the input only from 0 side of interface, it's the interior state, while 1 is ghost state
	  //Df_Dalpha_S(g) = dUgF(t, f, S, alpha, g)
          dUdx_in[f] = dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 0*D*M_G + 0*M_G + g]; //0 for dUdx
          dUdy_in[f] = dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 0*D*M_G + 1*M_G + g]; //1 for dUdy
	  //printf("\t\tfield = %d, dUdx = %f, dUdy = %f\n",f,dUdx_in[f],dUdy_in[f]);
        }
      //The x,y gradient has now been imported: now adjust it
      //The idea here is dUgF(dot)normal = 0, dUgF(dot)tangent = untouched
      //Note that for a reflective boundary, we set the normal component negative, so this is slightly different
      //For some details on the math here, see Phil's notebook, 01/20/2016 entry
      //printf("\tNow applying treatment\n");
      for (int f = 0; f < N_F; f++)
	{
	  //dUdx_out[f] = 0.0;
	  //dUdy_out[f] = dUdy_in[f];
	  //printf("\t\tfield = %d, dUdx = %f, dUdy = %f\n",f,dUdx_out[f],dUdy_out[f]);
          dUdx_out[f] = ny*(dUdx_in[f]*ny - dUdy_in[f]*nx);
          dUdy_out[f] = nx*(dUdy_in[f]*nx - dUdx_in[f]*ny);
        }
      //Now, relay the boundary-treated gradient to the dUgF array, both sides of interface
      for (int f = 0; f < N_F; f++)
	{
          dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 0*D*M_G + 0*M_G + g] = dUdx_out[f]; //side 0, x derivative
          dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 1*D*M_G + 0*M_G + g] = dUdx_out[f]; //side 1, x derivative
          dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 0*D*M_G + 1*M_G + g] = dUdy_out[f]; //side 0, y derivative
          dUgF[t*N_F*2*D*M_G + f*2*D*M_G + 1*D*M_G + 1*M_G + g] = dUdy_out[f]; //side 1, y derivative
        }

#endif //on dimensions
      
#ifdef USE_CPU
    }
  }
#endif
}

scalar getMew(scalar mew_ref, scalar T_ref, scalar Cvis, scalar T)
{
  return mew_ref * (T_ref+Cvis) / (T+Cvis) * pow ((T+0.0)/(T_ref+0.0),1.5);
}
//========================================================================================
arch_global void Anflw(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement A inflow boundary conditions
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of A inflow BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED
     
      //There is no imediate 1D application for this condition
      scalar rho = 1.0;
      scalar p = 1.0;
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar sos = sqrt(gamma*p/rho);
      scalar u = 2.15*sos; //Mach 2.15 inflow
      scalar Et = p/(gamma-1.0) + 0.5*rho*u*u;
      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = Et; //side 1, field 2

#elif TWOD
      /*
      //These conditions are based on inflow specifications
      //for HiOCFD4 workshop, problem BL2
      //Everything here can be considered a freestream variable
      scalar T = 288.15;
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar Rgas  = constants::GLOBAL_RGAS;
      scalar Mach = 2.15;
      scalar sos = sqrt(gamma*Rgas*T);
      scalar u = Mach*sos;
      scalar v = 0;
      //Use Sutherland to get freestream viscosity:
      scalar mew_ref = constants:: GLOBAL_MEWREF;
      scalar T_ref = constants::GLOBAL_TREF;
      scalar Cvis = constants::GLOBAL_CVIS;

      scalar mew = mew_ref * (T_ref+Cvis) / (T+Cvis) * pow(T/T_ref,1.5);

      //Establish shock position:
      scalar Xshock = 1.0;
      scalar Re = powf(10,5); //Reynolds number
      scalar rho = (Re*mew)/(u*Xshock);
      scalar p = rho*Rgas*T;
      scalar Et = p/(gamma-1) + 0.5*rho*(u*u + v*v);
      //printf("Boundary Values {mew, rho, u, v, p} = {%f, %f, %f, %f, %f}\n", mew, rho,u,v,p);

      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = rho*v; //side 1, field 2
      UF[((t*N_F+3)*2+1)*M_s+j] = Et; //side 1, field 3
      */


      //These conditions are inflow for Shock-Vortex interaction, HiOCFD5
      //Upstream condition:
      scalar Mshock = 1.5; //mach number of the standing shock
      scalar rsh = constants::GLOBAL_GAMMA; //ratio of specific heats
      scalar Rg  = 1.0; // gas constant, demanded by workshop people
      scalar rho = 1.0;
      scalar u = Mshock * sqrt(rsh);
      scalar v = 0.0;
      scalar p = 1.0;
      

      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = rho*v; //side 1, field 2
      UF[((t*N_F+3)*2+1)*M_s+j] = p/(rsh-1.0) + 0.5*rho*(u*u + v*v); //side 1, field 3
#ifdef RADSINGLEFLUID
      //I think it's better to leave C value at interior value
      //UF[((t*N_F+4)*2+1)*M_s+j] = 0.0; //side 1, field 4
      //UF[((t*N_F+5)*2+1)*M_s+j] = 0.0; //side 1, field 5
#endif


      //printf("\t{rho, mox, moy, Et} = {%f, %f, %f, %f}\n", rho,rho*u,rho*v,Et);
     
#endif //on dimensions
      
#ifdef USE_CPU
    }
    }
#endif
}
//========================================================================================
arch_global void KJet(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement Kushner Jet inflow boundary conditions
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of Kushner Jet inflow BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED
     
      //There is no imediate 1D application for this condition
      scalar rho = 1.0;
      scalar p = 1.0;
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar sos = sqrt(gamma*p/rho);
      scalar u = 2.15*sos; //Mach 2.15 inflow
      scalar Et = p/(gamma-1.0) + 0.5*rho*u*u;
      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = Et; //side 1, field 2

#elif TWOD
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar Rgas  = constants::GLOBAL_RGAS;
      scalar RhoRatio = 0.25; //indicates density of jet compared to ambient environment
      scalar rho = RhoRatio * 1.225; //density in jet
      scalar Tfree  = 300; //ambient temperature; the jet is hotter by RhoRatio
      scalar pres = 1.225*Rgas*Tfree; //freestream pressure (jet is same pressure)
      scalar Mach = 0.25; //mach number of the jet
      scalar sos = sqrt(gamma*Rgas*Tfree);
      scalar u = 0;
      scalar v = -Mach*sos; //jet is pointed down
      scalar Et = pres/(gamma-1.0) + 0.5*rho*(u*u + v*v);
      
      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+3)*2+1)*M_s+j] = Et; //side 1, field 3
      //UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      //UF[((t*N_F+2)*2+1)*M_s+j] = rho*v; //side 1, field 2
      //New idea: try countering any momentum inside the domain by making
      //rho*v (and also rho*u, actually the AVERAGE momentum at the boundary interface
      //Didn't work: setting rho to outer state, setting Et based on contender, momentums by contender (SLAU)
      //Didn't work: setting rho and Et by outer state, momentumns by contender (SLAU)
      //Didn't work: setting rho by contender, Et based on outer state, momentums by contender (SLAU)

      //Decent configuration, but not good:
      /*
      //The rho set:
      scalar TargetAvg = rho;
      scalar Contender = UF[((t*N_F+0)*2+0)*M_s+j]; //interior value for rho
      //UF[((t*N_F+0)*2+1)*M_s+j] = 2.0 * TargetAvg - Contender; //rho from exterior side of interface

      //The rho*u set:
      TargetAvg = rho*u;
      Contender = UF[((t*N_F+1)*2+0)*M_s+j]; //interior value for momentum in x direction
      UF[((t*N_F+1)*2+1)*M_s+j] = 2.0 * TargetAvg - Contender; //rho*u from exterior side of interface

      //The rho*v set:
      TargetAvg = rho*v;
      Contender = UF[((t*N_F+2)*2+0)*M_s+j]; //interior value for momentum in y direction
      UF[((t*N_F+2)*2+1)*M_s+j] = 2.0 * TargetAvg - Contender; //rho*v from exterior side of interface

      //The Eg set:
      TargetAvg = Et;
      Contender = UF[((t*N_F+3)*2+0)*M_s+j]; //interior value for Et
      //UF[((t*N_F+3)*2+1)*M_s+j] = 2.0 * TargetAvg - Contender; //rho*v from exterior side of interface
      */
      //Another configuration: Sets off a bommb at nozzle outlet, but code remains stable, result looks cool
      /*
      scalar rhoFace = RhoRatio*1.225;
      scalar pFace = 1.225*Rgas*Tfree;
      scalar vFace = -Mach*sos;
      //New idea: Force both face states to use my prescribed BC
      UF[((t*N_F+0)*2+0)*M_s+j] = rhoFace;
      UF[((t*N_F+0)*2+1)*M_s+j] = rhoFace;
      UF[((t*N_F+1)*2+0)*M_s+j] = 0.0;
      UF[((t*N_F+1)*2+1)*M_s+j] = 0.0;
      UF[((t*N_F+2)*2+0)*M_s+j] = rhoFace*vFace;
      UF[((t*N_F+2)*2+1)*M_s+j] = rhoFace*vFace;
      UF[((t*N_F+3)*2+0)*M_s+j] = pFace/(gamma-1.0) + 0.5*rhoFace*(vFace*vFace);
      UF[((t*N_F+3)*2+1)*M_s+j] = pFace/(gamma-1.0) + 0.5*rhoFace*(vFace*vFace);
      */
      //BEST CONFIGUREATION BELOW!!!! generates shock wave at initialization, but past that everything is perfect.
      scalar rhoFace = RhoRatio*1.225;
      scalar pFace = 1.225*Rgas*Tfree;
      scalar vFace = -Mach*sos;
      //New idea: Force both face states to use my prescribed BC
      UF[((t*N_F+0)*2+0)*M_s+j] = rhoFace;
      UF[((t*N_F+0)*2+1)*M_s+j] = rhoFace;
      UF[((t*N_F+1)*2+0)*M_s+j] = 0.0;
      UF[((t*N_F+1)*2+1)*M_s+j] = 0.0;
      UF[((t*N_F+2)*2+0)*M_s+j] = rhoFace*vFace;
      UF[((t*N_F+2)*2+1)*M_s+j] = rhoFace*vFace;
      //mess with rho and veloctiy but not energy:
      //UF[((t*N_F+3)*2+0)*M_s+j] = pFace/(gamma-1.0) + 0.5*rhoFace*(vFace*vFace);
      UF[((t*N_F+3)*2+1)*M_s+j] = UF[((t*N_F+3)*2+0)*M_s+j];
      

      //UF[((t*N_F+3)*2+1)*M_s+j] = Et; //side 1, field 3

      //printf("\t{rho, mox, moy, Et} = {%f, %f, %f, %f}\n", rho,rho*u,rho*v,Et);
     
#endif //on dimensions
      
#ifdef USE_CPU
    }
    }
#endif
}

//========================================================================================
arch_global void Homo(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement Homogenous (U=2 for scalar) boundary conditions
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of Homogeneous BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    //printf("Homogeneouds BC treatment: Attacking global interface %d\n", t);
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif
#ifdef SCALARAD
#ifdef ONED
     
      //Not yet coded
#elif TWOD

      UF[((t*N_F+0)*2+1)*M_s+j] = 2.0; //side 1, field 0

      //printf("\t{rho, mox, moy, Et} = {%f, %f, %f, %f}\n", rho,rho*u,rho*v,Et);
     
#endif //on dimensions
#endif //on SCALARAD
      //singlefluid and radsinglefluid: outflow condition for the shock-vortex interaction
#ifdef SINGLEFLUID
      scalar rsh = constants::GLOBAL_GAMMA; //ratio of specific heats
      scalar Mshock = 1.5; //mach number of the standing shock
      //Upstream condition:
      scalar rho0 = 1.0;
      scalar u0 = Mshock * sqrt(rsh);
      scalar v0 = 0.0;
      scalar p0 = 1.0;
      
      //Upstream vs downstream ratios:
      scalar ruby = rsh - 1.0;
      scalar chip = rsh + 1.0;
      scalar RhoRatio = (2.0 + ruby*Mshock*Mshock) / (chip*Mshock*Mshock);
      scalar pRatio = 1.0 + 2*rsh/chip*(Mshock*Mshock-1.0);
      scalar uRatio = 1.0 / RhoRatio;
      scalar vRatio = 1.0;
      
      //Downstream condition:
      scalar rho = rho0 / RhoRatio;
      scalar u   = u0   / uRatio;
      scalar v   = v0   / vRatio;
      scalar p   = p0   * pRatio;

      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = rho*v; //side 1, field 2
      UF[((t*N_F+3)*2+1)*M_s+j] = p/(rsh-1.0) + 0.5*rho*(u*u + v*v);; //side 1, field 3
#endif //for singlefluid
#ifdef RADSINGLEFLUID
      scalar rsh = constants::GLOBAL_GAMMA; //ratio of specific heats
      scalar Mshock = 1.5; //mach number of the standing shock
      //Upstream condition:
      scalar rho0 = 1.0;
      scalar u0 = Mshock * sqrt(rsh);
      scalar v0 = 0.0;
      scalar p0 = 1.0;
      
      //Upstream vs downstream ratios:
      scalar ruby = rsh - 1.0;
      scalar chip = rsh + 1.0;
      scalar RhoRatio = (2.0 + ruby*Mshock*Mshock) / (chip*Mshock*Mshock);
      scalar pRatio = 1.0 + 2*rsh/chip*(Mshock*Mshock-1.0);
      scalar uRatio = 1.0 / RhoRatio;
      scalar vRatio = 1.0;
      
      //Downstream condition:
      scalar rho = rho0 / RhoRatio;
      scalar u   = u0   / uRatio;
      scalar v   = v0   / vRatio;
      scalar p   = p0   * pRatio;

      UF[((t*N_F+0)*2+1)*M_s+j] = rho; //side 1, field 0
      UF[((t*N_F+1)*2+1)*M_s+j] = rho*u; //side 1, field 1
      UF[((t*N_F+2)*2+1)*M_s+j] = rho*v; //side 1, field 2
      UF[((t*N_F+3)*2+1)*M_s+j] = p/(rsh-1.0) + 0.5*rho*(u*u + v*v);; //side 1, field 3
      //Leave C at interior values
      //UF[((t*N_F+4)*2+1)*M_s+j] = 0.0; //side 1, field 4
      //UF[((t*N_F+5)*2+1)*M_s+j] = 0.0; //side 1, field 5
#endif //for radsinglefluid
#ifdef USE_CPU
    }
    }
#endif
}

arch_global void SubOut(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Implement subsonic outflow boundary conditions. MUST SET AMBIENT CONDITION HERE
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of subsonic outflow BC in boundaryMap
    \param[out] UF suitably modified for BCs
  */
#ifdef USE_CPU
  for(int k = 0; k < M_B; k++){
    int t = boundaryMap[start+k];
    for(int j = 0; j < M_s; j++){
#elif USE_GPU
      int t = boundaryMap[start+blockIdx.x];
      int j = threadIdx.x;
#endif

#ifdef ONED

#elif TWOD
      scalar gamma = constants::GLOBAL_GAMMA;
      scalar Rgas  = constants::GLOBAL_RGAS;

      //Coming from inside the domain:
      scalar rho_IN = UF[((t*N_F+0)*2+0)*M_s+j];
      scalar mox_IN = UF[((t*N_F+1)*2+0)*M_s+j];
      scalar u_IN = mox_IN / rho_IN;
      scalar moy_IN = UF[((t*N_F+2)*2+0)*M_s+j];
      scalar v_IN = moy_IN / rho_IN;
      scalar Eg_IN =  UF[((t*N_F+3)*2+0)*M_s+j];
      scalar p_IN = (gamma-1.0) * (Eg_IN - 0.5*rho_IN*(u_IN*u_IN + v_IN*v_IN));

      //ambient state:
      scalar rho_OUT = 1.225;
      scalar u_OUT = 0.0;
      scalar v_OUT = 0.0;
      scalar T_OUT = 300.0;
      scalar p_OUT = 1.225 * Rgas * T_OUT;

      //Configuration 1:
      //interface pressure = ambient
      //density and velocity come from inside domain
      scalar rho_FACE = rho_IN;
      scalar mox_FACE = mox_IN;
      scalar moy_FACE = moy_IN;
      scalar Eg_FACE = p_OUT/(gamma-1.0) + 0.5*rho_IN*(u_IN*u_IN + v_IN*v_IN);

      UF[((t*N_F+0)*2+0)*M_s+j] = rho_FACE;
      UF[((t*N_F+0)*2+1)*M_s+j] = rho_FACE;
      UF[((t*N_F+1)*2+0)*M_s+j] = mox_FACE;
      UF[((t*N_F+1)*2+1)*M_s+j] = mox_FACE;
      UF[((t*N_F+2)*2+0)*M_s+j] = moy_FACE;
      UF[((t*N_F+2)*2+1)*M_s+j] = moy_FACE;
      UF[((t*N_F+3)*2+0)*M_s+j] = Eg_FACE;
      UF[((t*N_F+3)*2+1)*M_s+j] = Eg_FACE;
      

      //UF[((t*N_F+3)*2+1)*M_s+j] = Et; //side 1, field 3

      //printf("\t{rho, mox, moy, Et} = {%f, %f, %f, %f}\n", rho,rho*u,rho*v,Et);
     
#endif //on dimensions
      
#ifdef USE_CPU
    }
    }
#endif
}

//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C"
void LrflctiveBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Host C function to launch rflctive kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of reflective BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

  rflctive arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}

extern "C"
void LnoslipBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Host C function to launch noslip kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of no-slip boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of noslip BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

 noslip arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}

extern "C"
void LnogradBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* dUgF){
  /*!
    \brief Host C function to launch nograd kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of zero-gradient boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of zero-gradient BC in boundaryMap
    \param[out] dUgF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

  nograd arch_args (M_s, M_B, boundaryMap, normals, start, dUgF);
}

extern "C"
void LAnflwBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Host C function to launch Anflw kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of A inflow boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of A inflow BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

 Anflw arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}

extern "C"
void LKJetBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Host C function to launch KJet kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of KJet inflow boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of KJet inflow BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

 KJet arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}

extern "C"
void LSubOutBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF){
  /*!
    \brief Host C function to launch SubOut kernel.
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of SubOut boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of SubOut BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

 SubOut arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}

extern "C"
void LHomoBoundary(int M_s, int M_B, int* boundaryMap, scalar* normals, int start, scalar* UF)
{
    /*!
    \brief Host C function to launch  Homo kernel (homogeneous boundary condition).
    \param[in] M_s number of nodes on an interface
    \param[in] M_B number of boundaries to treat; specifically, the number of Homogeneous dirichlet boundaries
    \param[in] boundaryMap index of interfaces to treat
    \param[in] normals normals to the interfaces
    \param[in] start start index of Homogeneous BC in boundaryMap
    \param[out] UF suitably modified for BCs
    \section Description
    In GPU mode, launches M_B blocks of M_s threads.
  */
#ifdef USE_GPU
  dim3 dimBlock(M_s,1,1);
  dim3 dimGrid(M_B,1);
#endif

 Homo arch_args (M_s, M_B, boundaryMap, normals, start, UF);
}
