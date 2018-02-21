/*!
  \file error_quant.cc
  \brief Initial condition function definitions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "error_quant.h"
#include <stdlib.h>

void Exact_sinphil(scalar* XYZ, scalar* sol, scalar Ldomain)
{
  /*!
    brief: return the sinphil solution at specified location
    \param[in] N_F number of fields
    \param[in] XYZ the location
    \param[out] sol the solution, one entry per field variable
  */
  
  for (int fc = 0; fc < N_F; fc++)
    {
#ifdef ONED
      sol[fc] = 2.0 + sin(2.0*M_PI/Ldomain * XYZ[0]);
#endif
#ifdef TWOD
      sol[fc] = 2.0 + sin(2.0*M_PI/Ldomain * XYZ[0]) * sin(2.0*M_PI/Ldomain * XYZ[1]);
#endif
#ifdef THREED
      sol[fc] = 2.0 + sin(2.0*M_PI/Ldomain * XYZ[0]) * sin(2.0*M_PI/Ldomain * XYZ[1]) * sin(2.0*M_PI/Ldomain * XYZ[2]);
#endif
    }
  
}

void Exact_zero(scalar* XYZ, scalar* sol)
{
  /*!
    brief: return the zerosolution at specified location
    \param[in] N_F number of fields
    \param[in] XYZ the location
    \param[out] sol the solution, one entry per field variable
  */
  
  for (int fc = 0; fc < N_F; fc++)
    {
      sol[fc] = 0.0;

    }
}

void Exact_normvtx(scalar* XYZ, scalar* sol)
{
  /*!
    brief 2D isentropic vortex
   */
  //Doreen vortex parameters
  scalar rsh = 1.4; //ratio of specific heats, needs to match with CONSTANTS:gamma
  scalar K_rsh = 1.0 / rsh; //isentropic constant
  scalar k_vort = 2.0;
  scalar u_Doreen = 1.0;
  scalar v_Doreen = 0.0;
  scalar rho_Doreen = 1.0;
  scalar exp_const = 2.71828182845904523536;
  scalar x_loc = XYZ[0];
  scalar y_loc = XYZ[1];
  scalar x_sq = x_loc*x_loc;
  scalar y_sq = y_loc*y_loc;
  scalar kxy = k_vort * (x_sq + y_sq);
  scalar prod1 = (rsh-1.0)*pow(exp_const, -2.0*kxy) / (4.0*rsh*k_vort*K_rsh);
  scalar diff1 = pow(rho_Doreen, rsh-1.0) - prod1;
  scalar rho = pow(diff1, 1.0/(rsh-1.0));
  scalar u = u_Doreen - y_loc * pow(exp_const, -kxy);
  scalar v = x_loc * pow(exp_const, -kxy);
  scalar p = K_rsh * pow(rho, rsh);
  scalar Et = p/(rsh-1.0) + 0.5*rho*(u*u + v*v);
  //printf("e=%d, i=%d, (x,y) = (%f, %f); (rho,u,v,p) = ( %f, %f, %f, %f)\n",e,i,x,y,rho,u,v,p);
  // set the initial fields

  //Again, this setup is built to handle 2D and 3D cases
  sol[0] = rho;
  sol[1] = rho*u;
  sol[2] = rho*v;
#ifdef THREED
  sol[3] = 0.0;
#endif
  sol[N_F-1] = Et;
}

void Exact_HiOWvtx(scalar* XYZ, scalar* sol)
{
  // For the High-Order Workshop problem VI1 in 2017
  // Vortex transport by uniform flow
  scalar exp_const = 2.71828182845904523536;
  //vortex parameters
  scalar Mach_free = 0.5;
  scalar p_free = 100000;
  scalar T_free = 300.0;
  scalar Rg = 287.15;
  scalar rsh = 1.4;
  scalar cp = Rg*(rsh)/(rsh-1.0);
  scalar rho_free = p_free/(Rg*T_free);
  scalar a_free = sqrt(rsh*p_free / rho_free);
  scalar u_free = Mach_free*a_free;
  scalar v_free = 0.0;
  scalar X_epi = 0.05;
  scalar Y_epi = 0.05;
  scalar R_vortex = 0.005;
  scalar Beta_vortex = 0.2;


  scalar x_loc = XYZ[0];
  scalar y_loc = XYZ[1];
  scalar x_sq = (x_loc-X_epi) * (x_loc-X_epi);
  scalar y_sq = (y_loc-Y_epi) * (y_loc-Y_epi);
  scalar r_epi = sqrt( x_sq + y_sq ) / R_vortex;
  scalar r_sq = r_epi*r_epi;
	      
  scalar u_pert = (-u_free*Beta_vortex) * (y_loc-Y_epi) / R_vortex * pow(exp_const,-0.5*r_sq);
  scalar v_pert = (u_free*Beta_vortex) * (x_loc-X_epi) / R_vortex * pow(exp_const,-0.5*r_sq);
  scalar T_pert =  0.5*pow(u_free*Beta_vortex , 2.0) / cp * pow(exp_const,-r_sq); 

  scalar T_out = T_free - T_pert;
  scalar u_out = u_free + u_pert;
  scalar v_out = 0.0 + v_pert;
  scalar rho_out = rho_free * pow(T_out/T_free , 1.0/(rsh-1.0));
  scalar p_out = rho_out * Rg * T_out;


  //conserved variables
  sol[0] = rho_out;
  sol[1] = rho_out*u_out;
  sol[2] = rho_out*v_out;
  sol[N_F-1] = p_out/(rsh-1.0) + 0.5*rho_out*(u_out*u_out + v_out*v_out);
#ifdef THREED
  sol[3] = 0.0;
#endif
}
