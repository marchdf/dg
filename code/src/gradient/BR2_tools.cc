/*!
  \file BR2_tools.cc
  \BR2 operations, including local lift functions and gradient calculation
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/

#include "BR2_tools.h"

void Jump_Residual(int M_T, int N_s, int M_G, int N_N, int* BR2_Map, scalar* UgF, scalar* phigF, scalar* normals, scalar* Jac_trace, scalar* JF, scalar* weight_trace, scalar* BR2_resi, int ruk)
{
  int print_jump = 0;
  /*!
    This function is interface-oriented, but relays cell-specific information
    UgF is solution along both sides of each interface, interface oriented
    phigF is basis trace from elements in each side of interface, interface oriented
    normals is outward normal from cell A to cell B, I hope.
    Lift_vec is the residual of the system to solve for the local lift on every face of every cell, cell oriented 
  */
  //printf("Entered Jump_Residual Routine\n");

#ifdef ONED
  {
    //int D = 1;
  
    for (int t = 0; t < M_T; t++) //interface
      {
	//The interface needs to know the bordering cells, to grab lift data
	int om_A, sig_A, om_B, sig_B;
	om_A = BR2_Map[t*4 + 0]; //cell on one side of interface
	sig_A = BR2_Map[t*4 + 1]; //reference element face (0<sig<2*N_N) linking om_A to the interface
	om_B = BR2_Map[t*4 + 2]; //cell on the other side of interface
	sig_B = BR2_Map[t*4 + 3]; //reference element face (0<sig<2*N_N) linking om_B to the interface;
		     
	scalar ndir = normals[t*D+0];    //the normal needs to run from cell A to cell B
	  
	scalar integrand_A;
	scalar integrand_B;
	scalar Lift_A;
	scalar Lift_B;
	scalar diff;
	for (int alpha = 0; alpha < D; alpha = alpha + 1) //Direction
	  {
	    for (int k = 0; k < N_s; k = k + 1) //test basis index
	      {
		for (int field = 0; field < N_F; field = field + 1) //field (for example, rho or rho*u)
		  {
		    Lift_A = 0.0;
		    Lift_B = 0.0;
		    //quadrature node = g = 0
		    //Use the UA-UB difference and normal from A to B
		    //0 is for cell A, 1 is for cell B
		    //Originally had negative sign here, might need a positive
		    diff = -0.5*(UgF[((t*N_F+field)*2+0)*M_G+0] - UgF[((t*N_F+field)*2+1)*M_G+0]) * ndir;
		  
		    integrand_A = phigF[k*D*N_N*M_G + sig_A*M_G + 0] * diff;
		    integrand_B = phigF[k*D*N_N*M_G + sig_B*M_G + 0] * diff;
		    Lift_A = integrand_A;
		    Lift_B = integrand_B;
		    
		    //Populate the output:
		    BR2_resi[om_A*N_N*N_F*D*N_s + (sig_A%N_N)*N_F*D*N_s + field*D*N_s + alpha*N_s + k] = Lift_A;
		    BR2_resi[om_B*N_N*N_F*D*N_s + (sig_B%N_N)*N_F*D*N_s + field*D*N_s + alpha*N_s + k] = Lift_B;
		  }
	      }
	  }
	//k = test basis number
	//KDG = dof per cell
	//D = geometric dimension
	//sig_ = index of the interface from cell's perspective
	//N_N = number of faces per element
      } 
  }
#elif TWOD
  //int D = 2;
  for (int t = 0; t < M_T; t++) //interface
    {
      //The interface needs to know the bordering cells, to grab lift data
      int om_A, sig_A, om_B, sig_B;
      om_A = BR2_Map[t*4 + 0]; //cell on one side of interface
      sig_A = BR2_Map[t*4 + 1]; //reference element face linking om_A to the interface
      om_B = BR2_Map[t*4 + 2]; //cell on the other side of interface
      sig_B = BR2_Map[t*4 + 3]; //reference element face linking om_B to the interface;
      scalar ndir_phys[2];
      ndir_phys[0] = normals[t*D+0]; //the normal needs to run from cell A to cell B
      ndir_phys[1] = normals[t*D+1];
      //In an attempt to get the 2D version to work properly, adding in 
      //some scaling factors for the face normal
      scalar ndir_br2[2];
      //plugging in a zero for g argument because Jacobian is uniform along a face, I hope
      scalar dx_dxi =  Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + 0*D*D + 0*D + 0];//0.05/2.0;
      scalar dy_deta = Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + 0*D*D + 1*D + 1];//0.33/2.0;
      scalar dx_deta = Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + 0*D*D + 1*D + 0];//0.0;
      scalar dy_dxi =  Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + 0*D*D + 0*D + 1];//0.0;
      //the 2 is because the reference element has length 2
      ndir_br2[0] = ndir_phys[0]*2.0*dx_dxi  + ndir_phys[1]*2.0*dy_dxi;
      ndir_br2[1] = ndir_phys[0]*2.0*dx_deta + ndir_phys[1]*2.0*dy_deta;
      scalar integrand_A;
      scalar integrand_B;
      scalar Lift_A;
      scalar Lift_B;
      scalar diff;
      //For 2D case, need to use the Jacobian of the face to proportion quadrature
      //Mark uses face Jacobian in the redirstribute_q function
      scalar Jac_edge_A = JF[t*2 + 0];
      scalar Jac_edge_B = JF[t*2 + 1]; //ask Marc abourt using these values. Indexing likely wrong
      if (ruk == 1 && print_jump == 1)
	{
	  printf("Jump_Residaul Attacking Interface %d: Jac_edgeA = %f, Jac_EdgeB = %f\n",t, Jac_edge_A, Jac_edge_B);
	  printf("More info: dx_dxi = %f, dx_deta = %f, dy_dxi = %f, dy_deta = %f\n",dx_dxi,dx_deta,dy_dxi,dy_deta);
	}
      for (int alpha = 0; alpha < D; alpha = alpha + 1) //Direction
	{
	  if (ruk == 1 && print_jump == 1)
	    {
	      printf("\talpha = %d\n",alpha);
	    }
	  //printf("Jump_Residaul Attacking Interface %d, alpha = %d\n",t,alpha);
	  for (int k = 0; k < N_s; k = k + 1) //test basis index
	    {
	      /*
	      if (ruk == 1 && print_jump == 1)
		{
		  printf("\t\tk = %d\n",k);
		}
	      */
	      //printf("Jump_Residaul Attacking Interface %d, alpha = %d, DOF = %d\n",t,alpha,k);
	      for (int field = 0; field < N_F; field = field + 1) //field (for example, rho or rho*u)
		{
		  //printf("Jump_Residaul Attacking Interface %d, alpha = %d, DOF = %d, Field = %d\n",t,alpha,k,field);
		  Lift_A = 0.0;
		  Lift_B = 0.0;
		  for (int g = 0; g < M_G; g = g + 1) //quadrature node
		    {
		      //printf("Jump_Residaul Attacking Interface %d, alpha = %d, DOF = %d, Field = %d, quad_node = %d\n",t,alpha,k,field,g);
		      //Use the UA-UB difference and normal from A to B
		      //0 is for cell A, 1 is for cell B
		      diff = -0.5*(UgF[((t*N_F+field)*2+0)*M_G+g] - UgF[((t*N_F+field)*2+1)*M_G+g]) * ndir_phys[alpha];
		      /*
		      if (ruk == 0 && fabs(diff) > pow(10,-10) && k == 0 && field == 0)
			{
			  printf("Substantial diff discovered: t=%d,om_A=%d,om_B=%d,field=%d,g=%d\n",t,om_A,om_B,field,g); 
			}
		      */
		      /*
		      if (ruk == 0 && ndir[1] == 1.0 && diff != 0)
			{
			  printf("Detected nonzero y jump, omA=%d,om_B=%d, field=%d, gnode=%d, Jump=%15.14f\n",om_A,om_B,field,g,diff);
			}
		      */
		      if (t == 2)
			{
			  //printf("Accessed UgF Succussfully, diff = %f\n", diff);
			}
		      integrand_A = phigF[k*D*N_N*M_G + sig_A*M_G + g] * diff;
		       if (t==2)
			{
			  //printf("Accessed phigF successfully\n");
			}
		      integrand_B = phigF[k*D*N_N*M_G + sig_B*M_G + g] * diff;
		     
		      //Lift_A += integrand_A * GQ_w[g] * Jac_edge_A;
		      //Lift_B += integrand_B * GQ_w[g] * Jac_edge_B;
		      Lift_A += integrand_A * weight_trace[sig_A*M_G + g] * Jac_edge_A;
		      Lift_B += integrand_B * weight_trace[sig_B*M_G + g] * Jac_edge_B;

		      
		      /*
		      if (ruk == 1 && print_jump == 1 && field == 3)
			{
			  printf("\t\t\t{g, diff, int_A, int_B} = {%d, %f, %f, %f}\n",g,diff,integrand_A, integrand_B);
			}
		      */
		    }
		  /*
		  //Setting some jumps to zero:
		  //if (fabs(ndir[0]) > 0.5) //sets x jumps to zero
		  if (fabs(ndir[1]) > 0.5) //sets y jumps to zero
		    {
		      Lift_A *= 0.0;
		      Lift_B *= 0.0;
		    }
		  */

		  BR2_resi[om_A*N_N*N_F*D*N_s + (sig_A%N_N)*N_F*D*N_s + field*D*N_s + alpha*N_s + k] = Lift_A;
		  BR2_resi[om_B*N_N*N_F*D*N_s + (sig_B%N_N)*N_F*D*N_s + field*D*N_s + alpha*N_s + k] = Lift_B;
		}
	    }
	}
      //k = test basis number
      //N_s = dof per cell
      //D = geometric dimension
      //sig_ = index of the interface from cell's perspective
      //N_N = number of faces per element
    }

#endif
  //printf("Finished Jump_Residual Routine\n");
}
//============================================================
//Lift_Solve uses the BR2_resi results to get the lift weights
void Lift_Solve(int N_E, int N_N, int N_s, scalar* h_Minv, scalar* BR2_resi, scalar* BR2_poly, int ruk)
{
  /*!
  N_E = element count
  N_N = faces per element
  N_s = approximation order
  h_Minv = the inverse mass matrix in single-column form
  BR2_resi: the residual in the br2 gradient correction system
  This function is cell-oriented
  BR2_poly is the column representation of lift weights, in physical xy coordinates
  the coordinate system matters greatly because br2_poly represents a gradient

  NOTE: If SSAM is the unit mass matrix, BR2_poly will be wrt reference coordinates
  if SSAM has already accounted for geometry, then BR2_poly is wrt physical coordinates
  01/04/2016: Learned today that the SSAM is wrt physical coordinates
  */

  //Wherever you see if ruk==, it's because I only want to see what's going on at one RK substep
  /*
  if (ruk == 1)
    {
      printf("Entered Lift_Solve\n");
    }
  */
  int print_lift = 0;
  for (int e = 0; e < N_E; e = e + 1)
    {
      
      if (ruk == 1 && print_lift == 1)
	{
	  printf("--Element %d:---\n",e);
	  /*
	  for (int k = 0; k < N_s; k++)
	    {
	      for (int j = 0; j < N_s; j++)
		{
		  printf("SSAM[%d][%d] = %f\n", k,j,h_Minv[e*N_s*N_s + j*N_s + k]);
		}
	    }
	  */
	  
	  for (int sig = 0; sig < N_N; sig++)
	    {
	      for (int alpha = 0; alpha < D; alpha++)
		{
		  for (int field = 0; field < N_F; field++)
		    {
		      for (int j = 0; j < N_s; j++)
			{
			  
			  if (field == 3 && print_lift == 1)// && BR2_resi[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + j] != 0)
			    {
			      printf("\t\tBR2_resi[e=%d][sig=%d][alpha=%d][field=%d][j=%d] = %f\n", e,sig,alpha,field,j , BR2_resi[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + j]);
			    }
			  
			}
		    }
		}
	    }
	}
      
      for (int sig = 0; sig < N_N; sig = sig + 1)
	{
	  for (int alpha = 0; alpha < D; alpha = alpha + 1)
	    {
	      for (int field = 0; field < N_F; field = field + 1)
		{
		  //k,j indexing is tricky here because of matrix multiplication
		  //In abstract terms, j is for the column of SSAM, row of Resi, while k is for the row of SSAM.
		  for (int k = 0; k < N_s; k = k + 1)
		    {
		      scalar sol = 0.0;
		      for (int j = 0; j < N_s; j = j + 1)
			{
			  //sol += h_Minv[e*N_s*N_s + j*N_s + k] * BR2_resi[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + j];
			  sol += h_Minv[e*N_s*N_s + j*N_s + k] * BR2_resi[e*N_N*N_F*D*N_s + sig*N_F*D*N_s + field*D*N_s + alpha*N_s + j];
			}
		      BR2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + k] = sol;
		       
		      if (field == 3 && ruk == 1 && print_lift == 1)// && BR2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + k] != 0)
			{
			  printf("\t\t\t\tBR2_poly[e=%d][sig=%d][alpha=%d][f=%d][k=%d] = %f\n",e,sig,alpha,field,k,BR2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + alpha*N_F*N_s + field*N_s + k]);
			}
		      
		    }
		}
	    }
	}
    }
  /*
  if (ruk == 1)
    {
      printf("Left Lift_Solve\n");
    }
  */
}

void populate_CellGrad(int N_E, int N_G, int N_s, int N_N, scalar* U, scalar* invJac,  scalar* dphig, scalar* phig, scalar* br2_poly, scalar* dUg, int ruk)
{
  /*!
  This function is cell-oriented
  U is cell solution weights
  N_s is degrees of freedom per element
  dphig is  gradient of shape function traces wrt xi, eta
  phig is solution basis at quadrature nodes
  br2_poly is the weights for the gradient correction, wrt coordinates determined by SSAM
  dUg is the output, corrected gradient at interior quadrature points wrt xi,eta
  */
  //printf("Entered Populate_CellGrad\n");
  //12/11/2015: I'm still thinking that the output from this function needs to be a gradient wrt physical coordinates
  //01/01/2016: Definitely thinking that the output from here should be wrt physical coordinates. ATTEND TO PHYSICS.CU
  //01/04/2016: Attended to Physics.cu.
#ifdef ONED
  {
    //int D = 1;
    scalar corr_mult = 1.0;
    scalar corr_dUdx;
    scalar corr_dUdxi;
    scalar naive_dUdxi;
    scalar dUdxi;
    scalar naive_dUdx;
    scalar dUdx;
    scalar dxi_dx;
    
    for (int e = 0; e < N_E; e++)
      {
	for (int field = 0; field < N_F; field++)
	  {
	    /*
	      if (ruk == 1 && field == 2)
	      {
	      printf("-----For element %d, field = %d:-----\n", e, field);
	      }
	    */
	    for (int g = 0; g < N_G; g++)
	      {
		corr_dUdx = 0;
		corr_dUdxi = 0;
		naive_dUdxi = 0;
		//dx_dxi = Jac(e, g); //I'm grabbing in actual matrix form, because the column-major form is not populated in main
		//dxi_dx = invJac(e*D+0, g*D+0);
		dxi_dx = invJac[e*N_G*D*D + 0 + 0 + 0];
		//Populate naive cell gradient
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    naive_dUdxi += U[e*N_F*N_s + field*N_s + k] * dphig[k*N_G*D + g*D + 0];
		  }
		//Populate Correction
		for (int sig = 0; sig < N_N; sig = sig + 1)
		  {
		    for (int k = 0; k < N_s; k = k + 1)
		      {
			//The zero in br2_poly index is for dUdx, as opposed to dU in other directions
			corr_dUdx += br2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phig[k*N_G + g];
			//corr_dUdx += br2_poly[e*N_N*N_F*D*N_s + (sig%N_N)*(N_F*D*N_s) + field*(D*N_s) + 0*N_s + k] * phig[k*N_G + g];
		      }
		  }
		corr_dUdx = corr_dUdx * corr_mult;
		/*
		//Move correction to reference coordinates
		corr_dUdxi =  corr_dUdx * dx_dxi;
		naive_dUdx = naive_dUdxi * dxi_dx;
		
		//Now, calculate the corrected gradient
		dUdxi =  naive_dUdxi  + corr_dUdxi;
		*/
		//MOve naive gradient tp physical coordinates, then make correction
		naive_dUdx = naive_dUdxi * dxi_dx;
		dUdx = naive_dUdx + corr_dUdx;
		/*  
		    if (field == 2 && ruk == 1)
		    {
		    printf("node[%d]: naive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",g, naive_dUdx, corr_dUdx, dUdx);
		    }
		*/
		//dUdx = naive_dUdx + corr_dUdx*dxi_dx*0.01;
		dUg[e*N_F*N_G*D + field*N_G*D + g*D + 0] = dUdx;
	      //dUg[e*(N_F*N_G*D) + field*(N_G*D) + g*D + 0] = dUdxi;
	      }
	  }
      }
  }
#elif TWOD
  {
    //int D = 2;
    int print_cell = 0;
    scalar corr_mult = 1.0;
    scalar corr_dUdx;
    scalar corr_dUdy;
    scalar corr_dUdxi;
    scalar corr_dUdeta;
    scalar naive_dUdxi;
    scalar naive_dUdeta;
    scalar naive_dUdx;
    scalar naive_dUdy;
    scalar dUdxi;
    scalar dUdeta;
    /*
      scalar dx_dxi;
      scalar dy_dxi;
      scalar dx_deta;
      scalar dy_deta;
    */
    for (int e = 0; e < N_E; e++)
      {
	
	if (ruk == 3 && print_cell == 1)
	  { 
	    printf("---Element %d Interior:---\n", e);
	  }
	
	for (int field = 0; field < N_F; field++)
	  {
	    /*
	    if (ruk == 3 && print_cell == 1)
	      {
		printf("\tField = %d:\n", field);
	      }
	    */
	    for (int g = 0; g < N_G; g++)
	      {
		corr_dUdx = 0;
		corr_dUdy = 0;
		corr_dUdxi = 0;
		corr_dUdeta = 0;
		naive_dUdxi = 0;
		naive_dUdeta = 0;
	      
		scalar dxi_dx  = invJac[e*N_G*D*D + g*D*D + 0 + 0];
		scalar deta_dx = invJac[e*N_G*D*D + g*D*D + 0 + 1];
		scalar dxi_dy  = invJac[e*N_G*D*D + g*D*D + D + 0];
		scalar deta_dy = invJac[e*N_G*D*D + g*D*D + D + 1];
		
		if (ruk == 3 && field == 0 && print_cell == 1)
		  {
		    printf("\tnode[%d]: dxi_dx = %f, deta_dx = %f, dxi_dy = %f, deta_dy = %f\n",g,dxi_dx,deta_dx,dxi_dy,deta_dy);
		  }
		
		//Populate naive cell gradient
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    naive_dUdxi  += U[e*N_F*N_s + field*N_s + k] * dphig[k*N_G*D + g*D + 0];
		    naive_dUdeta += U[e*N_F*N_s + field*N_s + k] * dphig[k*N_G*D + g*D + 1];
		  }
		//Populate Correction
		for (int sig = 0; sig < N_N; sig = sig + 1)
		  {
		    for (int k = 0; k < N_s; k = k + 1)
		      {
			corr_dUdx += br2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phig[k*N_G + g];// / dxi_dx;
			corr_dUdy += br2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + 1*N_F*N_s + field*N_s + k] * phig[k*N_G + g];// / deta_dy;
			//corr_dUdx += br2_poly[e*(N_N*N_F*D*N_s) + sig*(N_F*D*N_s) + field*(D*N_s) + 0*N_s + k] * phig[k*N_G + g];
			//corr_dUdy += br2_poly[e*(N_N*N_F*D*N_s) + sig*(N_F*D*N_s) + field*(D*N_s) + 1*N_s + k] * phig[k*N_G + g];
			//sum_dU += br2_poly[e*(Chi*N_F*D*KDG) + sig*(N_F*D*KDG) + field*(D*KDG) + alpha*KDG + k] * phig[g,k];
		      }
		  }
		/*
		//Move correction to reference coordinates
		corr_dUdxi =  corr_dUdx * dx_dxi  + corr_dUdy * dy_dxi;
		corr_dUdeta = corr_dUdx * dx_deta + corr_dUdy * dy_deta;

		//Now, calculate the corrected gradient
	      
		dUdxi =  naive_dUdxi  + corr_dUdxi;
		dUdeta = naive_dUdeta + corr_dUdeta;
		dUg[e*(N_F*N_G*D) + field*(N_G*D) + g*D + 0] = dUdxi;
		dUg[e*(N_F*N_G*D) + field*(N_G*D) + g*D + 1] = dUdeta;
		*/

		//MOve naive gradient to physical coordinates, then apply correction
		corr_dUdx = corr_dUdx * corr_mult;
		corr_dUdy = corr_dUdy * corr_mult;
		/*
		if (ruk == 1 && corr_dUdy != 0.0)
		  {
		    printf("Nonzero dUdy Cell correction encountered: e=%d,f=%d,g=%d\n", e,field, g);
		  }
		*/
		naive_dUdx = naive_dUdxi*dxi_dx + naive_dUdeta*deta_dx;
		naive_dUdy = naive_dUdxi*dxi_dy + naive_dUdeta*deta_dy;
		scalar dUdx = naive_dUdx + corr_dUdx;
		scalar dUdy = naive_dUdy + corr_dUdy;
		
		  if (field == 0 && ruk == 3 && print_cell == 1)
		  {
		  printf("\t\tnaive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",naive_dUdx, corr_dUdx, dUdx);
		  printf("\t\tnaive_dUdy = %f, corr_dUdy = %f, dUdy = %f\n",naive_dUdy, corr_dUdy, dUdy);
		  }
		
		dUg[e*N_F*N_G*D + field*N_G*D + g*D + 0] = dUdx;
		dUg[e*N_F*N_G*D + field*N_G*D + g*D + 1] = dUdy;
	      }
	  }
      }
  }
#endif
  //printf("Left Populate_CellGrad\n");
}


//======================================================== 
void populate_FaceGrad(int M_T, int M_G, int N_s, int N_N, int* BR2_Map, scalar* U,  scalar* inv_Jac_trace, scalar* phigF, scalar* dphigF, scalar* br2_poly, scalar* dUgF, int ruk)
{
  int print_face = 0;
  /*!
  Interface-Oriented function to populate the corrected gradient along cell interfaces.
  Because of two cells, there is a big jump across the array when switching from A
  gradient population to B gradient population
  N_s is degrees of freedom per element
  U is solution weights, I hope
  phigF is basis traces from the cells at the interface quadrature points
  dphigF is the derivatives of phigF wrt xi,eta
  br2_poly is the column vector of the weights of the local lift functions, wrt x,y
  dUgF is the solution gradient on each side of the interface wrt xi,eta
  */
  //printf("Entered populate_FaceGrad\n");
#ifdef ONED
  {
    //int D = 1;
    scalar corr_mult = 1.0;
    scalar Chi_Mult = N_N*1.0;
    /*
      if (ruk == 1)
      {
      printf("Entered populate_FaceGrad, Chi_Mult = %f\n",Chi_Mult);
      }
    */
    int om_A;
    int sig_A;
    int om_B;
    int sig_B;

    scalar corr_dUdx;
    scalar corr_dUdxi;
    scalar naive_dUdxi;
    scalar naive_dUdx;
    scalar dUdxi;					
    scalar dx_dxi;
    scalar dxi_dx;
    scalar dUdx;
 
   
    for (int t = 0; t < M_T; t = t + 1) //interface
      {
	//printf("---For interface %d:---\n",t);
	//The interface needs to know the bordering cells, to grab lift data
	/*
	//See dg/dg_functions.cc , dg_mappings for how to get om_A, om_B
	const simpleInterface &face = interfaces[t];   // get the interface
	const simpleElement *el0 = face.getElement(0); // get the element to the left
	const simpleElement *el1 = face.getElement(1); // get the element to the right
	om_A = ElementMap.at(el0->getId()); //cell 0 
	om_B = ElementMat.at(el1->getId()); //cell 1
	//The ClosureID call may return any of 2*N_N numbers based on interface direction; must mod down to appropriate face designation
	sig_A = face.getClosureId(0) % N_N; //interface index from cell A perspective
	sig_B = face.getClosureId(1) % N_N; //interface index from cell B perspective
	*/
	om_A = BR2_Map[t*4 + 0]; //cell on one side of interface
	sig_A = BR2_Map[t*4 + 1]; //reference element face linking om_A to the interface
	om_B = BR2_Map[t*4 + 2]; //cell on the other side of interface
	sig_B = BR2_Map[t*4 + 3]; //reference element face linking om_B to the interface;
	for (int field = 0; field < N_F; field = field + 1)
	  {
	    /*
	      if (ruk == 1 && field == 2)
	      {
	      printf("-----For interface %d, field = %d:-----\n",t,field);
	      }
	    */
	    for (int g = 0; g < M_G; g = g + 1) //quadrature node on interface
	      {
	      
		//___________________________________
		//A Cell: //get the dx_dxi transformation
		//dx_dxi = Jac_trace(om_A*D + 0, sig_A*M_G*D + g*D + 0); //Must create this array
		dxi_dx = inv_Jac_trace[om_A*N_N*M_G*D*D + sig_A*M_G*D*D + g*D*D + 0*D + 0];
	      
		naive_dUdxi = 0.0;
		corr_dUdx = 0.0;
		//Populate naive cell gradient
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    naive_dUdxi += U[om_A*N_F*N_s + field*N_s + k] * dphigF[k*N_N*M_G*D + sig_A*M_G*D + g*D + 0]; //Must create dphigF
		  }
		//Populate Correction
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    //Must create phigF
		    // br2_poly[e*N_N*D*N_F*N_s + sig*D*N_F*N_s + 0*N_F*N_s + field*N_s + k]
		    //phigF[k*D*N_N*_M_G + sig*_M_G + g]
		    corr_dUdx += Chi_Mult * br2_poly[om_A*N_N*D*N_F*N_s + (sig_A%N_N)*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phigF[k*D*N_N*M_G + sig_A*M_G + g];
		  }
		corr_dUdx = corr_dUdx * corr_mult;
		//corr_dUdx = 0;
		/*
		//Move correction to reference coordinates
		corr_dUdxi =  corr_dUdx * dx_dxi;
	      
		//Now, calculate the corrected gradient
		dUdxi =  naive_dUdxi  + corr_dUdxi;
	     
		dUgF[((t*N_F+field)*2 + 0)*M_G + g] = dUdxi; //the zero means element A
		*/
		naive_dUdx = naive_dUdxi * dxi_dx;
		dUdx = naive_dUdx + corr_dUdx;
		//dUdx = naive_dUdx + corr_dUdx*dxi_dx*0.01;
		/*  
		    if (ruk == 1 && field == 2)
		    {
		    printf("Side 0, node[%d]: naive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",g, naive_dUdx, corr_dUdx, dUdx);
		    }
		*/
		dUgF[t*N_F*2*M_G + field*2*M_G + 0*M_G + g] = dUdx; //the zero means element A
		//dUgF[((t*N_F+field)*2 + 0)*M_G + g] = dUdx; //the zero means element A
		//_____________________________________
		//B Cell: repopulate the dx_dxi transformation, also zero some summations
		//dx_dxi = Jac_trace(om_B*D + 0, sig_B*M_G*D + g*D + 0); //Must create this array
		dxi_dx = inv_Jac_trace[om_B*N_N*M_G*D*D + sig_B*M_G*D*D + g*D*D + 0*D + 0];
		naive_dUdxi = 0.0;
		corr_dUdx = 0.0;
		//Populate naive cell gradient
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    naive_dUdxi += U[om_B*N_F*N_s + field*N_s + k] * dphigF[k*N_N*M_G*D + sig_B*M_G*D + g*D + 0];
		  }
		//Populate Correction
		for (int k = 0; k < N_s; k = k + 1)
		  {
		    corr_dUdx += Chi_Mult * br2_poly[om_B*N_N*D*N_F*N_s + (sig_B%N_N)*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phigF[k*N_N*M_G + sig_B*M_G + g];
		  }
		corr_dUdx = corr_dUdx * corr_mult;
		//corr_dUdx = 0;
		/*
		//Move correction to reference coordinates
		corr_dUdxi =  corr_dUdx * dx_dxi;
	      
		//Now, calculate the corrected gradient
		dUdxi =  naive_dUdxi  + corr_dUdxi;
	     
		dUgF[((t*N_F+field)*2 + 1)*M_G + g] = dUdxi; //the one means element B
		*/
		naive_dUdx = naive_dUdxi * dxi_dx;
		dUdx = naive_dUdx + corr_dUdx;
		//dUdx = naive_dUdx + corr_dUdx*dxi_dx*0.01;
		/*
		  if (ruk == 1 && field == 2)
		  {
		  printf("Side 1, node[%d]: naive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",g, naive_dUdx, corr_dUdx, dUdx);
		  }
		*/
		dUgF[t*N_F*2*M_G + field*2*M_G + 1*M_G + g] = dUdx; //the 1 means element B
		//dUgF[((t*N_F+field)*2 + 1)*M_G + g] = dUdx; //the one means element B
		//_____________________________________
	      }
	  }
      }
  }
#elif TWOD
  //int D = 2;
  scalar corr_mult = 1.0;
  scalar Chi_Mult = N_N*0.5;
  int om_A;
  int sig_A;
  int om_B;
  int sig_B;

  scalar corr_dUdx;
  scalar corr_dUdy;
  scalar corr_dUdxi;
  scalar corr_dUdeta;
  scalar naive_dUdxi;
  scalar naive_dUdeta;
  scalar naive_dUdx;
  scalar naive_dUdy;
  scalar dUdxi;
  scalar dUdeta;

  scalar dxi_dx;
  scalar dxi_dy;
  scalar deta_dx;
  scalar deta_dy;
  
  for (int t = 0; t < M_T; t = t + 1) //interface
    {
      
     //The interface needs to know the bordering cells, to grab lift data
      /*
      //See dg/dg_functions.cc , dg_mappings for how to get om_A, om_B
      const simpleInterface &face = interfaces[t];   // get the interface
      const simpleElement *el0 = face.getElement(0); // get the element to the left
      const simpleElement *el1 = face.getElement(1); // get the element to the right
      om_A = ElementMap.at(el0->getId()); //cell 0
      om_B = ElementMat.at(el1->getId()); //cell 1
      sig_A = face.getClosureId(0); //interface index from cell A perspective
      sig_B = face.getClosureId(1); //interface index from cell B perspective
      */
      om_A = BR2_Map[t*4 + 0]; //cell on one side of interface
      sig_A = BR2_Map[t*4 + 1]; //reference element face linking om_A to the interface
      om_B = BR2_Map[t*4 + 2]; //cell on the other side of interface
      sig_B = BR2_Map[t*4 + 3]; //reference element face linking om_B to the interface;
      if (ruk == 3 && print_face == 1)
	{
	  printf("---For interface %d: omA = %d, omB = %d---\n",t,om_A,om_B);
	}
      for (int field = 0; field < N_F; field = field + 1)
	{
	  for (int g = 0; g < M_G; g = g + 1) //quadrature node on interface
	    {
	      //___________________________________
	      //A Cell:

	      dxi_dx  = inv_Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + g*D*D + 0*D + 0];
	      dxi_dy  = inv_Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + g*D*D + 1*D + 0];
	      deta_dx = inv_Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + g*D*D + 0*D + 1];
	      deta_dy = inv_Jac_trace[om_A*D*N_N*M_G*D*D + sig_A*M_G*D*D + g*D*D + 1*D + 1];
	      /*
	      if (field == 3 && ruk == 3 && print_face == 1)
		{
		  printf("\tSide 0 node[%d]: dxi_dx = %f, deta_dx = %f, dxi_dy = %f, deta_dy = %f\n",g,dxi_dx,deta_dx,dxi_dy,deta_dy);
		}
	      */
	      naive_dUdxi = 0.0;
	      naive_dUdeta = 0.0;
	      corr_dUdx = 0.0;
	      corr_dUdy = 0.0;
	      //Populate naive cell gradient, using DG DOF
	      for (int k = 0; k < N_s; k = k + 1)
		{
		  naive_dUdxi  += U[om_A*N_F*N_s + field*N_s + k] * dphigF[k*D*N_N*M_G*D + sig_A*M_G*D + g*D + 0];
		  naive_dUdeta += U[om_A*N_F*N_s + field*N_s + k] * dphigF[k*D*N_N*M_G*D + sig_A*M_G*D + g*D + 1];
		}
	      //Populate Correction
	      for (int k = 0; k < N_s; k = k + 1)
		{
		  //zero, one stands for x directionm y direction, respectively
		  corr_dUdx += Chi_Mult * br2_poly[om_A*N_N*D*N_F*N_s + (sig_A%N_N)*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phigF[k*D*N_N*M_G + sig_A*M_G + g];// / dxi_dx;
		  corr_dUdy += Chi_Mult * br2_poly[om_A*N_N*D*N_F*N_s + (sig_A%N_N)*D*N_F*N_s + 1*N_F*N_s + field*N_s + k] * phigF[k*D*N_N*M_G + sig_A*M_G + g];// / deta_dy;
		}

	      naive_dUdx = naive_dUdxi*dxi_dx + naive_dUdeta*deta_dx;
	      naive_dUdy = naive_dUdxi*dxi_dy + naive_dUdeta*deta_dy;

	      corr_dUdx = corr_dUdx * corr_mult;
	      corr_dUdy = corr_dUdy * corr_mult;
	      /*
	      if (ruk == 3 && corr_dUdy != 0.0)
		  {
		    printf("Nonzero dUdy Face correction encountered: t=%d,side=0,f=%d,g=%d\n", t,field, g);
		  }
	      */
	      scalar dUdx = naive_dUdx + corr_dUdx;
	      scalar dUdy = naive_dUdy + corr_dUdy;
	      if (field == 0 && ruk == 3 && print_face == 1)
		{
		  printf("\t\tomA, g=%d: naive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",g,naive_dUdx, corr_dUdx, dUdx);
		  printf("\t\tomA, g=%d: naive_dUdy = %f, corr_dUdy = %f, dUdy = %f\n",g,naive_dUdy, corr_dUdy, dUdy);
		}
	      //the 2 is for two sides per interface
	      dUgF[t*N_F*2*D*M_G + field*2*D*M_G + 0*D*M_G + 0*M_G + g] = dUdx; //the 0 means element A, 0 means x direction
	      dUgF[t*N_F*2*D*M_G + field*2*D*M_G + 0*D*M_G + 1*M_G + g] = dUdy; //the 0 means element A, 1 means y direction

	      //_____________________________________
	      //B Cell:

	      dxi_dx  = inv_Jac_trace[om_B*2*N_N*M_G*D*D + sig_B*M_G*D*D + g*D*D + 0*D + 0];
	      dxi_dy  = inv_Jac_trace[om_B*2*N_N*M_G*D*D + sig_B*M_G*D*D + g*D*D + 1*D + 0];
	      deta_dx = inv_Jac_trace[om_B*2*N_N*M_G*D*D + sig_B*M_G*D*D + g*D*D + 0*D + 1];
	      deta_dy = inv_Jac_trace[om_B*2*N_N*M_G*D*D + sig_B*M_G*D*D + g*D*D + 1*D + 1];
	      /*
	      if (field == 3 && ruk == 1 && print_face == 1)
		{
		  printf("\tSide 1 node[%d]: dxi_dx = %f, deta_dx = %f, dxi_dy = %f, deta_dy = %f\n",g,dxi_dx,deta_dx,dxi_dy,deta_dy);
		}
	      */
	      naive_dUdxi = 0.0;
	      naive_dUdeta = 0.0;
	      corr_dUdx = 0.0;
	      corr_dUdy = 0.0;
	      //Populate naive cell gradient
	      for (int k = 0; k < N_s; k = k + 1)
		{
		  naive_dUdxi +=  U[om_B*N_F*N_s + field*N_s + k] * dphigF[k*D*N_N*M_G*D + sig_B*M_G*D + g*D + 0];
		  naive_dUdeta += U[om_B*N_F*N_s + field*N_s + k] * dphigF[k*D*N_N*M_G*D + sig_B*M_G*D + g*D + 1];
		}
	      //Populate Correction
	      for (int k = 0; k < N_s; k = k + 1)
		{
		  corr_dUdx += Chi_Mult * br2_poly[om_B*N_N*D*N_F*N_s + (sig_B%N_N)*D*N_F*N_s + 0*N_F*N_s + field*N_s + k] * phigF[k*D*N_N*M_G + sig_B*M_G + g];// / dxi_dx;
		  corr_dUdy += Chi_Mult * br2_poly[om_B*N_N*D*N_F*N_s + (sig_B%N_N)*D*N_F*N_s + 1*N_F*N_s + field*N_s + k] * phigF[k*D*N_N*M_G + sig_B*M_G + g];// / deta_dy;
		}

	      naive_dUdx = naive_dUdxi*dxi_dx + naive_dUdeta*deta_dx;
	      naive_dUdy = naive_dUdxi*dxi_dy + naive_dUdeta*deta_dy;

	      corr_dUdx = corr_dUdx * corr_mult;
	      corr_dUdy = corr_dUdy * corr_mult;
	      /*
	      if (ruk == 1 && corr_dUdy != 0.0)
		{
		  printf("Nonzero dUdy Face correction encountered: t=%d,side=1,f=%d,g=%d\n", t,field, g);
		}
	      */
	      dUdx = naive_dUdx + corr_dUdx;
	      dUdy = naive_dUdy + corr_dUdy;
	     
	      if (field == 0 && ruk == 3 && print_face == 1)
		{
		  printf("\t\tomB, g=%d: naive_dUdx = %f, corr_dUdx = %f, dUdx = %f\n",g,naive_dUdx, corr_dUdx, dUdx);
		  printf("\t\tomB, g=%d: naive_dUdy = %f, corr_dUdy = %f, dUdy = %f\n",g,naive_dUdy, corr_dUdy, dUdy);
		}

	      //the 2 is for two sides per interface
	      dUgF[t*N_F*2*D*M_G + field*2*D*M_G + 1*D*M_G + 0*M_G + g] = dUdx; //the 1 means element B, 0 means x direction
	      dUgF[t*N_F*2*D*M_G + field*2*D*M_G + 1*D*M_G + 1*M_G + g] = dUdy; //the 1 means element B, 1 means y direction
	    }
	}
    }
  #endif
  //printf("Left Populate_FaceGrad\n");
}


