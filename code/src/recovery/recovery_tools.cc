/*!
  \file recovery_tools.cc
  \Recovery operations, for both full and biased recovery
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/

#include "recovery_tools.h"

void COPY_get_element_types(const int order, int &msh_hex, int &msh_qua, int &msh_tri, int &msh_lin){
  if      (order==0)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==1)  {msh_hex = MSH_HEX_8; msh_qua = MSH_QUA_4;    msh_tri = MSH_TRI_3;    msh_lin = MSH_LIN_2;  }
  else if (order==2)  {msh_hex = MSH_HEX_27; msh_qua = MSH_QUA_9;    msh_tri = MSH_TRI_6;    msh_lin = MSH_LIN_3;  }
  else if (order==3)  {msh_hex = MSH_HEX_64; msh_qua = MSH_QUA_16;   msh_tri = MSH_TRI_10;   msh_lin = MSH_LIN_4;  }
  else if (order==4)  {msh_hex = MSH_HEX_125; msh_qua = MSH_QUA_25;   msh_tri = MSH_TRI_15;   msh_lin = MSH_LIN_5;  }
  else if (order==5)  {/*msh_hex = MSH_HEX_216;*/ msh_qua = MSH_QUA_36;   msh_tri = MSH_TRI_21;   msh_lin = MSH_LIN_6;  }
  else if (order==6)  {/*msh_hex = MSH_HEX_343; */msh_qua = MSH_QUA_49;   msh_tri = MSH_TRI_28;   msh_lin = MSH_LIN_7;  }
  else if (order==7)  {/*msh_hex = MSH_HEX_512; */msh_qua = MSH_QUA_64;   msh_tri = MSH_TRI_36;   msh_lin = MSH_LIN_8;  }
  else if (order==8)  {/*msh_hex = MSH_HEX_729; */msh_qua = MSH_QUA_81;   msh_tri = MSH_TRI_45;   msh_lin = MSH_LIN_9;  }
  else if (order==9)  {/*msh_hex = MSH_HEX_1000; */msh_qua = MSH_QUA_100;  msh_tri = MSH_TRI_55;   msh_lin = MSH_LIN_10; }
  else if (order==10) {/*msh_hex = MSH_HEX_1331; */msh_qua = MSH_QUA_121;  msh_tri = MSH_TRI_66;   msh_lin = MSH_LIN_11; }
  else {printf("Invalid order number.");}
}

//First subroutine: 1D Legendre polynomials
scalar Leg_L2(int index, scalar x) //Legendre Basis on [-1,1] domian
{
  /*!
    \brief get 1D Legendre polynomial of a given index on bi-unit domain
    \param[in] index degree of Legendre polynomial
    \param[in] x location on bi-unit domain
  */
  scalar P0 = 1.0;
  switch (index)
    {
    case 0: 
      {
	return 1.0; 
	break;
      }
    default:
      {
	scalar P1 = x;
	scalar Pnew = P1;
	for (int i = 2; i <= index; i = i + 1)
	  {
	    Pnew = (2.0*i - 1.0) / (i + 0.0) * x*P1 - (i - 1.0) / (i + 0.0) * P0;
	    P0 = P1;
	    P1 = Pnew;
	  }
	return Pnew;
	break;
      }
    }
}

//Recovery basis setup, 2D simplex elements
void RecoveryIndices_2DSimplex(int N_s, int p, int N_basis, int* subIndex_relay)
{
  /*!
    \brief index the recovery/icb basis for 2D triangular(simplex) elements
    \param[in] N_s solution points per element
    \param[in] p solution polynomial order per element
    \param[in] N_basis dimension of the recovery/icb basis
    \param[out] subIndex_relay the set of pr/ps tensor product indices for recovery/icb basis
   */
  /*
    Subroutine to get indexing for Recovery on 2D simplex elements
    inputs: N_s=solution DOF per element, p=polynomial order, N_basis=dimension of recovered solution
    output: subIndex_relay
   */

  //Output of this subroutine is the subIndex structure, which is populated
  //in the deg loop above.
  int px = -1;
  int py = -1;
  for (int deg = 0; deg < N_basis; deg++)
    {
      switch (p)
	{
	case 0:
	  {
	    px = deg;
	    py = 0;
	    break;
	  }
	case 1:
	  {
	    switch(deg)
	      {
		//modes appearing in icb solution and full recovery:
	      case 0:{px=0; py=0; break;}
	      case 1:{px=0; py=1; break;}
	      case 2:{px=1; py=0; break;}

		//These go on top of Ns original modes:
	      case 3:{px=2; py=0; break;}
	      case 4:{px=1; py=1; break;}
		
		//Additional modes for full recovery:
	      
	      case 5:{px=3; py=0; break;}
		
	      default:
		{
		  printf("deg too high in Recovery_tri_index\n");
		  px = -1;
		  py = -1;
		  break;
		}
	      }
	    break;
	  }
	case 2:
	  {
	    switch(deg)
	      {
		//modes appearing in icb solution and full recovery:
	      case 0:{px=0; py=0; break;} //(0,0)
	      case 1:{px=0; py=1; break;} //(0,1)
	      case 2:{px=0; py=2; break;} //(0,2)
	      case 3:{px=1; py=0; break;} //(1,0)
	      case 4:{px=1; py=1; break;} //(1,1)
	      case 5:{px=2; py=0; break;} //(2,0)

		//These go on top of Ns original modes:
	      case 6:{px=1; py=2; break;} //(1,2)
	      case 7:{px=2; py=1; break;} //(2,1)
	      case 8:{px=3; py=0; break;} //(3,0)
		
		//Additional modes for full recovery full recovery:
	      case 9: {px=3; py=1; break;} //(3,1)
	      case 10:{px=4; py=0; break;} //(4,0)
	      case 11:{px=5; py=0; break;} //(5,0)
	      default:
		{
		  printf("deg too high in Recovery_tri_index\n");
		  px = -1;
		  py = -1;
		  break;
		}
	      }
	    break;
	  }
	case 3:
	  {
	    switch(deg)
	      {
		//Starter modes from Ns=10 DG basis functions:
	      case 0:{px=0; py=0; break;} //(0,0)
	      case 1:{px=0; py=1; break;} //(0,1)
	      case 2:{px=0; py=2; break;} //(0,2)
	      case 3:{px=0; py=3; break;} //(0,3)
	      case 4:{px=1; py=0; break;} //(1,0)
	      case 5:{px=1; py=1; break;} //(1,1)
	      case 6:{px=1; py=2; break;} //(1,2)
	      case 7:{px=2; py=0; break;} //(2,0)
	      case 8:{px=2; py=1; break;} //(2,1)
	      case 9:{px=3; py=0; break;} //(3,0)

		//The 4 ICB/Recovery modes on top of Ns DG modes:
	      case 10:{px=3; py=1; break;} //(1,3)
	      case 11:{px=2; py=2; break;} //(2,2)
	      case 12:{px=4; py=0; break;} //(3,1)
	      case 13:{px=1; py=3; break;} //(4,0)
		
		//The remaining 6 modes to complete recovery
	      case 14:{px=5; py=0; break;} //(5,0)
	      case 15:{px=6; py=0; break;} //(6,0)
	      case 16:{px=7; py=0; break;} //(7,0)
	      case 17:{px=5; py=1; break;} //(5,1)	
	      case 18:{px=4; py=1; break;} //(4,1)
	      case 19:{px=3; py=2; break;} //(3,2)
		
	  default:
	    {
	      printf("deg too high in Recovery_tri_index\n");
	      px = -1;
	      py = -1;
	      break;
	    }
	  }
	break;
      }
    case 4:
      {
	switch(deg)
	  {
	    //modes appearing in icb:
	  case 0:{px=0; py=0; break;}
	  case 1:{px=0; py=1; break;}
	  case 2:{px=0; py=2; break;}
	  case 3:{px=0; py=3; break;}
	  case 4:{px=0; py=4; break;}
	  case 5:{px=1; py=0; break;}
	  case 6:{px=1; py=1; break;}
	  case 7:{px=1; py=2; break;}
	  case 8:{px=1; py=3; break;}
	  case 9:{px=1; py=4; break;}

	  case 10:{px=2; py=0; break;}
	  case 11:{px=2; py=1; break;}
	  case 12:{px=2; py=2; break;}
	  case 13:{px=2; py=3; break;}
	  case 14:{px=3; py=0; break;}
	  case 15:{px=3; py=1; break;}
	  case 16:{px=3; py=2; break;}
	    
	    
	  case 17:{px=4; py=0; break;}
	  case 18:{px=4; py=1; break;}
	    
	  case 19:{px=5; py=0; break;}
	    
	    //modes also appearing in full recovery:
	  case 20:{px=3; py=3; break;}
	  case 21:{px=4; py=2; break;}
	  
	  case 22:{px=5; py=1; break;}
	  case 23:{px=5; py=2; break;}
	    
	  case 24:{px=6; py=0; break;}
	  case 25:{px=6; py=1; break;}
	  case 26:{px=7; py=0; break;}
	  case 27:{px=7; py=1; break;}

	  case 28:{px=8; py=0; break;}
	  case 29:{px=9; py=0; break;}

	    {
	      printf("deg too high in Recovery_tri_index\n");
	      px = -1;
	      py = -1;
	      break;
	    }
	  }
	break;
      }
	case 5: //21 DG modes, 27 icb modes, 42 recovery modes
      {
	switch(deg)
	  {
	    //modes appearing in icb:
	  case 0:{px=0; py=0; break;}
	  case 1:{px=0; py=1; break;}
	  case 2:{px=0; py=2; break;}
	  case 3:{px=0; py=3; break;}
	  case 4:{px=0; py=4; break;}
	  case 5:{px=0; py=5; break;}
	  case 6:{px=1; py=0; break;}
	  case 7:{px=1; py=1; break;}
	  case 8:{px=1; py=2; break;}
	  case 9:{px=1; py=3; break;}
	  case 10:{px=1; py=4; break;}
	  case 11:{px=1; py=5; break;}

	  case 12:{px=2; py=0; break;}
	  case 13:{px=2; py=1; break;}
	  case 14:{px=2; py=2; break;}
	  case 15:{px=2; py=3; break;}
	  case 16:{px=2; py=4; break;}
	  case 17:{px=3; py=0; break;}
	  case 18:{px=3; py=1; break;}
	  case 19:{px=3; py=2; break;}
	  case 20:{px=3; py=3; break;}


	  case 21:{px=4; py=0; break;}
	  case 22:{px=4; py=1; break;}
	  case 23:{px=4; py=2; break;}

	  case 24:{px=5; py=0; break;}
	  case 25:{px=5; py=1; break;}
	  case 26:{px=6; py=0; break;}

	    //modes also appearing in full recovery:
	  case 28:{px=3; py=4; break;}
	  
	  case 27:{px=4; py=3; break;}
	  
	  case 29:{px=5; py=2; break;}
	  case 30:{px=5; py=3; break;}


	  case 31:{px=6; py=1; break;}
	  case 32:{px=6; py=2; break;}
	  case 33:{px=7; py=0; break;}
	  case 34:{px=7; py=1; break;}
	  case 35:{px=7; py=2; break;}

	  case 36:{px=8; py=0; break;}
	  case 37:{px=8; py=1; break;}
	  case 38:{px=9; py=0; break;}
	  case 39:{px=9; py=1; break;}

	  case 40:{px=10; py=0; break;}
	  case 41:{px=11; py=0; break;}

	  default:
	    {
	      printf("deg too high in Recovery_tri_index\n");
	      px = -1;
	      py = -1;
	      break;
	    }
	  }
	break;
      }
    default:
      {
	printf("N too high for simplex recovery\n");
	px = -1;
	py = -1;
	break;
      }
    }
      subIndex_relay[deg*D + 0] = px;
      subIndex_relay[deg*D + 1] = py;
    }
}
 
//Recovery basis given a set of points
void PsiAtPoints(int N_s, int p, int N_basis, int N_N, int Npoints, scalar* x, scalar* PsiOut)
{
  /*!
    \brief Get the recovery/icb basis functions at a discrete set of points
    \param[in] N_s solution points per element
    \param[in] p Solution polynomial order per element; p=N_s-1 in 1D
    \param[in] N_basis the dimension of the recovery/icb basis (2*N_s for recovery, less for icb) 
    \param[in] N_N sides per element
    \param[in] Npoints the count of points at which the recovery basis must be populated
    \param[in] x recovery coordinates at which to populate the recovery basis
    \param[out] PsiOut the recovery basis, all modes at all points
   */

  //input: N_s = solution DOF per element
  //input: p = discretization order
  //input: N_basis = dimension of the recovery basis
  //input: Npoints = number of points where I need to populate Psi
  //input: x = block of physical coordinates: N_d entries per point
  //output: PsiOut = the Psi basis for indices 0 thru (2K-1) at all points

  //Setting up index system for non-simplex elements, because fuck it:
  //printf("Called PsiAtPoints with N_s=%d, p=%d, N_basis=%d, Npoints=%d\n",N_s,p,N_basis,Npoints);
  int subIndex[N_basis][D];
  switch (D)
    {
    case 1:
      {
	for (int n = 0; n < N_basis; n++)
	  {
	    subIndex[n][0] = n;
	  }
	break;
      }
    case 2:
      {
	if (N_N == 4)
	  {
	    for (int n = 0; n < N_basis; n++)
	      {
		subIndex[n][0] = n / (p+1);
		subIndex[n][1] = n % (p+1);
		//printf("subIndex[%d][%d]=%d\n",n,0,subIndex[n][0]);
		//printf("subIndex[%d][%d]=%d\n",n,1,subIndex[n][1]);
	      }
	  }
	else if (N_N == 3)
	  {
	    //I haven't yet discovered a general formula for this.
	    //The code is nasty, so for cleanliness in this subroutine,
	    //the indexing is outsourced
	    int* subIndex_relay = new int[N_basis*D];
	    RecoveryIndices_2DSimplex(N_s, p, N_basis, subIndex_relay);
	    for (int n = 0; n < N_basis; n++)
	      {
		for (int a = 0; a < D; a++)
		  {
		    //	    printf("subIndex_relay[%d]=%d, ",n*D+a,subIndex_relay[n*D+a]);

		    subIndex[n][a] = subIndex_relay[n*D + a];
		    //printf("subIndex[%d][%d]=%d,\n ",n,a,subIndex[n][a]);
		  }
	      }
	  }
	else 
	  {
	    printf("\n\n\nCATASTROPHE!!! D=2, but N_N is not 3 or 4\n\n\n");
	  }
	break;
      }
    case 3:
      {
	if (N_N == 6)
	  {
	    for (int n = 0; n < N_basis; n++)
	      {
		subIndex[n][0] = n / pow(p+1,2);
		subIndex[n][1] = (n/(p+1)) % (p+1);
		subIndex[n][2] = n % (p+1);
	      }
	  }
	else
	  {
	    printf("\n\n\nCATASTROPHE!!! D=3, but N_N is not 6\n\n\n");
	  }
	break;
      }
    default:
      {
	printf("Trouble in recovery_tools.cc, D out of bounds\n");
	break;
      }
    }

  //Now, use the subIndex setup to build the tensor product recovery basis
  for (int n = 0; n < N_basis; n++)
    {
      //printf("index=%d, px=%d, py=%d:\n",n, subIndex[n][0], subIndex[n][1]);
      //printf("---Psi, degree=%d:---\n",n);
      for (int j = 0; j < Npoints; j++)
	{

	  int slot = n*Npoints + j;
	  scalar product = 1.0;
	  for (int a = 0; a < D; a++)
	    {
	      //NOTE: Not sure If I'm properly calling location here
	      product = product * Leg_L2(subIndex[n][a], x[j*D + a]);
	      //      printf("j=%d, a=%d: Legendre component = %f\n", j, a, Leg_L2(subIndex[n][a], x[j*D + a]));
	    }
	  PsiOut[slot] = product;
	  //	  printf("Psi[%d][%d] = %f\n",n,j,PsiOut[slot]);
	}
    }
}

//coordinate adjustment for dealing with periodic boundaries
void PeriodicityFix(int N_s, int GQresFace, scalar* xNodes_A, scalar* xNodes_B, scalar* xFace, int* FlagPeri_out, int* sign_shift_relay, scalar* Ldomain_relay)
{
  /*!
    \brief Physical coordinate adjustment for recovery over periodic boundary interfaces
    \param[in] N_s solution points per element
    \param[in] GQresFace quadrature points per interface, same as M_G
    \param[in] xNodes_A the physical solution point locations in element A
    \param[in] xNodes_B the physical solution point locations in element B
    \param[in] xFace quadrature node phyiscal locations on shared interface
    \param[out] FlagPeri_out tells the interface what direction it is periodic in, if any
    \param[out] sign_shift_relay sign of the directional correction to link the interfaces
    \param[out] Ldomain_relay distance of periodicity correction
   */
  //Mess with physical coordinates for interfaces residing on
  //periodic-pairing interfaces
  /*
    Objective: if the interface is detected to be periodic,
    then modify physical node locations of B element
    to make recovery work properly.
    Specifically, need CompShift and LDomain out of here,
    take care of the rest in the master subroutine.

    Also output FlagPeri, which tells which direction periodicity is in for the interface.
    Also output sign_shift
   */
  //Get face centroid location
  int verbose = 0;
  if (verbose > 0)
    {
      printf("xNodes_B on entry to Periodicity Fix\n");
      for (int k = 0; k < N_s; k++)
	{
	  printf("node %d: X=(",k);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",xNodes_B[k*D+a]);
	    }
	  printf(")\n");
	}
    }
  scalar* Fc = new scalar[D];
  for (int a = 0; a < D; a++)
    {
      Fc[a] = 0.0;
    }
  for (int j = 0; j < GQresFace; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Fc[a] += xFace[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Fc[a] = Fc[a] / (GQresFace + 0.0);
    }

  //Get centroids of the two elements
  scalar* Ac = new scalar[D]; //centroid of element A
  scalar* Bc = new scalar[D]; //centroid of element B
  //This isn't a strict centroid, I'm just taking arithmetic average of the nodes
  for (int a = 0; a < D; a++)
    {
      Ac[a] = 0.0;
      Bc[a] = 0.0;
    }
  for (int j = 0; j < N_s; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Ac[a] += xNodes_A[j*D + a];
	  Bc[a] += xNodes_B[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Ac[a] = Ac[a] / (N_s + 0.0);
      Bc[a] = Bc[a] / (N_s + 0.0);
    }
 
  //Make the periodicity adjustment here:
  scalar* pathFacetoBc = new scalar[D];
  scalar* pathFacetoAc = new scalar[D];
  for (int a = 0; a < D; a++)
    {
      pathFacetoBc[a] = Bc[a] - Fc[a];
      pathFacetoAc[a] = Ac[a] - Fc[a];
      if (verbose == 1)
	{
          printf("pathFacetoAc[%d] = %f\n", a, pathFacetoAc[a]);
	  printf("pathFacetoBc[%d] = %f\n", a, pathFacetoBc[a]);
	}
    }
  //Have to specify some tolerance. I say, if centroid of B is more than
  //2.5x distance as centroid of A, we have a periodic edge.
  //But what direction do we apply periodicity in?
  int FlagPeri[D];
  for (int a = 0; a < D; a++)
    {
      FlagPeri[a] = 0;
      //Cool coding issue: If you use abs instead of fabs,
      //this routine fails to detect the periodic edges
      //for meshes with small dimensions. I'm thinking that
      //abs is maybe meant for integers and thus truncates
      //things to zero.

      //Another issue: if both the Ac and Bc paths are zero for a given component,
      //a bit of numerical pollution may trigger the case structure, so I assign
      //a minimum tolerance for Bc to trigger the case.
      scalar epsLocal = pow(10,-10);
      if (fabs(pathFacetoBc[a]) > 2.5*fabs(pathFacetoAc[a]) && fabs(pathFacetoBc[a])>epsLocal)
	{
	  if (verbose == 1){printf("Inside PeriodicityFix: Edge is periodic in direction %d\n", a);}
	  FlagPeri[a] = 1;
	}
    }
  //  printf("No segfault 5\n");
  //Whatever direction the edge is periodic in, that component of pysical coordinates must be adjusted
  //11/07/2017: For tough mesh, both FlagPeri components might be nonzero,
  //and I need to choose the more severe one. Specifically a 2D HiOCFD5 vortex issue
#ifdef TWOD
  if (FlagPeri[0] == 1 && FlagPeri[1] == 1)
    {
      //Whichever pathFacetoBc[a] magnitude is larger, that one should get FlagPeri
      scalar mag0 = fabs(pathFacetoBc[0]);
      scalar mag1 = fabs(pathFacetoBc[1]);
      FlagPeri[1] = 0;
      if (mag1 > mag0)
	{
	  FlagPeri[1] = 1;
	  FlagPeri[0] = 0;
	}
    }
#endif
  
  //The element nodes give information on problem domain: since the edge is periodic, max path
  //between element nodes will be the domain size in that direction
  scalar Ldomain = 0.0;
  int CompShift = -1;
  /*
    A cool lesson: first time, I left CompShift uniitialized, so if periodicity fix
    was not needed, there would be junk in CompShift. However, you see a few lines
    down where xNodes_A and xNodes_B are called, and Compshift is part of the index argument.
    Depending on the value of Compshift, this can give a segfault. So, ser CompShift=0 by default
   */
  for (int a = 0; a < D; a++)
    {
      if (FlagPeri[a] > 0)
	{
	  CompShift = a; 
	  for (int kA = 0; kA < N_s; kA++)
	    {
	      for (int kB = 0; kB < N_s; kB++)
		{
		  Ldomain = fmax(Ldomain, fabs(xNodes_A[kA*D + a] - xNodes_B[kB*D + a]));
		}
	    }
	}
    }
  // printf("No segfault 6\n");
  if (verbose == 1) {printf("Inside PeriodicyFix. For periodicity adjustment: Ldomain = %f. Also, CompShift(pre-correction) = %d\n", Ldomain,CompShift);}
   CompShift = fmax(0, CompShift); //to prevent segfault below
  //Now, adjust XrelB and XrelNodesB by the domain size, then recalculate a few things
  int FlagPeriSum = 0;
  for (int a = 0; a < D; a++)
    {
      FlagPeriSum += FlagPeri[a];
    }
  int sign_shift = 1; //assume that B coordinates need to be increaesd
  //Check if B coordinates need to be decreased instead:
  if (xNodes_B[0 + CompShift] > xNodes_A[0 + CompShift]) //Need to be careful about compshit to avoid overreach here
    {
      //element B is ahead of A and will need to have its coordinates decreadesd instead
      sign_shift = -1;
    }

  //Populate the output
  sign_shift_relay[0] = sign_shift;
  Ldomain_relay[0] = Ldomain;
  for (int a = 0; a < D; a++)
    {
      FlagPeri_out[a] = FlagPeri[a];
    }

  if (verbose == 1){  printf("PeriodicityFix subroutine: Preparing for pointer deletion\n");}
  delete[] Fc;
  delete[] Ac;
  delete[] Bc;
  delete[] pathFacetoBc;
  delete[] pathFacetoAc;
  
   if (verbose > 0)
    {
      printf("xNodes_B on exit from Periodicity Fix\n");
      for (int k = 0; k < N_s; k++)
	{
	  printf("node %d: X=(",k);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",xNodes_B[k*D+a]);
	    }
	  printf(")\n");
	}
    }
  if (verbose == 1) {printf("PeriodicityFix subroutine: Finished with pointer deletion\n");}


} //End periodicityFix subroutine

//Distance function, accounts for number of spatial dimenstions:
scalar distance(scalar* X1, scalar* X2)
{
  /*!
    \brief get the distance between two points
    \param[in] X1 coordinate of first point
    \param[in] X2 coordinate of second point
  */
  scalar ssq = 0.0;
  for (int a = 0; a < D; a++){
    ssq += pow(X1[a]-X2[a] , 2);}
  
  //scalar output = pow(ssq,0.5);
  scalar output = sqrt(ssq);
  return output;
}

void get_CrossProduct(scalar* U, scalar* V, scalar* perp)
{
  /*
    \brief get the cross-prodict of two vectors
    \param[in] U first vector
    \param[in] V the other vector
    \param[out] perp the perpendicular vector
  */
  int verbose = 0;
#ifdef ONED
  printf("Cross-product subroutine not ready for 1D\n");
#endif
#ifdef TWOD
  printf("Cross-product subroutine not ready for 2D\n");
#endif
#ifdef THREED
  perp[0] = U[1]*V[2] - U[2]*V[1];
  perp[1] = U[2]*V[0] - U[0]*V[2];
  perp[2] = U[0]*V[1] - U[1]*V[0];
  if (verbose == 1)
    {
      scalar mag = sqrt(perp[0]*perp[0] + perp[1]*perp[1] + perp[2]*perp[2]);
      printf("magnitude check on cross-product vector: magnitude = %f\n",mag);
    }
#endif
}

//Get recovery coordinates given physical coordinates
void TransformPhysical2Recovery(int N_s, int GQresA, int GQresFace, scalar* xNodes_A, scalar* xNodes_B, scalar* xGQ_A, scalar* xGQ_B, scalar* xFace, scalar* FaceNormal, scalar* rA_out, scalar* rB_out, scalar* rFace_out)
{
  /*!
    \brief Get the recovery coordinates (r,s) from the physical coordinates (x,y).
    GQresA can be just a few points when I'm applying biased recovery procedure
    \param[in] N_s = solution points per element, assumed same in A,B
    \param[in] GQresA = count of quadrature points in element A (assumed same in element B)
    \param[in] GQresFace = count of quadrature points on the interface, same as M_G
    \param[in] xNodes_A = gmsh nodes of element A
    \param[in] xNodes_B = gmsh nodes of element B
    \param[in] xGQ_A = physical quadrature node locations element A
    \param[in] xGQ_B = physoica; quadratre node locations element B
    \param[in] xFace = physoca; quadrature node locations on the interface
    \param[out] rA_out = recovery coordinates quadrature nodes element A
    \param[out] rB_out = recovery coordinates quadrature nodes element B
    \param[out] rFace_out = recovery coordinates quadrature nodes on face
  */
  //Establish the face centroid location
  //int GQresFace = xFace.size() / (D+0.0);
  //int GQresA = xGeo_A.size() / (D+0.0);
  // printf("No segfault 1\n");
  int verbose = 0;
  int FaceNormStyle = 0; //0 for centroid-centroid normal, 1 for actual face normal
  int FcStyle = 0; //0 for centroid of the two-element union, 1 for actual face centroid
  if (verbose > 0)
    {
  printf("Qudarature node locations imported by Tran....Recovery\n");
  printf("xGQ, element A:\n");
  for (int g = 0; g < GQresA; g++)
    {
      printf("node %d: (r1,r2,r3) = (",g);
      for (int a = 0; a < D; a++)
	{
	  printf("%f, ",xGQ_A[g*D + a]);
	}
      printf(")\n");
    }
  printf("xGQ, element B:\n");
  for (int g = 0; g < GQresA; g++)
    {
      printf("node %d: (r1,r2,r3) = (",g);
      for (int a = 0; a < D; a++)
	{
	  printf("%f, ",xGQ_B[g*D + a]);
	}
      printf(")\n");
    }
    }
  //Fc = face centroid
  scalar* Fc = new scalar[D];
  for (int a = 0; a < D; a++)
    {
      Fc[a] = 0.0;
    }
  //Periodicity fix uses xFace
  for (int j = 0; j < GQresFace; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Fc[a] += xFace[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Fc[a] = Fc[a] / (GQresFace + 0.0);
    }
  /*
  if (FcStyle == 1)
    {
      for (int j = 0; j < GQresFace; j++)
	{
	  for (int a = 0; a < D; a++)
	    {
	      Fc[a] += xFace[j*D + a];
	    }
	}
      for (int a = 0; a < D; a++)
	{
	  Fc[a] = Fc[a] / (GQresFace + 0.0);
	}
    }
  else if (FcStyle == 0)
    {
      for (int k = 0; k < 
    }
  else
    {
      printf("FcStyle out of range in TransformPhysical2Recovery\n");
      exit(1);
    }
  */
  if (verbose == 1)
    {
#ifdef THREED
      printf("Face centroid Location: (%f, %f, %f)\n", Fc[0], Fc[1], Fc[2]);
#endif
    }
  // printf("No segfault 2\n");
  //The face 'centroid' has been established. Get relative locations of all quadrature nodes
  scalar* XrelA = new scalar[GQresA * D];
  scalar* XrelB = new scalar[GQresA * D];
  scalar* XrelNodesA = new scalar[N_s * D];
  scalar* XrelNodesB = new scalar[N_s * D];
  scalar* XrelFace = new scalar[GQresFace * D];
  for (int g = 0; g < GQresA; g++)
    {
      for (int a = 0; a < D; a++)
	{
	  XrelA[g*D + a] = xGQ_A[g*D + a] - Fc[a];
	  XrelB[g*D + a] = xGQ_B[g*D + a] - Fc[a];
	}
    }
  for (int g = 0; g < GQresFace; g++)
    {
      for (int a = 0; a < D; a++)
	{
	  XrelFace[g*D + a] = xFace[g*D + a] - Fc[a];
	}
    }
  //Also, get the relative locations of all solution nodes
  for (int k = 0; k < N_s; k++)
    {
      for (int a = 0; a < D; a++)
	{
	  XrelNodesA[k*D + a] = xNodes_A[k*D + a] - Fc[a];
	  XrelNodesB[k*D + a] = xNodes_B[k*D + a] - Fc[a];
	}
    }
  //Get centroids of the two elements
  scalar* Ac = new scalar[D]; //centroid of element A
  scalar* Bc = new scalar[D]; //centroid of element B
  //This isn't a strict centroid, I'm just taking arithmetic average of the nodes
  for (int a = 0; a < D; a++)
    {
      Ac[a] = 0.0;
      Bc[a] = 0.0;
    }
  for (int j = 0; j < N_s; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Ac[a] += xNodes_A[j*D + a];
	  Bc[a] += xNodes_B[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Ac[a] = Ac[a] / (N_s + 0.0);
      Bc[a] = Bc[a] / (N_s + 0.0);
    }

  int* FlagPeri = new int[D];
  int* sign_shift_relay = new int[1];
  scalar* Ldomain_relay = new scalar[1];
  PeriodicityFix(N_s, GQresFace, xNodes_A, xNodes_B, xFace, FlagPeri, sign_shift_relay, Ldomain_relay);
  int sign_shift = sign_shift_relay[0];
  scalar Ldomain = Ldomain_relay[0];
  if (verbose > 0){printf("In TransformPhysical2Recovery. Exited Periodicity fix. sign_shift=%d, Ldomain=%f\n",sign_shift, Ldomain);
    for (int a = 0; a < D; a++){
      printf("FlagPeri[%d] = %d\n", a, FlagPeri[a]);}
  }

  //Use information from periodicity fix to adjust XrelB and XrelNodesB
  for (int a = 0; a < D; a++)
    {
      if (FlagPeri[a] > 0) //make periodicity adjustment in this direction
	{
	  for (int k = 0; k < N_s; k++)
	    {
	      XrelNodesB[k*D + a] += sign_shift*Ldomain;
	    }
	  for (int g = 0; g < GQresA; g++)
	    {
	      XrelB[g*D + a] += sign_shift*Ldomain;
	    }
	}
    }
  //From this high, the only necessary correction is to adjust centroid of B
  for (int a = 0; a < D; a++)
    {
      if (FlagPeri[a] > 0)
	{
	  Bc[a] += sign_shift*Ldomain;
	}
    }
  //printf("No segfault 7\n");
  /*
  for (int a = 0; a < D; a++)
    {
      Bc[a] = 0.0;
    }
  for (int j = 0; j < N_s; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Bc[a] += xNodes_B[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Bc[a] = Bc[a] / (N_s + 0.0);
    }
  */
  if (verbose > 0)
    {
      printf("After periodicity correction:\n");
      scalar pathFacetoBc[D];
      scalar pathFacetoAc[D];
      for (int a = 0; a < D; a++)
	{
	  pathFacetoBc[a] = Bc[a] - Fc[a];
	  pathFacetoAc[a] = Ac[a] - Fc[a];
	  printf("pathFacetoAc[%d] = %f\n", a, pathFacetoAc[a]);
	  printf("pathFacetoBc[%d] = %f\n", a, pathFacetoBc[a]);
	}
    }
  
  if (FcStyle == 0)
    {
      //Fc is centroid of the bi-element union. Changing this requires that I alter
      //everything that was just calculated
      for (int a = 0; a < D; a++)
	{
	  Fc[a] = 0.5*(Ac[a]+Bc[a]);
	  /*
	  Fc[a] = 0.0;
	  for (int g = 0; g < GQresA; g++)
	    {
	      Fc[a] += xGQ_A[g*D+a] + (xGQ_B[g*D+a] + sign_shift*Ldomain*FlagPeri[a]);
	    }
	  Fc[a] = Fc[a] / (2.0*GQresA);
	  */
	}
      //now, adjust XrelA, XrelB, XrelNodesA, and XrelNodesB and XrelFace based on different Fc
      for (int a = 0; a < D; a++)
	{
	  for (int k = 0; k < N_s; k++)
	    {
	      XrelNodesA[k*D + a] = xNodes_A[k*D + a] - Fc[a];
	      XrelNodesB[k*D + a] = (xNodes_B[k*D + a] + sign_shift*Ldomain*FlagPeri[a]) - Fc[a];
	    }
	  for (int g = 0; g < GQresA; g++)
	    {
	      XrelA[g*D + a] = xGQ_A[g*D + a] - Fc[a];
	      XrelB[g*D + a] = (xGQ_B[g*D + a] + sign_shift*Ldomain*FlagPeri[a]) - Fc[a];
	    }
	  for (int g = 0; g < GQresFace; g++)
	    {
	      XrelFace[g*D + a] = xFace[g*D + a] - Fc[a];
	    }
	}
    }

  //Now, calculate the face's approximate normal vector:
  scalar* facenorm = new scalar[D];
  if (FaceNormStyle == 0)
    {
      for (int a = 0; a < D; a++)
	{
	  facenorm[a] = Bc[a] - Ac[a];
	  //printf("facenorm[%d]=%f\n",a,facenorm[a]);
	}
      //It's not actually a normal yet: need to divide all D components by magnitude
      scalar MagNorm = 0.0;
      for (int a = 0; a < D; a++)
	{
	  MagNorm += facenorm[a]*facenorm[a];
	}
      //printf("MagNorm=%f\n",MagNorm);
      MagNorm = pow(MagNorm, 0.5);
      //printf("MagNorm=%f\n",MagNorm);
      for (int a = 0; a < D; a++)
	{
	  facenorm[a] = facenorm[a] / MagNorm;
	  if (verbose > 0) {printf("In TransformPhysical2Recovery. Centroid->Centroid facenorm[%d]=%f\n",a,facenorm[a]);}
	}
    }
  else
    {
      for (int a = 0; a < D; a++){facenorm[a] = FaceNormal[a];}
    }
  // printf("No segfault 8\n");
  //Now, I have the face's normal vector. Apply coordinate rotation.
  //See yellow notebook 05/18/2017 for the reasoning here
  scalar Mrotate[D][D];// = new scalar[D][D];
  scalar theta;
  switch (D)
    {
    case 1:
      {
	//1D case: we just have a +/-1 here for which way the interface is facing
	Mrotate[0][0] = facenorm[0];
	break;
      }
    case 2:
      {
	//here, must acount for the angle of normal.
	//theta is the counterclockwise angle of face normal, measured from positive x axis
	theta = atan2(facenorm[1], facenorm[0]); //atan2(dely, delx)
	
	Mrotate[0][0] = cos(theta);
	Mrotate[0][1] = sin(theta);
	Mrotate[1][0] = -sin(theta);
	Mrotate[1][1] = cos(theta);
	break;
      }
    case 3:
      {
	//I'm not sure that the math here is correct
	//there are two angles: theta is same as 2D case.
	//phi = counterclokwise angle of face normal in xz plane, measured from z axis positve towards z axis
	//MATRICES ARE INVERTING, BUT THE ROTATION MATRIX LOOKS FISHY TO ME. MUST DO SOME MATH
	scalar epsLocal = 2*pow(10,-9);
	//theta = atan(ny/nx)
	scalar thetaLocal = atan2(0,1); //assumption: face normal aligned with +x axis
	if (fabs(facenorm[1])+fabs(facenorm[0]) > epsLocal)
	  {
	    thetaLocal = atan2(facenorm[1], facenorm[0]);
	  }
	/*
	//phi = atan(nx/nz)
	scalar phiLocal = atan2(1,0); //assumption: face normal aligned with +x axis
	if (fabs(facenorm[0])+fabs(facenorm[2]) > epsLocal)
	  {
	    phiLocal = atan2(facenorm[0], facenorm[2]);
	  }
	printf("thetaLocal(deg) = %f, phiLocal(deg)=%f\n", thetaLocal*180.0/M_PI, phiLocal*180.0/M_PI);
	*/
	//scalar phiLocal = atan2(1,0);
	
	/*
	Mrotate[0][0] = cos(theta_local) * cos(phi_local);
	Mrotate[0][1] = sin(theta_local);
	Mrotate[0][2] = cos(theta_local) * sin(phi_local);
	Mrotate[1][0] = -sin(theta_local) * cos(phi_local);
	Mrotate[1][1] = cos(theta_local);
	Mrotate[1][2] = -sin(theta_local) * sin(phi_local);
	Mrotate[2][0] = -sin(phi_local);
	Mrotate[2][1] = 0.0;
	Mrotate[2][2] = cos(phi_local);
	*/
	/*
	Mrotate[0][0] = cos(phiLocal)*cos(thetaLocal);
	Mrotate[0][1] = cos(phiLocal)*sin(thetaLocal);
	Mrotate[0][2] = -sin(phiLocal);
	Mrotate[1][0] = -sin(thetaLocal);
	Mrotate[1][1] = cos(thetaLocal);
	Mrotate[1][2] = 0.0;
	Mrotate[2][0] = sin(phiLocal)*cos(thetaLocal);
	Mrotate[2][1] = sin(phiLocal)*sin(thetaLocal);
	Mrotate[2][2] = cos(phiLocal);
	*/
	/*
	for (int a = 0; a < D; a++){
	  for (int b = 0; b < D; b++){
	    Mrotate[a][b] = 0.0;}}

	Mrotate[0][0] = cos(thetaLocal);
	Mrotate[0][1] = sin(thetaLocal);
	Mrotate[1][0] = -sin(thetaLocal);
	Mrotate[1][1] = cos(thetaLocal);
	*/
	/*
	scalar MA[3][3];
	scalar MB[3][3];
	MA[0][0] =  cos(thetaLocal);
	MA[0][1] = sin(thetaLocal);
	MA[0][2] = 0.0;
	MA[1][0] = -sin(thetaLocal);
	MA[1][1] = cos(thetaLocal);
	MA[1][2] = 0.0;
	MA[2][0] = 0.0;
	MA[2][1] = 0.0;
	MA[2][2] = 1.0;

	//Now, identify the phi anngle in the once-rotated coordinate frame
	scalar nx_prime = facenorm[0]*(cos(thetaLocal) + sin(thetaLocal));
	scalar nz_prime = facenorm[2];
	scalar phiLocal = atan2(1,0); //assumption: face normal aligned with +xprime axis
	if (fabs(facenorm[0])+fabs(facenorm[2]) > epsLocal)
	  {
	    phiLocal = atan2(nx_prime, nz_prime);
	  }
	printf("With respect to original coordinate frame:\n");
	printf("thetaLocal(deg) = %f, phiLocal(deg)=%f\n", thetaLocal*180.0/M_PI, phiLocal*180.0/M_PI);
	printf("Instead with respect to theta-rotated coordinates: phiLocal(deg) = %f\n",phiLocal*180.0/M_PI);
	*/
	/*
	MB[0][0] = cos(phiLocal);
	MB[0][1] = 0.0;
	MB[0][2] = -sin(phiLocal);//*(-1.0);
	MB[1][0] = 0.0;
	MB[1][1] = 1.0;
	MB[1][2] = 0.0;
	MB[2][0] = sin(phiLocal);//*(-1.0);
	MB[2][1] = 0.0;
	MB[2][2] = cos(phiLocal);
	*/
	/*
	MB[0][0] = cos(phiLocal);
	MB[0][1] = 0.0;
	MB[0][2] = -sin(phiLocal);//*(-1.0);
	MB[1][0] = 0.0;
	MB[1][1] = 1.0;
	MB[1][2] = 0.0;
	MB[2][0] = sin(phiLocal);//*(-1.0);
	MB[2][1] = 0.0;
	MB[2][2] = cos(phiLocal);
	for (int row = 0; row < 3; row++)
	  {
	    for (int col = 0; col < 3; col++)
	      {
		Mrotate[row][col] = 0.0;
		for (int k = 0; k < 3; k++)
		  {
		    Mrotate[row][col] += MB[row][k]*MA[k][col];
		  }}}
	printf("MA:\n");
	for (int row = 0; row < D; row++)
	  {
	    printf("row %d: (, ",row);
	    for (int col = 0; col < D; col++)
	      {
		printf("%f, ",MA[row][col]);
	      }
	    printf(")\n");
	  }
	printf("\n");
	printf("MB:\n");
	for (int row = 0; row < D; row++)
	  {
	    printf("row %d: (, ",row);
	    for (int col = 0; col < D; col++)
	      {
		printf("%f, ",MB[row][col]);
	      }
	    printf(")\n");
	  }
	printf("\n");
	*/
	//Fuck all that, I'm taking a vector-based approach
	scalar* dir_r = new scalar[3]; //physical direction of r-coordinate
	scalar* dir_s = new scalar[3]; //physical direction of s-coordinate
	scalar* dir_t = new scalar[3]; //physical direction of t-coordinate
	//r is in the face-normal direction.
	//Also, it needs to be the actual facenormal, not the centroid-centroid path
	//across the element pair, which I foolishly label facenorm above.
	
	//	for (int a = 0; a < D; a++) {dir_r[a] = facenorm[a];}
	for (int a = 0; a < D; a++) {dir_r[a] = FaceNormal[a];}
	//s is one of face-tangential vectors: to make things pretty, use one of the edges
	//of the face.
	//But how do I know if two pointes define a face edge? I know that
	//of all the four-point combinations, 4 will yield diagonals across the face.
	//For example, might be (13,31,24, and 42). So do this: find the pair for which
	//the centroid-to-(edge midpoint) distance is maximized, and that vertex pair
	//wil definitely NOT be forming a diagonal across the face. This approach
	//would require me to know which element nodes correspond to the face.
	//So, that's going to happen here.
	int markerA[4];
	//	int markerB[4];
	for (int j = 0; j < 4; j++)
	  {
	    markerA[j] = -1;
	    //  markerB[j] = -1;
	  }
	//scalar epsFind = pow(10,-8);
	//To begin: get distances of all element A face centroid possiblilities
	scalar distShape[8*8*8*8];
	for (int j1 = 0; j1 < 8; j1++)
	  {
	    for (int j2 = 0; j2 < 8; j2++)
	      {
		for (int j3 = 0; j3 < 8; j3++)
		  {
		    for (int j4 = 0; j4 < 8; j4++)
		      {
			scalar* average = new scalar[3];
			for (int a = 0; a < D; a++)
			  {
			    average[a] = 0.25*(xNodes_A[j1*D+a] + xNodes_A[j2*D+a] + xNodes_A[j3*D+a] + xNodes_A[j4*D+a]);
			    
			  }
			//	printf("average = (%f,%f,%f)\n", average[0], average[1], average[2]);
			distShape[j1*8*8*8 + j2*8*8 + j3*8 + j4] = distance(average, Fc);
			delete[] average;
		      }
		  }
	      }
	  }
	//minimum distance yields the proper node indices
	scalar distMin = 2*distShape[0] + 100.0;
	for (int j1 = 0; j1 < 8; j1++)
	  {
	    for (int j2 = 0; j2 < 8; j2++)
	      {
		for (int j3 = 0; j3 < 8; j3++)
		  {
		    for (int j4 = 0; j4 < 8; j4++)
		      {
			if (j1 != j2 && j1 != j3 && j1 != j4 && j2 != j3 && j2 != j4 && j3 != j4)
			  { 
			    if (distShape[j1*8*8*8 + j2*8*8 + j3*8 + j4] < distMin)
			      {
				distMin = distShape[j1*8*8*8 + j2*8*8 + j3*8 + j4];

				if (verbose == 1) {
				  printf("identified distMin=%f with vertices %d, %d, %d, %d\n",distMin, j1,j2,j3,j4);
				}
				
				markerA[0] = j1;
				markerA[1] = j2;
				markerA[2] = j3;
				markerA[3] = j4;
			      }
			  }
		      }
		  }	
	      }
	  }
	if (verbose == 1)
	  {
	    printf("Vertex decision for Mrotate: omA 4 face vertices are %d, %d, %d, %d\n",markerA[0], markerA[1], markerA[2], markerA[3]);
	    printf("The physical locations of these four vertices:\n");
	    for (int j = 0; j < 4; j++)
	      {
		printf("node %d: X = (%f, %f, %f)\n", markerA[j], xNodes_A[markerA[j]*D + 0], xNodes_A[markerA[j]*D + 1], xNodes_A[markerA[j]*D + 2]);
	      }
	  }
	//Now, identify some path of markerA nodes that maximizes their midpoint's distance from face centroid
	//The thing I want to prevent is the s vector being a face diagonal
	scalar distMax = 0.0;
	int node1 = 0;
	int node2 = 3;
	for (int j1 = 0; j1 < 4; j1++)
	  {
	    for (int j2 = 0; j2 < 4; j2++) 
	      {
		if (j1 != j2)
		  {
		    scalar* average = new scalar[3];
		    for (int a = 0; a < D; a++)
		      {
			average[a] = 0.5*(xNodes_A[markerA[j1]*D+a] + xNodes_A[markerA[j2]*D+a]);
		      }
		    scalar distLocal = distance(average, Fc);
		    delete[] average;
		    if (distLocal > distMax)
		      {
			node1 = markerA[j1];
			node2 = markerA[j2];
			distMax = distLocal;
			if (verbose == 1) {printf("svector nodes are presently %d and %d, distMax = %f\n", j1,j2,distMax);}
		      }
		  }
	      }
	  }
	//Okay, node 1 and node 2 lie in the plane of the interface, and the segment between
	//them forms one of the edges of the face. So, I will now use the direction of node1->node2 to orient the s coordinate
	scalar path_12[3];
	for (int a = 0; a < D; a++)
	  {
	    path_12[a] = xNodes_A[node2*D + a] - xNodes_A[node1*D + a];
	  }
	scalar magPath = sqrt(path_12[0]*path_12[0] + path_12[1]*path_12[1] + path_12[2]*path_12[2]);
	for (int a = 0; a < D; a++)
	  {
	    dir_s[a] = path_12[a]/magPath;
	  }
	if (verbose == 1){printf("the dir_s vector that is hopefully in face plane: (%f, %f, %f)\n", dir_s[0], dir_s[1], dir_s[2]);}
	//Next up: dir_t is the cross-product of dir_r and dir_s, forming the third vector in the orthogonal set.
	get_CrossProduct(dir_r, dir_s, dir_t);

	//Now, use the three directional vectors to form the rotation matrix
	for (int a = 0; a < D; a++){
	  Mrotate[0][a] = dir_r[a];
	  Mrotate[1][a] = dir_s[a];
	  Mrotate[2][a] = dir_t[a]; }
	
	
	/*
	scalar dist2Edge[4*4];
	for (int j1 = 0; j1 < 4; j1++){
	  for (int j2 = 0; j2 < 4; j2++){
	    scalar lineCent[3];
	    for (int a = 0; a < D; a++){
	      lineCent[a] = 0.5*(
	    dist2Edge[
	break;
	*/
	delete[] dir_r;
	delete[] dir_s;
	delete[] dir_t;
	break;
      }
    }
  if (verbose > 0)
    {
      printf("In TransformPhysical2RecoveryThe Mrotate matrix for the interface:\n");
      for (int row = 0; row < D; row++)
	{
	  printf("row %d: (, ",row);
	  for (int col = 0; col < D; col++)
	    {
	      printf("%f, ",Mrotate[row][col]);
	    }
	  printf(")\n");
	}
      printf("\n");
    }
  //printf("No segfault 9\n");
  //Get the rotated, physical coordinates. This means multiply Mrotate by the Xrel coordinates
  scalar* XrelrotA = new scalar[GQresA * D];
  scalar* XrelrotB = new scalar[GQresA * D];
  scalar* XrelrotFace = new scalar[GQresFace * D];
  scalar* XrelrotNodesA = new scalar[N_s * D];
  scalar* XrelrotNodesB = new scalar[N_s * D];
  //Rotated relative coordinates of element quadrature points:
  for (int g = 0; g < GQresA; g++)
    {
      /*
      XrelrotA[g*D + 0] = Mrotate[0][0]*XrelA[g*D + 0] + Mrotate[0][1]*XrelA[g*D+1];
      XrelrotB[g*D + 0] = Mrotate[0][0]*XrelB[g*D + 0] + Mrotate[0][1]*XrelB[g*D+1];

      XrelrotA[g*D + 1] = Mrotate[1][0]*XrelA[g*D + 0] + Mrotate[1][1]*XrelA[g*D+1];
      XrelrotB[g*D + 1] = Mrotate[1][0]*XrelB[g*D + 0] + Mrotate[1][1]*XrelB[g*D+1];
      */
      
      for (int a = 0; a < D; a++)
	{
	  XrelrotA[g*D + a] = 0.0;
	  XrelrotB[g*D + a] = 0.0;
	  for (int col = 0; col < D; col++)
	    {
	      
	      XrelrotA[g*D + a] += Mrotate[a][col] * XrelA[g*D + col];
	      XrelrotB[g*D + a] += Mrotate[a][col] * XrelB[g*D + col];
	      //   printf("XrelrotA(node %d, comp %d) = %f\n",g,a,XrelrotA[g*D + a]);
	      //   printf("XrelrotB(node %d, comp %d) = %f\n",g,a,XrelrotB[g*D + a]);
	    }
	}
      
    }
  //Rotated relative coordinates of interface quadrature points:
  for (int g = 0; g < GQresFace; g++)
    {
      //XrelrotFace[g*D + 0] = Mrotate[0][0]*XrelFace[g*D + 0] + Mrotate[0][1]*XrelFace[g*D+1];
      //XrelrotFace[g*D + 1] = Mrotate[1][0]*XrelFace[g*D + 0] + Mrotate[1][1]*XrelFace[g*D+1];
      
      for (int a = 0; a < D; a++)
	{
	  XrelrotFace[g*D + a] = 0.0;
	  for (int col = 0; col < D; col++)
	    {
	      XrelrotFace[g*D + a] += Mrotate[a][col] * XrelFace[g*D + col];
	    }
	}
      
    }
  //Rotated relative coordinates of element solution nodes:
  for (int k = 0; k < N_s; k++)
    {
      for (int a = 0; a < D; a++)
	{
	  XrelrotNodesA[k*D + a] = 0.0;
	  XrelrotNodesB[k*D + a] = 0.0;
	  for (int col = 0; col < D; col++)
	    {
	      XrelrotNodesA[k*D + a] += Mrotate[a][col] * XrelNodesA[k*D + col];
	      XrelrotNodesB[k*D + a] += Mrotate[a][col] * XrelNodesB[k*D + col];
	    }
	}
    }
  //printf("No segfault 10\n");
  //First, get scaling coordinates in each direction.
  //For the scaling, I'm looking at the maximum distance of the solution nodes.
  //Why? because I want to scale recovery coordinate based on the element boundaries,
  //not just the location of the quadrature points.
  scalar* MaxRel = new scalar[D];
  for (int a = 0; a < D; a++)
    {
      MaxRel[a] = 0.0;
      for (int k = 0; k < N_s; k++)
	{
	  MaxRel[a] = fmax(MaxRel[a], fabs(XrelrotNodesA[k*D + a]));
	  MaxRel[a] = fmax(MaxRel[a], fabs(XrelrotNodesB[k*D + a]));
	}
    }
  //printf("MaxRel (r1,r2,r3) = (\n");
  for (int a = 0; a < D; a++)
    {
      //printf("%f, ",MaxRel[a]);
    }
  //printf(")\n");
  //The very last step: scale the recovery coordinates
  for (int g = 0; g < GQresA; g++)
    {
      for (int a = 0; a < D; a++)
	{
	  //printf("XrelrotA(node %d, comp %d) = %f\n",g,a,XrelrotA[g*D + a]);
	  //printf("XrelrotB(node %d, comp %d) = %f\n",g,a,XrelrotB[g*D + a]);
	  rA_out[g*D + a] = XrelrotA[g*D + a] / MaxRel[a];
	  rB_out[g*D + a] = XrelrotB[g*D + a] / MaxRel[a];
	}
    }
  //Lastly, the face coordinates:
  for (int g = 0; g < GQresFace; g++)
    {
      for (int a = 0; a < D; a++)
	{
	  rFace_out[g*D + a] = XrelrotFace[g*D + a] / MaxRel[a];
	}
    }
  // printf("No segfault 11\n");
  if (verbose > 0)
    {
      printf("The recovery coordinates:\n");
      printf("r, element A:\n");
      for (int g = 0; g < GQresA; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rA_out[g*D + a]);
	    }
	  printf(")\n");
	}
      printf("r, element B:\n");
      for (int g = 0; g < GQresA; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rB_out[g*D + a]);
	    }
	  printf(")\n");
	}
      printf("r, interface:\n");
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rFace_out[g*D + a]);
	    }
	  printf(")\n");
	}
    
      printf("theta deg= %f\n",theta*180.0/3.141592);
      printf("End of TransformPhysical2Recovery subroutine, preparing for pointer deletion\n");
  }
  //outputs are rA_out, rB_out, and rFace_out, in pointer form. This routine is finished

  //Time to delete some stuff
  
  delete[] Fc;
  delete[] XrelA;
  delete[] XrelB;
  delete[] XrelNodesA;
  delete[] XrelNodesB;
  delete[] XrelFace;
  delete[] Ac;
  delete[] Bc;
  delete[] FlagPeri;
  delete[] sign_shift_relay;
  delete[] Ldomain_relay;
  delete[] facenorm;
  delete[] XrelrotA;
  delete[] XrelrotB;
  delete[] XrelrotFace;
  delete[] XrelrotNodesA;
  delete[] XrelrotNodesB;
  delete[] MaxRel;
  
  if (verbose > 0) {printf("End of TransformPhysical2Recovery subroutine, done with pointer deletion\n");}

}

//Build an interface's discrete recovery operator
void PSIxR_shell(deck inputs, TIMERS &timers, int N_N, int N_s, int p, /*simpleMesh &m, */int M_G, int Style, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, int BoundaryTag, fullMatrix<scalar> &psi, int* Better_InvMap, int t, int M_s, int omA, int omB,scalar* PSIxR, scalar* PSIxR_A, scalar* PSIxR_B)
{
  /*!
    \brief get interface solution operator for one interface; can be recovery (Style=0), average (Style=1), or averaged biased recovery (Style=2)
    \param[in] inputs the deck information
    \param[in] N_s solution points per elements
    \param[in] p solution order 
    \param[in] m the mesh, all of it.
    \param[in] M_G quadrature points per interfaces
    \param[in] Style method indicator
    \param[in] xGeo_A solution node physical coordinates, element A
    \param[in] xGeo_B solution node physical coordinates, element B
    \param[in] xFace quadrature node physical coordinates, interface
    \param[in] FaceNormal the actual face-normal, imported from normals in main.cc
    \param[in] BoundaryTag 1 indicates that the interface is a physical boundary
    \param[in] psi Not the recovery basis, but the nodal interface basis
    \param[in] Better_InvMap holds indices of supported nodes on each face of each element
    \param[in] t the interface address
    \param[in] M_s supported nodes per face of element
    \param[in] omA element A index
    \param[in] omB element B index
    \param[in] PSIxR the interfaces (GQ_s x 2N_s) interface solution operator
    \param[out] PSIxR_A left half of interface solution operator
    \param[out] PSIxR_B right half of interface solution operator
   */
  int verbose = 0;
  int YesBlend = 0; //choice of whether or not to blend recovered solution with averaging
  scalar Blend2Avg = 0.0; //For methods {0,2}, blends Recovered solution w/ average.
  
  if (BoundaryTag == 1)
    {
      //It's a boundary interface, we just want the trace of DG solution
       if (verbose > 0)
	{
	  printf("In PSIxR_shell: This is a boundary interface, so setting PSIxR to pick up trace\n");
	}
       //It's a boundary interface. The right-hand side (element B contribution)
       //is set to zero; looking in DG kernels, you can see that both 
       //left and right halves of recovery matrix will be multiplied
       //by same element's DOF, so only one nonzero half is needed.
       //Sending full value (not half) of trace.
       //Also, there's something not exactly right with BinarySide address
       //because sometimes the boundary element is side 1 of the boundary interface.
       //So, until I fix that, the trace needs to show up from left and right halves of PSIxR_matrix
       
       fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
       for (int g = 0; g < M_G; g++){
	for (int col = 0; col < 2*N_s; col++){
	  PSIxR_matrix(g,col) = 0.0;}}
      //Now, for each g, there are M_s nonzero entries in PSIxR_matrix,
      //each corresponding to the value of a supported basis function
      for (int m = 0; m < M_s; m++)
	{
	  int fc = 0;
	  int side = 0; //omA gets side 0
	  int face = ((t*N_F+fc)*2+side)*M_s+m;
	  int Masterf0Index = Better_InvMap[face];
	  int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
	  for (int g = 0; g < M_G; g++)
	    {
	      PSIxR_matrix(g,row) = psi(g,m);
	      PSIxR_matrix(g,row+N_s) = psi(g,m);
	    }
	} 
      if (verbose > 0)
	{
	  printf("The PSixR matrix:\n");
	  for (int row = 0; row < M_G; row++)
	    {
	      printf("row [g=%d]: ",row);
	      for (int col = 0; col < 2*N_s; col++)
		{
		  printf("%f, ",PSIxR_matrix(row,col));
		}
	      printf("\n");
	    }
	}
      //Communicate PSIxR_matrix to global storage
      //Split the PSIxR matrix
      fullMatrix<scalar> PSIxR_A_matrix (M_G, N_s);
      fullMatrix<scalar> PSIxR_B_matrix (M_G, N_s);
      for (int g = 0; g < M_G; g++)
	{
	  for (int col = 0; col < N_s; col++)
	    {
	      PSIxR_A_matrix(g,col) = PSIxR_matrix(g,col + 0  );
	      PSIxR_B_matrix(g,col) = PSIxR_matrix(g,col + N_s);
	    }
	}
      
      
      //Serialize the PSIxR matrix; PSIxR, PSIxR_A, and PSIxR_B are what main.cc needs
      //The organization here must pair properly with use of PSIxR in src/dg/kernels_phil.cu
      for (int g = 0; g < M_G; g++)
	{
	  for (int col = 0; col < 2*N_s; col++)
	    {
	      int slot = g*2*N_s + col;
	      PSIxR[slot] = PSIxR_matrix(g,col);
	    }
	}
      //Serialize the split matrix components
      for (int g = 0; g < M_G; g++)
	{
	  for (int col = 0; col < N_s; col++)
	    {
	      int slot = g*N_s + col;
	      PSIxR_A[slot] = PSIxR_A_matrix(g,col);
	      PSIxR_B[slot] = PSIxR_B_matrix(g,col);
	    }
	}
    }
  else if (BoundaryTag == 0)
    {
      //Interior interface
      switch(Style)
	{
	case 0:
	  {
	    //Recovery approach: Interface solution operator is the recovery operator.
	    //This routine takes care of grabbing PSIxR, PSIxR_A, and PSIxR_B.
	    //The 0 meanns "Recovery"
	    int N_basis = 2*N_s;
	    GetPSIxR_oneface(N_N, timers, N_s, p, N_basis, M_G, 0, N_superG, phiRef, dphi, points, weight, xGeo_A, xGeo_B, xFace, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR, PSIxR_A, PSIxR_B);
	    
	    if (YesBlend == 1)
	      {
		//Use blending to modify PSIxR_A and PSIxR_B
	    //Try blending for robustness:
	    scalar* PSIxR_dummy = new scalar[M_G*2*N_s];
	    scalar* PSIxR_A_dummy = new scalar[M_G*N_s];
	    scalar* PSIxR_B_dummy = new scalar[M_G*N_s];
	
	    //The "1" means average
	    GetPSIxR_oneface(N_N, timers, N_s, p, N_basis, M_G, 1, N_superG, phiRef, dphi, points, weight, xGeo_A, xGeo_B, xFace, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR_dummy, PSIxR_A_dummy, PSIxR_B_dummy);
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < 2*N_s; col++)
		  {
		    PSIxR[g*2*N_s+col] = (1.0-Blend2Avg)*PSIxR[g*2*N_s+col] + Blend2Avg*PSIxR_dummy[g*2*N_s + col];
		  }
		for (int col = 0; col < N_s; col++)
		  {
		    PSIxR_A[g*N_s+col] = (1.0-Blend2Avg)*PSIxR_A[g*N_s+col] + Blend2Avg*PSIxR_A_dummy[g*N_s+col];
		    PSIxR_B[g*N_s+col] = (1.0-Blend2Avg)*PSIxR_B[g*N_s+col] + Blend2Avg*PSIxR_B_dummy[g*N_s+col];
		  }
	      }

	    
	    delete[] PSIxR_dummy;
	    delete[] PSIxR_A_dummy;
	    delete[] PSIxR_B_dummy;
	      }
	    break;
	  }
	case 1:
	  {
	    //BR2 approach: interface solution operator is arithmetic average
	    //This routine takes care of grabbing PSIxR, PSIxR_A, and PSIxR_B.
	    //The 1 meanns "Arithmetic average"
	    int N_basis = 2*N_s;
	    GetPSIxR_oneface(N_N, timers, N_s, p, N_basis, M_G, 1, N_superG, phiRef, dphi, points, weight, xGeo_A, xGeo_B, xFace, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR, PSIxR_A, PSIxR_B);
	    break;
	  }
	case 2:
	  {
	    //AR approach: take average of biased recoveries, unless of course we are on a boundary
	    //Step 1: Get the biased recovery operators, both A-dominant and B-dominant
	    int N_EP;  //nodal equivalence points per element (may need to change for 3D)
	    if     (inputs.getElemType() == "lin") {N_EP=1;}
	    else if(inputs.getElemType() == "tri") {N_EP=p+1;}
	    else if(inputs.getElemType() == "qua") {N_EP=p+1;}
	    else if(inputs.getElemType() == "hex") {N_EP=pow(p+1,2);}
	    int N_icb = N_EP + N_s;
	    //N_icb = N_icb - 1;
	    if (verbose > 0) {printf("In PSIxR_shell: order=%d, N_icb=%d, bout to assemble biased recovery operator for t=%d\n",p, N_icb, t);}
	    scalar* PSIxR_biA_local = new scalar[M_G*2*N_s];
	    scalar* PSIxR_biA_A_local = new scalar[M_G*N_s];
	    scalar* PSIxR_biA_B_local = new scalar[M_G*N_s];
	    scalar* PSIxR_biB_local = new scalar[M_G*2*N_s];
	    scalar* PSIxR_biB_A_local = new scalar[M_G*N_s];
	    scalar* PSIxR_biB_B_local = new scalar[M_G*N_s];
	    int order = p;
	    
	    printf("CATASTROPHE!!! PSIxR_shell is presently not ready for GradMethod == 2\n");
	    //GetPSIxR_Biased_oneface_Nov2017(inputs, N_s, order, /*m,*/ N_icb, M_G, xGeo_A, xGeo_B, xFace, FaceNormal, PSIxR_biA_local, PSIxR_biB_local, PSIxR_biA_A_local, PSIxR_biA_B_local, PSIxR_biB_A_local, PSIxR_biB_B_local);

	    //Step 2: Average the biased recovery operators to get interface solution operator
	    fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < 2*N_s; col++)
		  {
		    PSIxR_matrix(g,col) = 0.5*(PSIxR_biA_local[g*2*N_s + col] + PSIxR_biB_local[g*2*N_s + col]);
		  }
	      }

	    //Step 3b: Average that with the standard averaging operation
	    
	    // This doesn't help, even when I change to 99% average, 1% averaged recovery
	    scalar* PSIxR_dummy = new scalar[M_G*2*N_s];
	    scalar* PSIxR_A_dummy = new scalar[M_G*N_s];
	    scalar* PSIxR_B_dummy = new scalar[M_G*N_s];
	    int N_basis = 2*N_s;
	    GetPSIxR_oneface(N_N, timers, N_s, p, N_basis, M_G, 1, N_superG, phiRef, dphi, points, weight, xGeo_A, xGeo_B, xFace, FaceNormal, BoundaryTag, psi, Better_InvMap, t, M_s, omA, omB, PSIxR_dummy, PSIxR_A_dummy, PSIxR_B_dummy);
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < 2*N_s; col++)
		  {
		    PSIxR_matrix(g,col) = (1.0-Blend2Avg)*PSIxR_matrix(g,col) + Blend2Avg*PSIxR_dummy[g*2*N_s + col];
		  }
	      }
	    
	    //Step 3: Communicate back to main.cc
	    //Communicate PSIxR_matrix to global storage
	    //Split the PSIxR matrix
	    fullMatrix<scalar> PSIxR_A_matrix (M_G, N_s);
	    fullMatrix<scalar> PSIxR_B_matrix (M_G, N_s);
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < N_s; col++)
		  {
		    PSIxR_A_matrix(g,col) = PSIxR_matrix(g,col + 0  );
		    PSIxR_B_matrix(g,col) = PSIxR_matrix(g,col + N_s);
		  }
	      }
      
      
	    //Serialize the PSIxR matrix; PSIxR, PSIxR_A, and PSIxR_B are what main.cc needs
	    //The organization here must pair properly with use of PSIxR in src/dg/kernels_phil.cu
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < 2*N_s; col++)
		  {
		    int slot = g*2*N_s + col;
		    PSIxR[slot] = PSIxR_matrix(g,col);
		  }
	      }
	    //Serialize the split matrix components
	    for (int g = 0; g < M_G; g++)
	      {
		for (int col = 0; col < N_s; col++)
		  {
		    int slot = g*N_s + col;
		    PSIxR_A[slot] = PSIxR_A_matrix(g,col);
		    PSIxR_B[slot] = PSIxR_B_matrix(g,col);
		  }
	      }
      
	    //Delete some stuff
	    delete[] PSIxR_biA_local;
	    delete[] PSIxR_biA_A_local;
	    delete[] PSIxR_biA_B_local;
	    delete[] PSIxR_biB_local;
	    delete[] PSIxR_biB_A_local;
	    delete[] PSIxR_biB_B_local;

	    delete[] PSIxR_dummy;
	    delete[] PSIxR_A_dummy;
	    delete[] PSIxR_B_dummy;
	    break;
	  }
	}
    }
  else
    {
      printf("CATASTROPHE!!! Unsopported BOundaryTag value in PSIxR_shell\n");
      exit(1);
    }
}

void GetPSIxR_oneface(int N_N, TIMERS &timers, int N_s, int p, int N_basis, int M_G, int GradMethod, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight,  scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, int BoundaryTag, fullMatrix<scalar> &psi, int* Better_InvMap, int t, int M_s, int omA, int omB, scalar* PSIxR, scalar* PSIxR_A, scalar* PSIxR_B)
{
  /*!
    \brief Build the Recovery operator for a given interface, non-boundary (or if it is a boundary, it had better be a periodic boundary)
    \param[in] inputs stuff from the deck, I need mesh type and solution order
    \param[in] N_s = DOF per element
    \param[in] p = polynomial order
    \param[in] N_basis = dimension of recovery basis
    \param[in] M_G = quadrature points per interface
    \param[in] GradMethod tells us how to build recovery operator; in some options, I substitue in average instead (BR2 vs CGR)
    \param[in] xGeo_A = gmsh nodes for element A
    \param[in] xGeo_B = gmsh nodes for element B
    \param[in] xFace = quadrature node physical locations on the interface
    \param[in] FaceNormal = the interface's normal vector
    \param[in] BoundaryTag= 1 or 0, saying if the interface is a boundary interface
    \param[in] psi face basis
    \param[in] Better_InvMap tells each element which of its nodes are supported on a given interface
    \param[in] t the global interface index
    \param[in] M_s supported basis nodes per element face
    \param[in] omA the A element of the interface
    \param[in] omB the B element of the interface
    \param[out] PSIxR: the interface's discrete recovery operator
    \param[out] PSUxR_A: the half of PSIxR corresponding to element A
    \param[out] PSIxR_B: the half of PSIxR corresponding to element B
   */

  //To get physical node locations: multiply nodal DG shape functions
  //by nodal coordinate values. First step is to populate
  //the DG basis on reference element

  /*
    Musings of a grad student: First draft of this subroutine did not use the 
    det(jacobians) of the matrices, which according to all known theory,
    should be part of an integral calculation over an element. I made this decision
    for a mix of reasons. First, Marc's code assumes a constant |J| over each element.
    So, including |J| couldn't help me account for element shape. Also, I thought to 
    myself, since recovery involves integral relationships over the two elements
    separately, recovery could work properly by disregarding |J| completely and integrating
    as if each element were the reference element itself.
    
    There's just one possible hole in this reasoning: What if the two elements are different sizes?
    Would |J| allow one element to be more influential than another, as maybe it should be? Here, I am conflicted.
    I want to include |J|, but since no integral ever crosses the interface, the recovery procedure should
    be unaffected by the respective sizes of the elements; the necessary N_s constraints WILL be satisfied
    over each element whether I include |J| or not.

    This begs the question, in the case of arbitrary element geometery, can the recovered solution be
    non-unique? Perhaps I could have 2 recovered solutions that each satisfy the N_s prescribed constraints
    over each element but are not identical. I don't think this non-uniqueness is possible because the constraint
    count perfectly matches the dimension of the recovered solution (2*N_s).

    I will say, if this code is to run perturbed quad elements, I think |J| definitely needs to be
    part of the recovery process.
   */

  int verbose = 0;
 
  timers.start_timer(59);
  //Get the detJxW distribution for each element
  scalar* detJxW_A = new scalar[N_superG];
  scalar* detJxW_B = new scalar[N_superG];
  dg_detJ_OneElement(N_superG, xGeo_A, dphi, detJxW_A);
  dg_detJ_OneElement(N_superG, xGeo_B, dphi, detJxW_B);
  //That was just the detJ distribution; multiply by quadrature weights for det_JxW
  for (int g = 0; g < N_superG; g++){
    detJxW_A[g] = detJxW_A[g] * weight(g,0);
    detJxW_B[g] = detJxW_B[g] * weight(g,0); }
  //Also, this is an appropriate spot to scale the det_JxW values; they can
  //be quite small, and I don't want an overflow error on the LHS inversion
  
  scalar detNorm = 0.0;
  for (int g = 0; g < N_superG; g++)
    {
      detNorm += detJxW_A[g] + detJxW_B[g];}
  detNorm = detNorm / (2.0*N_superG);
  for (int g = 0; g < N_superG; g++){
    detJxW_A[g] = detJxW_A[g] / detNorm;
    detJxW_B[g] = detJxW_B[g] / detNorm; }
  

  if (verbose > 0)
    {
      printf("Quadrature nodes, ref element:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",points(g,a));
	    }
	  printf(")\n");
	}
      printf("Normalized detJxW, element A:\n");
      for (int g = 0; g < N_superG; g++){
	printf("\tdetJ(g=%d) = %f\tdetJxW(g=%d) = %f\n",g, detJxW_A[g]/weight(g,0), g, detJxW_A[g]); } 
      printf("Normalized detJxW, element B:\n");
      for (int g = 0; g < N_superG; g++){
	printf("\tdetJ(g=%d) = %f\tdetJxW(g=%d) = %f\n",g, detJxW_B[g]/weight(g,0), g, detJxW_B[g]); } 
    }

  //  printf("No Segfault yet A\n");
  //The solution basis has been populated on reference element: relay it to phiA, phiB storage
  fullMatrix<scalar> phi_A = phiRef;
  fullMatrix<scalar> phi_B = phiRef;

  scalar* xGQ_A = new scalar[N_superG*D];
  scalar* xGQ_B = new scalar[N_superG*D];
  //   printf("No Segfault yet B\n");
  //Got physical GQ locations in each element.
  //This procedure is particluar to a nodal basis; xGeo must be the solution node locations
  for (int g = 0; g < N_superG; g++)
    {
      for (int a = 0; a < D; a++)
	{
	  xGQ_A[g*D + a] = 0.0;
	  xGQ_B[g*D + a] = 0.0;
	  for (int k = 0; k < N_s; k++)
	    {
	      xGQ_A[g*D + a] += phi_A(g,k) * xGeo_A[k*D + a];
	      xGQ_B[g*D + a] += phi_B(g,k) * xGeo_B[k*D + a];
	    }
	}
    }

  timers.stop_timer(59);
  timers.start_timer(60);
  //  printf("No Segfault yet C\n");
  /*
  printf("Qudarature node locations in GetPSIxR_oneface\n");
  printf("xGQ, element A:\n");
  for (int g = 0; g < N_superG; g++)
    {
      printf("node %d: (r1,r2,r3) = (",g);
      for (int a = 0; a < D; a++)
	{
	  printf("%f, ",xGQ_A[g*D + a]);
	}
      printf(")\n");
    }
  printf("xGQ, element B:\n");
  for (int g = 0; g < N_superG; g++)
    {
      printf("node %d: (r1,r2,r3) = (",g);
      for (int a = 0; a < D; a++)
	{
	  printf("%f, ",xGQ_B[g*D + a]);
	}
      printf(")\n");
      }*/

  //Transform the physical quadrature locations to recovery coordinates:
  scalar* rGQ_A = new scalar[N_superG*D];
  scalar* rGQ_B = new scalar[N_superG*D];
  scalar* rFace = new scalar[M_G*D];
  // printf("No Segfault yet D\n");
  TransformPhysical2Recovery(N_s, N_superG, M_G, xGeo_A, xGeo_B, xGQ_A, xGQ_B, xFace, FaceNormal, rGQ_A, rGQ_B, rFace);
  // printf("No Segfault yet E\n");
  if (verbose > 0)
    {
      printf("Got rGQ_A, rGQ_B, and rFace from TrasformPhysical2Recovery, check return:");
      printf("The recovery coordinates:\n");
      printf("r, element A:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rGQ_A[g*D + a]);
	    }
	  printf(")\n");
	}
      printf("r, element B:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rGQ_B[g*D + a]);
	    }
	  printf(")\n");
	}
   
      printf("r, interface:\n");
      for (int g = 0; g < M_G; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rFace[g*D + a]);
	    }
	  printf(")\n");
	}
    }
  //Populate the recovery basis in both elements
  scalar* Psi_A = new scalar[N_basis*N_superG];
  scalar* Psi_B = new scalar[N_basis*N_superG];
  //printf("N_s=%d, p=%d, N_basis=%d, N_superG=%d\n",N_s,p,N_basis,N_superG);
  PsiAtPoints(N_s, p, N_basis, N_N, N_superG, rGQ_A, Psi_A);
  //return 0;
  PsiAtPoints(N_s, p, N_basis, N_N, N_superG, rGQ_B, Psi_B);
  //Also, go ahead and grab the basis along the interface
  scalar* PsiFace = new scalar[M_G*N_basis]; //recovery basis along interface quadrature points
  PsiAtPoints(N_s, p, N_basis, N_N, M_G, rFace, PsiFace);

  timers.stop_timer(60);

  /*
  printf("Psi from PsiAtPoints:\n");
  printf("---Element A:---\n");
  for (int k = 0; k < N_basis; k++)
    {
      printf("index %d:\n",k);
      for (int g = 0; g < N_superG; g++)
	{
	  printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rGQ_A[g*D+0], rGQ_A[g*D+1],Psi_A[k*N_superG + g]);
	}
    }

  printf("\n---Element B:---\n");
  for (int k = 0; k < N_basis; k++)
    {
      printf("index %d:\n",k);
      for (int g = 0; g < N_superG; g++)
	{
	  printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rGQ_B[g*D+0], rGQ_B[g*D+1],Psi_B[k*N_superG + g]);
	}
    }

  printf("\n---The Interface:---\n");
  for (int k = 0; k < N_basis; k++)
    {
      printf("index %d:\n",k);
      for (int g = 0; g < M_G; g++)
	{
	  printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rFace[g*D+0], rFace[g*D+1],PsiFace[k*M_G + g]);
	}
    }
  */

  timers.start_timer(61);
  fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
  if (BoundaryTag == 0)
    {
      //boundarytag==0 means that it is an interior interface, so use typical recovery approach.
      //Otherwise, we will substitude DG trace for the recovery operator, directly in PSIxR
      //Get the DG basis functions at the quadrature points (done above)
      //Get the quadrature weights for super-resolution recovery
      //Pretend the cell geometry jacobians don't matter, because that's how this code works.
      //10/18/2017: Adding in an arithmetic average operation here; makes recovery operator behave like trace average
      //Build recovery relationsips:
      if (GradMethod == 0)
	{
	  //the zero means that recovery should be recovery, not average
	 
	  //Build the local detJxW structure (one array for omA, one for omB)
	  //scalar det_JxW_A[N_superG];
	  //scalar det_JxW_B[N_superG];

	  //combining phi, det, and weight helps a little in terms of speed
	  fullMatrix<scalar> PHIxdetJxW_A(N_s, N_superG);
	  fullMatrix<scalar> PHIxdetJxW_B(N_s, N_superG);
	  for (int k = 0; k < N_s; k++){
	    for (int g = 0; g < N_superG; g++){
	      PHIxdetJxW_A(k,g) = phi_A(g,k)*detJxW_A[g];
	      PHIxdetJxW_B(k,g) = phi_B(g,k)*detJxW_B[g];
	    }}

	  fullMatrix<scalar> LHS (N_basis, N_basis);
	  fullMatrix<scalar> RHS (N_basis, 2*N_s);
	  //First N_s rows: equivalence over element A:
	  for (int row = 0; row < N_s; row++)
	    {
	      //LHS: DG test basis vs Recovery solution basis:
	      for (int col = 0; col < N_basis; col++)
		{
		  LHS(row,col) = 0.0;
		  for (int g = 0; g < N_superG; g++)
		    {
		      //phi call is quadrature node, then index
		      //psi call is different, because it is.
		      int slot = col*N_superG + g;
		      //LHS(row,col) += phi_A(g, row) * Psi_A[slot] * weight(g,0);
		      //LHS(row,col) += phi_A(g, row) * Psi_A[slot] * detJxW_A[g];
		      LHS(row,col) += Psi_A[slot] * PHIxdetJxW_A(row,g);
		    }
		}
	      //RHS: DG test basis vs DG solution basis:
	      for (int col = 0; col < N_s; col++)
		{
		  RHS(row,col) = 0.0;
		  for (int g = 0; g < N_superG; g++)
		    {
		      //phi call is quadrature node, then index,
		      // RHS(row,col) += phi_A(g, row) * phi_A(g, col) * weight(g,0);
		      //RHS(row,col) += phi_A(g, row) * phi_A(g, col) * detJxW_A[g];
		      RHS(row,col) += phi_A(g, col) * PHIxdetJxW_A(row,g);
		    }
		}
	    }
	  //Second N_s rows: equivalence over element B:
	  for (int row = N_s; row < N_basis; row++)
	    {
	      int testDeg = row - N_s;
	      //LHS: DG test basis vs Recovery solution basis:
	      for (int col = 0; col < N_basis; col++)
		{
		  LHS(row,col) = 0.0;
		  for (int g = 0; g < N_superG; g++)
		    {
		      //phi call is quadrature node, then index
		      //psi call is different, because it is.
		      int slot = col*N_superG + g;
		      //LHS(row,col) += phi_B(g, testDeg) * Psi_B[slot] * weight(g,0);
		      //LHS(row,col) += phi_B(g, testDeg) * Psi_B[slot] * detJxW_B[g];
		      LHS(row,col) += Psi_B[slot] * PHIxdetJxW_B(testDeg,g);
		    }
		}
	      //RHS: DG test basis vs DG solution basis:
	      for (int col = N_s; col < 2*N_s; col++)
		{
		  int solDeg = col - N_s;
		  RHS(row,col) = 0.0;
		  for (int g = 0; g < N_superG; g++)
		    {
		      //phi call is quadrature node, then index,
		      //RHS(row,col) += phi_B(g, testDeg) * phi_B(g, solDeg) * weight(g,0);
		      //RHS(row,col) += phi_B(g, testDeg) * phi_B(g, solDeg) * detJxW_B[g];
		      RHS(row,col) += phi_B(g, solDeg) * PHIxdetJxW_B(testDeg,g);
		    }
		}
	    }

	  //PEJ 11/7/2017: Also do some matrix scaling here.
	  //Each row of LHS:RHS represents a separate constraint, I'm going to
	  //scale such that the sum of each row of LHS+RHS is N_s
	  /*
	  for (int row = 0; row < N_basis; row++){
	    scalar sum_equation = 0.0;
	    for (int col = 0; col < N_basis; col++){
	      sum_equation += fabs(LHS(row,col));}
	    
	    for (int col = 0; col < 2*N_s; col++){
	      sum_equation += fabs(RHS(row,col));}
	    
	    for (int col = 0; col < N_basis; col++){
	      LHS(row,col) = LHS(row,col) * (N_s+0.0)/sum_equation;}
	    
	    for (int col = 0; col < 2*N_s; col++){
	      RHS(row,col) = RHS(row,col) * (N_s+0.0)/sum_equation;}
	    
	    if (verbose > 0)
	      {
	      printf("Scaling work, row=%d: sum_equation=%f\n",row,sum_equation);
	    }
	  }
	  */
	       
	  
	  //Now, invert the system for the Recovery matrix
	  //found this inversion command in dg/dg_functions.cc/dg_inverse_mass_matrix.
	  //Brought in french as well to match style of code
	  // Inverser la matrice de masse.
	  fullMatrix<scalar> LHS_inv(N_basis,N_basis);
	  LHS.invert(LHS_inv);
	  
	  fullMatrix<scalar> Hope_Iden(N_basis, N_basis);
	  if (verbose > 0){
	    //Only need this checking matrix if verbose>0
	    Hope_Iden.gemm(LHS,LHS_inv);}
	  
	  //Recovery matrix is LHS^-1 * RHS
	  //found the matrix multiplication command in main.cc, subroutine  LagMono2DTransformsCartesian, see last line
	  fullMatrix<scalar> MatReco(N_basis, 2*N_s);
	  MatReco.gemm(LHS_inv, RHS);

	  //Alternative hand-coded matrix multiplication:
	  //slower than the gemm call.
	  /*
	  for (int row = 0; row < N_basis; row++){
	    for (int col = 0; col < 2*N_s; col++){
	      scalar sum = 0.0;
	      for (int j = 0; j < N_basis; j++){
		sum += LHS_inv(row,j)*RHS(j,col);
	      }
	      MatReco(row, col) = sum;
	    }
	  }
	  */
	  //Now, MatReco is the interface's recovery matrix
	  
	  
	  //Now, multiply the recovery matrix by the recovery basis (at interface) for the PSIxR structure.
	  //How to use PSIxR: f(r=0) = PSIxR times (DOF_A; DOF_B)
	  
	  //fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
	  for (int g = 0; g < M_G; g++)
	    {
	      for (int col = 0; col < 2*N_s; col++)
		{
		  scalar sum = 0;
		  for (int jj = 0; jj < N_basis; jj++)
		    {
		      sum += PsiFace[jj*M_G + g] * MatReco(jj,col);
		    }
		  PSIxR_matrix(g,col) = sum;
		}
	    }
	  if (verbose > 0)
	    {
	      printf("The LHS\n");
	      for (int row = 0; row < N_basis; row++)
		{
		  printf("row %d:  (",row);
		  for (int col = 0; col < N_basis; col++)
		    {
		      printf("%6.4f, ",LHS(row,col));
		    }
		  printf("\n");
		}
	      printf("The RHS\n");
	      for (int row = 0; row < N_basis; row++)
		{
		  printf("row %d:  (",row);
		  for (int col = 0; col < 2*N_s; col++)
		    {
		      printf("%6.4f, ",RHS(row,col));
		    }
		  printf("\n");
		}
	      printf("The LHS_inverse\n");
	      for (int row = 0; row < N_basis; row++)
		{
		  printf("row %d:  (",row);
		  for (int col = 0; col < N_basis; col++)
		    {
		      printf("%6.4f, ",LHS_inv(row,col));
		    }
		  printf("\n");
		}
	      printf("LHS x LHS_inverse\n");
	      for (int row = 0; row < N_basis; row++)
		{
		  printf("row %d:  (",row);
		  for (int col = 0; col < N_basis; col++)
		    {
		      printf("%6.4f, ",Hope_Iden(row,col));
		    }
		  printf("\n");
		}
	      printf("The Recovery matrix:\n");
	      for (int row = 0; row < N_basis; row++)
		{
		  printf("row %d: ",row);
		  for (int col = 0; col < 2*N_s; col++)
		    {
		      printf("%f, ",MatReco(row,col));
		    }
		  printf("\n");
		}
	      printf("The PSixR matrix:\n");
	      for (int row = 0; row < M_G; row++)
		{
		  printf("row [g=%d]: ",row);
		  for (int col = 0; col < 2*N_s; col++)
		    {
		      printf("%f, ",PSIxR_matrix(row,col));
		    }
		  printf("\n");
		}
	    }
	}
      else if (GradMethod == 1)
	{
	  //In this case, we take the arithmetic average of the traces instead
	  for (int g = 0; g < M_G; g++){
	    for (int col = 0; col < 2*N_s; col++){
	      PSIxR_matrix(g,col) = 0.0;}}
	  //Now, for each g, there are M_s nonzero entries in PSIxR_matrix,
	  //each corresponding to the value of a supported basis function
	  for (int m = 0; m < M_s; m++)
	    {
	      int fc = 0;
	      int side = 0; //omA gets side 0
	      int face = ((t*N_F+fc)*2+side)*M_s+m;
	      int Masterf0Index = Better_InvMap[face];
	      int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
	      for (int g = 0; g < M_G; g++)
		{
		  PSIxR_matrix(g,row) = 0.5*psi(g,m);
		}
	      side = 1; //omB gets side 1
	      face = ((t*N_F+fc)*2+side)*M_s+m;
	      Masterf0Index = Better_InvMap[face];
	      row = Masterf0Index - omB*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
	      for (int g = 0; g < M_G; g++)
		{
		  PSIxR_matrix(g,row+N_s) = 0.5*psi(g,m);
		}
	      
	    }
	  if (verbose > 0)
	    {
	      printf("The PSixR matrix, GradMethod == %d:\n", GradMethod);
	      for (int row = 0; row < M_G; row++)
		{
		  printf("row [g=%d]: ",row);
		  for (int col = 0; col < 2*N_s; col++)
		    {
		      printf("%f, ",PSIxR_matrix(row,col));
		    }
		  printf("\n");
		}
	    }
	} //end case for GradMethod == 1
    } //end case for interior interface
  else
    {
      if (verbose > 0)
	{
	  printf("This is a boundary interface, so setting PSIxR to pick up trace\n");
	}
      //It's a boundary interface. The right-hand side (element B contribution)
      //is set to zero; looking in DG kernels, you can see that both 
      //left and right halves of recovery matrix will be multiplied
      //by same element's DOF, so only one nonzero half is needed.
      //Sending full value (not half) of trace.
      //Also, there's something not exactly right with BinarySide address
      //because sometimes the boundary element is side 1 of the boundary interface.
      //So, until I fix that, the trace needs to show up from left and right halves of PSIxR_matrix
      for (int g = 0; g < M_G; g++){
	for (int col = 0; col < 2*N_s; col++){
	  PSIxR_matrix(g,col) = 0.0;}}
      //Now, for each g, there are M_s nonzero entries in PSIxR_matrix,
      //each corresponding to the value of a supported basis function
      for (int m = 0; m < M_s; m++)
	{
	  int fc = 0;
	  int side = 0; //omA gets side 0
	  int face = ((t*N_F+fc)*2+side)*M_s+m;
	  int Masterf0Index = Better_InvMap[face];
	  int row = Masterf0Index - omA*(N_F*N_s) - fc*N_s; //gives k = element index corresponding to the node
	  for (int g = 0; g < M_G; g++)
	    {
	      PSIxR_matrix(g,row) = psi(g,m);
	      PSIxR_matrix(g,row+N_s) = psi(g,m);
	    }
	} 
      if (verbose > 0)
	{
	  printf("The PSixR matrix:\n");
	  for (int row = 0; row < M_G; row++)
	    {
	      printf("row [g=%d]: ",row);
	      for (int col = 0; col < 2*N_s; col++)
		{
		  printf("%f, ",PSIxR_matrix(row,col));
		}
	      printf("\n");
	    }
	}
    }
   //Split the PSIxR matrix
   fullMatrix<scalar> PSIxR_A_matrix (M_G, N_s);
   fullMatrix<scalar> PSIxR_B_matrix (M_G, N_s);
   for (int g = 0; g < M_G; g++)
     {
       for (int col = 0; col < N_s; col++)
	 {
	   PSIxR_A_matrix(g,col) = PSIxR_matrix(g,col + 0  );
	   PSIxR_B_matrix(g,col) = PSIxR_matrix(g,col + N_s);
	 }
     }


   //Serialize the PSIxR matrix
   //The organization here must pair properly with use of PSIxR in src/dg/kernels_phil.cu
   for (int g = 0; g < M_G; g++)
     {
       for (int col = 0; col < 2*N_s; col++)
	 {
	   int slot = g*2*N_s + col;
	   PSIxR[slot] = PSIxR_matrix(g,col);
	 }
     }
   //Serialize the split matrix components
   for (int g = 0; g < M_G; g++)
     {
       for (int col = 0; col < N_s; col++)
	 {
	   int slot = g*N_s + col;
	   PSIxR_A[slot] = PSIxR_A_matrix(g,col);
	   PSIxR_B[slot] = PSIxR_B_matrix(g,col);
	 }
     }
   /*
   //Now, matrix-multiply for the PSIxR structure
   for (int g = 0; g < M_G; g++)
     {
       for (int col = 0; col < 2*N_s; col++)
	 {
	   int row_of_PSIxR = g;
	   int col_of_PSIxR = col;
	   int slot = row_of_PSIxR*2*N_s + col_of_PSIxR;
	   PSIxR[slot] = 0.0;
	   for (int jj = 0; jj < N_basis; jj++)
	     {
	       PSIxR[slot] += PsiFace[row_of_PSIxR*N_basis + jj] * MatReco(jj,col_of_PSIxR);
	     }
	 }
     }
   */
   timers.stop_timer(61);
   //Done. The output from this subroutine is the PSIxR matrix
   if (verbose > 0) {printf("GetPSIxR_oneface: Preparing for pointer deletion\n");}
   delete[] xGQ_A;
   delete[] xGQ_B;
   delete[] rGQ_A;
   delete[] rGQ_B;
   delete[] rFace;
   delete[] Psi_A;
   delete[] Psi_B;
   delete[] PsiFace;
   delete[] detJxW_A;
   delete[] detJxW_B;
   if (verbose > 0) {printf("GetPSIxR_oneface: Finished with pointer deletion\n");}
}



scalar Dist2Edge_2D(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2)
{
  /*!
    \brief Distance from center to an edge
    \param x0 x-coordinate of cell center
    \param y0 y-coordinate of cell center
    \param x1 x-coordinate of first node defining the edge
    \param y1 y-coordinate of first node defining the edge
    \param x2 x-coordinate of secon node defining the edge
    \param y2 y-coordinate of second node defining the edge
    \return distance
    \section Description
    from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
  */
  return fabs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
scalar Dist2Plane_3D(scalar x1, scalar y1, scalar z1, 
		     scalar xA, scalar yA, scalar zA, 
		     scalar xB, scalar yB, scalar zB, 
		     scalar xC, scalar yC, scalar zC)
{
  /*!
    \brief get distance from a point to a plane; plane is defined by two line segemnts, AB and AC
    \param[in] xzy_1 = out-of-plane point
    \param[in] xyz_A = in-plane point
    \param[in] xyz_B = another in-plane point
    \param[in] xyz_c = another in-plane point
    Formula adapted from http://mathinsight.org/distance_point_plane
    normal vector calculation is the cross-product of AB and AC
   */

  //Step 1: Establish line segments AB and AC
  scalar AB[3]; AB[0] = xB-xA; AB[1] = yB-yA; AB[2] = zB-zA;
  scalar AC[3]; AC[0] = xC-xA; AC[1] = yC-yA; AC[2] = zC-zA;
  
  //Step 2: Take cross product of AB and AC
  scalar cross[3];
  cross[0] = AB[1]*AC[2] - AB[2]*AC[1];
  cross[1] = AB[2]*AC[0] - AB[0]*AC[2];
  cross[2] = AB[0]*AC[1] - AB[1]*AC[0];

  //Step 3: To get the normal, divide cross by its magnitude
  scalar magnitude = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
  scalar normal[3];
    for (int a = 0; a < 3; a++) {normal[a] = cross[a] / magnitude;}
  
  //Step 4: Get the path from point A to point 1
  scalar path[3];
  path[0] = x1 - xA;
  path[1] = y1 - yA;
  path[2] = z1 - zA;

  //Step 5: distance is magnitude of [dot product of (normal and path)]
  scalar distance = fabs(normal[0]*path[0] + normal[1]*path[1] + normal[2]*path[2]);
  return distance;
  
} 

void GrabPlaneVertices_3D(int GQresFace, scalar* xGeo_A, scalar* xFace, scalar* XV1, scalar* XV2, scalar* XV3, scalar* XV4)
{
  /*!brief Identify the four solution nodes that define the quadrilateral face populated by xFace
    \param[in] number of quadrature nodes on the interface
    \param[in] xGeo_A the solution nodes of element A on the interface
    \param[in] xFace the quadrature node locations on the interface
    \param[out] XV1,XV2,XV3,XV4 the solution node locations defining the interface.
  */
  int verbose = 0;
  scalar xA = xFace[0*D+0]; scalar yA = xFace[0*D+1]; scalar zA = xFace[0*D+2];
  scalar xB = xFace[1*D+0]; scalar yB = xFace[1*D+1]; scalar zB = xFace[1*D+2];
  int index_C = GQresFace-1;
  scalar xC = xFace[index_C*D+0]; scalar yC = xFace[index_C*D+1]; scalar zC = xFace[index_C*D+2];
  //Based on how the quadrature grid is laid out on each face, I assume
  //that points A,B and C are not coincident
  /*

  int LineYes = 1; //assume th three points lie on a line
  while (LineYes == 1)
    {
      scalar xAvg = 1.0/3.0 * (xA + xB + xC);
      scalar yAvg = 1.0/3.0 * (yA + yB + yC);
      scalar zAvg = 1.0/3.0 * (zA + zB + zC);
  */
  //We now seek the first four members of xGeo_A that lie in the plane.
  int bank[4];
  int count = 0;
  int k = 0;
  scalar tol = pow(10,-8);
  while (count < 4)
    {
      scalar x1 = xGeo_A[k*D+0]; scalar y1 = xGeo_A[k*D+1]; scalar z1 = xGeo_A[k*D+2];
      scalar distance = Dist2Plane_3D(x1,y1,z1 , xA,yA,zA , xB,yB,zB , xC,yC,zC);
      //if x1 is in the plane, then we have found one of the four face-defining vertices.
      if (distance < tol)
	{
	  bank[count] = k;
	  count++; //increment bank counter
	}
      k++; //search next solution point
    }
  //With the bank array populated, fill in the XV1...XV4 vectors
  for (int a = 0; a < D; a++)
    {
      XV1[a] = xGeo_A[bank[0]*D + a];
      XV2[a] = xGeo_A[bank[1]*D + a];
      XV3[a] = xGeo_A[bank[2]*D + a];
      XV4[a] = xGeo_A[bank[3]*D + a];
    }
  if (verbose > 0)
    {
      printf("Finishing GrabPlaneVertices_3D:\n");
      printf("bank = {%d, %d, %d, %d}\n", bank[0], bank[1], bank[2], bank[3]);
      for (int a = 0; a < D; a++)
    {
      printf("XV1[a=%d] = %f | XV2[a=%d] = %f | XV3[a=%d] = %f | XV4[a=%d] = %f\n", a, XV1[a], a, XV2[a], a, XV3[a], a, XV4[a]);
    }
    }
  //The coordinates of the face vertices are populated, the subroutine's work is done.
}
void Reco3DVector_From_Vertices(scalar* X1, scalar* X2, scalar* X3, scalar* X4, scalar* vr, scalar* vs, scalar* vt)
{
  /*!brief get the 3-component recovery coordinate vector from the four vertices defining the face of a hex
    \param[in] X1 physical cords of a point
    \param[in] X2 physical cords of a point
    \param[in] X3 physical cords of a point
    \param[in] X4 physical cords of a point
    \param[out] vr,vs,vt the 3-component normal+2tangent vector system
  */
  //First, get the two vectors that define the interface
  //plane. Assuming Cartesian mesh.
  int verbose = 0;
  scalar path12[D];
  scalar path14[D];
    if (verbose > 0)
      {
  printf("In Reco3DVector_from_Vertices:\n");
  for (int a = 0; a < D; a++)
    {
      printf("XV1[a=%d] = %f | XV2[a=%d] = %f | XV3[a=%d] = %f | XV4[a=%d] = %f\n", a, X1[a], a, X2[a], a, X3[a], a, X4[a]);
    }
      }
  //Grab the paths 14 and 12 first; hopefully, these are perpendicular
  for (int a = 0; a < D; a++)
    {
      path12[a] = X2[a] - X1[a];
      path14[a] = X4[a] - X1[a];
    }
  //get dot product to be sure these are perpendicular
  scalar Dot = path12[0]*path14[0] + path12[1]*path14[1] + path12[2]*path14[2];
  //If not perpendicular, adjust the longer one.
  if (fabs(Dot) > pow(10,-10))
    {
      //Whatever the longer vector is, that one is reaching across the diagonal. Must
      //Pick another point.
      scalar LocalMag12 = path12[0]*path12[0] + path12[1]*path12[1] + path12[2]*path12[2];
      scalar LocalMag14 = path14[0]*path14[0] + path14[1]*path14[1] + path14[2]*path14[2];
      if (LocalMag12 > LocalMag14)
	{
	  for (int a = 0; a < D; a++)
	    {
	      path12[a] = X3[a] - X1[a];
	    }
	}
      else
	{
	  for (int a = 0; a < D; a++)
	    {
	      path14[a] = X3[a] - X1[a];
	    }
	}
    }
  //12 and 14 are the t, s vectors. Now, get the r vector;
  //it is perpendicular to t and s, so cross product is approproate
  scalar path_r[D];
  //path12 cross path14 = det {   i           j        k     }
  //                          { path12[0] path12[1] path12[2]}
  //                          { path14[0] path14[1] path14[2]}
  path_r[0] =  path12[1]*path14[2] - path12[2]*path14[1];
  path_r[1] = -path12[0]*path14[2] + path12[2]*path14[0];
  path_r[2] =  path12[0]*path14[1] - path12[1]*path14[0];
  //Now, normalize all three vectors to get the recovery coordinate directions.
  scalar mag_r = 0.0;
  scalar mag_s = 0.0;
  scalar mag_t = 0.0;
  for (int a = 0; a < D; a++)
    {
      mag_r += path_r[a]*path_r[a];
      mag_s += path12[a]*path12[a];
      mag_t += path14[a]*path14[a];
    }
  mag_r = sqrt(mag_r);
  mag_s = sqrt(mag_s);
  mag_t = sqrt(mag_t);
  for (int a = 0; a < D; a++)
    {
      vr[a] = path_r[a] / mag_r;
      vs[a] = path12[a] / mag_s;
      vt[a] = path14[a] / mag_t;
    }
  if (verbose > 0)
    {
      printf("Leaving Reco3DVector_From_Vertices:\n");
      for (int a = 0; a < D; a++)
	{
	  printf("direction a=%d: vr=%f | vs=%f | vt=%f\n", a, vr[a], vs[a], vt[a]);
	}
    }
  //quantities of interest are vr, vs, vt. This subroutine is done.
}

//Get the recovery and reference coordinates of the nodal equivalence locations (ICB-N)
void NodalEqCord(deck inputs, int N_s, int p, int M_G, int N_EP, int N_N, const polynomialBasis* RefBasis, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, scalar* xirgA, scalar* xirgB, scalar* RrgA, scalar* RrgB, scalar* XrgA, scalar* XrgB)
{
  /*!
    \brief Get the coordinates of the non-dominant nodal equivalence points for the ICB-N method and possibly some HAG/CGR schemes that make use of biased recovery
    \param[in] inputs setup information from the deck
    \param[in] N_s solution points per element
    \param[in] p solution polynomial order
    \param[in] M_G quadrature points per interface
    \param[in] N_N neigbors (or sides) per element
    \param[in] RefBasis the DG basis and associated information on ref element
    \param[in] xGeo_A solution point locations, el A
    \param[in] xGeo_B solution point locations, el B
    \param[in] xFace quadrature point locations on the interface
    \param[in] FaceNormal the outward normal vector from the interface
    \param[out] xirgA reference coordinates, equivalence nodes, el A
    \param[out] xirgA reference coordinates, equivalence nodes, el B
    \param[out] RrgA recovery coordinates, equivalence nodes, el A
    \param[out] RrgB recovery coordinates, equivalence nodes, el B
    \param[out] XrgA physical coordinates, equivalence nodes, el A
    \param[out] XrgB physical coordinates, equivalence nodes, el B
  */
  int verbose = 0;
  //Step 1: Use Periodicity Fix subroutine to see if xGeo_B needs to be shifted
  //before this operation begins
  int* FlagPeri = new int[D];
  int* sign_shift_relay = new int[1];
  scalar* Ldomain_relay = new scalar[1];
  int GQresFace = M_G;
  
  PeriodicityFix(N_s, GQresFace, xGeo_A, xGeo_B, xFace, FlagPeri, sign_shift_relay, Ldomain_relay);
  int sign_shift = sign_shift_relay[0];
  scalar Ldomain = Ldomain_relay[0];
  if (verbose > 0)
    {
      printf("Inside NodalEqCord, Exited Periodicity fix. sign_shift=%d, Ldomain=%f\n",sign_shift, Ldomain);
      printf("The FlagPeri structure:\n");
      for (int a = 0; a < D; a++)
	{
	  printf("FlagPeri[%d]=%d\n", a, FlagPeri[a]);
	}
    }
  //delete[] sign_shift_relay;
  //delete[] Ldomain_relay;
#ifdef ONED
  {
    //1) for each element, identify the solution point number
    //corresponding to the interface
    scalar eps = pow(10,-10);
    int markerA[1]; markerA[0] = -1;
    int markerB[1]; markerB[0] = -1;
    scalar* FaceCent = new scalar[D];
    FaceCent[0] = xFace[0];
    //Execute periodicity fux on xGeo_B; if the interface
    //is interior, then correction will be set to zero
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, commencing periodicity fix\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] += sign_shift*Ldomain;
	      }
	  }
      }

    //Need REFERENCE coordinates of the face vertex
    //coordinates in each element. Massive pain.

    //Find interface node index in element A:
    for (int k = 0; k < N_s; k++)
      {
	scalar dist = fabs(FaceCent[0] - xGeo_A[k*D + 0]);
	if (dist < eps)
	  {
	    //the node corresponds to the interface location
	    markerA[0] = k;
	    if (verbose > 0) {printf("kA=%d, dist=%f\n", k, dist);}
	  }
      }
    //Find interface node in element B:
    for (int k = 0; k < N_s; k++)
      {
	scalar dist = fabs(FaceCent[0] - xGeo_B[k*D + 0]);
	if (dist < eps)
	  {
	    //the node corresponds to the interface location
	    markerB[0] = k;
	    if (verbose > 0) {printf("kB=%d, dist=%f\n", k, dist);}
	  }
      }
    if (verbose > 0)
      {
	printf("FaceCent(phys) = %f\n", FaceCent[0]);
	printf("markerA=%d, markerB=%d\n",markerA[0],markerB[0]);
      }


    //Given the solution node number, I can get its reference coordinates
    fullMatrix<double> RefPoints = RefBasis->points;
    scalar xiFace_A[1*D];
    scalar xiFace_B[1*D];
    xiFace_A[0*D + 0] = RefPoints(markerA[0], 0);

    xiFace_B[0*D + 0] = RefPoints(markerB[0], 0);

    if (verbose > 0) {printf("xiFace_A=%f, xiFace_B=%f\n", xiFace_A[0], xiFace_B[0]);}

    //Identify the p+1 gauss quadrature location on reference element
    fullMatrix<double> GQpoints, weight;
    gaussIntegration::getLine(p*2+1 , GQpoints, weight);

    if (verbose > 0)
      {
	printf("The %d quadrature point collection for choosing nodal equivalence:\n",int(GQpoints.size1()));
	for (int g = 0; g < GQpoints.size1(); g++)
	  {
	    printf("xi_gq(%d) = %f\n", g, GQpoints(g,0));
	  }
      }

    //Next: Find the N_EP quadrature nodes closest to the interface. These will be
    //the equivalence nodes.
    //Use the distance_to_edge subroutine from simplemesh.cc for this task in 2D.
    scalar EqDistA[weight.size1()]; 
    scalar EqDistB[weight.size1()];
    for (int g = 0; g < weight.size1(); g++)
      {
	//This is all done wrt reference coordinates:
	EqDistA[g] = fabs(GQpoints(g,0) - xiFace_A[0]);
	EqDistB[g] = fabs(GQpoints(g,0) - xiFace_B[0]);
      }

    //The sorting algorithm to identify the closest N_EP gaussian quadrature nodes
    int sortA[weight.size1()];
    int sortB[weight.size1()];
    for (int j = 0; j < weight.size1(); j++)
      {
	sortA[j] = weight.size1()-1;
	sortB[j] = weight.size1()-1;
      }
    for (int g1 = 0; g1 < weight.size1(); g1++)
      {
	for (int g2 = 0; g2 < weight.size1(); g2++)
	  {
	    if (EqDistA[g1] < EqDistA[g2] && g1 != g2)
	      {
		sortA[g1] -= 1;
	      }
	    if (EqDistB[g1] < EqDistB[g2] && g1 != g2)
	      {
		sortB[g1] -= 1;
	      }
	  }
      }

    int tiebreak = 1;
    while(tiebreak == 1)
      {
        for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 < g2)
		  {
		    sortA[g1] -= 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 < g2)
		  {
		    sortB[g1] -= 1;
		  }
	      }
	  }
	tiebreak = 0;
	for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
	      }
	  }
      }
    
    int FinalMapA[N_EP];
    int FinalMapB[N_EP];
    for (int j = 0; j < N_EP; j++)
      {
	for (int g = 0; g < weight.size1(); g++)
	  {
	    if (sortA[g] == j)
	      {
		FinalMapA[j] = g;
	      }
	    if (sortB[g] == j)
	      {
		FinalMapB[j] = g;
	      }
	  }
      }

    for (int j = 0; j < N_EP; j++)
      {
	for (int a = 0; a < D; a++)
	  {
	    xirgA[j*D + a] = GQpoints(FinalMapA[j], a);
	    xirgB[j*D + a] = GQpoints(FinalMapB[j], a);
	  }
      }
  
    //Okay, I now have xi coordinates of equivalence nodes. Get physical coordinates
    int GQres = weight.size1();
    fullMatrix<scalar> phi (GQres,N_s); 
    fullMatrix<double> phiD (GQres,N_s); 
    //RefBasis->f (RefPoints, phiD);
    RefBasis->f (GQpoints, phiD);
    for(int g = 0; g < GQres; g++){
      for(int i = 0; i < N_s; i++){
	phi(g,i) = (scalar)phiD(g,i);
      }
    }    

    //physical coordinates of equivalence nodes:
    //scalar* XrgA = new scalar[N_EP*D];
    //scalar* XrgB = new scalar[N_EP*D];
    

    //use basis functions to get the physical locations of equivalence nodes
    for (int g = 0; g < N_EP; g++){
      for (int a = 0; a < D; a++){
	XrgA[g*D + a] = 0.0;
	XrgB[g*D + a] = 0.0;
	for (int k = 0; k < N_s; k++)
	  {
	    //Basis functions are multiplied by physical solution point coordinates
	    //to get physical equivalence node coordinates.
	    XrgA[g*D + a] += phi(FinalMapA[g], k) * xGeo_A[k*D + a];
	    XrgB[g*D + a] += phi(FinalMapB[g], k) * xGeo_B[k*D + a];
	  } //end k summation
      }} //end equivalence node and dimension loop

    /*
    printf("Physical coordinates, equivalence nodes\n");
    printf("omA:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, XrgA[j*D + 0], XrgA[j*D + 1]);
      }
    printf("omB:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, XrgB[j*D + 0], XrgB[j*D + 1]);
      }
    */

    scalar* RrgFace = new scalar[D*M_G];
    //Now get the recovery coordinates corresponding to xirgA and xirgB
    //TransformPhysical2Recovery(N_s, GQres, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, RrgA, RrgB, RrgFace);
    TransformPhysical2Recovery(N_s, N_EP, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, FaceNormal, RrgA, RrgB, RrgFace);
    //TransformPhysical2Recovery(N_s, N_EP, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, RrgA, RrgB, RrgFace);
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, removing periodicity fix\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] -= sign_shift*Ldomain;
	      }
	  }
      }
    if (verbose > 0) {printf("End of 1D NodalEqCord after calling Recovery coordinates\n");}

    delete[] FlagPeri;
    delete[] sign_shift_relay;
    delete[] Ldomain_relay;
    delete[] FaceCent;
    //delete[] XrgA;
    //delete[] XrgB;
    delete[] RrgFace;

    //xirgA, xirgB, RrgA, RrgB have been populated. This subroutine's work is done

  } //end ONED case
#endif
#ifdef TWOD
  {
    //Execute periodicity fix on xGeo_B; if the interface
    //is interior, then correction will be set to zero
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, commencing periodicity fix on xGeo_B\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] += sign_shift*Ldomain;
	      }
	  }
      }

    //1) for each element, identify the solution point numbers
    //corresponding to the interface
    scalar eps = pow(10,-10);
    
    //marker indicates the 2 solution nodes of each element corresponding to interface
    int markerA[2];
    int markerB[2];
    for (int j = 0; j < 2; j++){
      markerA[j] = 0;
      markerB[j] = 0; }
    
    //Need REFERENCE coordinates of the face vertex
    //coordinates in each element. Massive pain.

    //Face centroid:
    scalar* FaceCent = new scalar[D];
    for (int a = 0; a < D; a++){
      FaceCent[a] = 0.0; }
    
    for (int g = 0; g < GQresFace; g++){
      for (int a = 0; a < D; a++){
	FaceCent[a] += 1.0/(GQresFace+0.0) * xFace[g*D+a]; }}
    
    //Now, for a given element, pick any two vertices
    // (nodes 0,1,2,(3) in gmsh ordering). The combination
    //that yields a centroid closest to FaceCent is the pair of
    //solution points corresponding to the face.
    //First, take care of element A
    scalar dist_local[N_N*N_N];
    for (int j1 = 0; j1 < N_N; j1++) 
      {
	for (int j2 = 0; j2 < N_N; j2++) 
	  {
	    //midpoint between the two solution points:
	    scalar* midpoint = new scalar[D];
	    for (int a = 0; a < D; a++)
	      {
		midpoint[a] = 0.5*(xGeo_A[j1*D + a] + xGeo_A[j2*D + a]);
	      }
	    dist_local[j1*N_N + j2] = distance(FaceCent, midpoint);
	    //   printf("omA: dist_local(%d,%d)=%f\n", j1,j2, dist_local[j1*N_N + j2]);
	    delete[] midpoint;
	  }
      }
    //Whichever dist_local is the smallest, that is the pair of vertices associated with the face
    scalar dist_min = dist_local[0];
    //   printf("omA: original dist_min=%f\n", dist_min);
    for (int j1 = 0; j1 < N_N; j1++)
      {
	for (int j2 = 0; j2 < N_N; j2++)
	  {
	    //use leq so that markerA gets set properly of zeroth interface is winner
	    if (dist_local[j1*N_N + j2] <= dist_min)
	      {
		dist_min = dist_local[j1*N_N + j2];
		markerA[0] = j1;
		markerA[1] = j2;
		//	printf("Set dist_min=%f, markerA[0]=%d, markerA[1]=%d\n", dist_min, markerA[0], markerA[1]);
	      }
	  }
      }
    //Now, markerA holds the indices of the face vertices from omA.
    //Repeat for omB:
    for (int j1 = 0; j1 < N_N; j1++) 
      {
	for (int j2 = 0; j2 < N_N; j2++) 
	  {
	    //midpoint between the two solution points:
	    scalar* midpoint = new scalar[D];
	    for (int a = 0; a < D; a++)
	      {
		midpoint[a] = 0.5*(xGeo_B[j1*D + a] + xGeo_B[j2*D + a]);
	      }
	    dist_local[j1*N_N + j2] = distance(FaceCent, midpoint);
	    //  printf("omB: dist_local(%d,%d)=%f\n", j1,j2, dist_local[j1*N_N + j2]);
	    delete[] midpoint;
	  }
      }
    //Whichever dist_local is the smallest, that is the pair of vertices associated with the face
    dist_min = dist_local[0];
    //   printf("omB: original dist_min=%f\n", dist_min);
    for (int j1 = 0; j1 < N_N; j1++)
      {
	for (int j2 = 0; j2 < N_N; j2++)
	  {
	    //use leq so that markerB gets set properly of zeroth interface is winner
	    if (dist_local[j1*N_N + j2] <= dist_min)
	      {
		dist_min = dist_local[j1*N_N + j2];
		markerB[0] = j1;
		markerB[1] = j2;
		//	printf("Set dist_min=%f, markerB[0]=%d, markerB[1]=%d\n", dist_min, markerB[0], markerB[1]);
	      }
	  }
      }
    //Now, markerB holds the indices of the face vertices from omB.
    if (verbose > 0)
      {
	printf("Idetified vertex nodes in both elements.\n");
	printf("omA: vertex indices are %d and %d\n", markerA[0], markerA[1]);
	printf("omB: vertex indices are %d and %d\n", markerB[0], markerB[1]);
      }

    //Given A solution node number, I can get its reference coordinates
    fullMatrix<double> RefPoints = RefBasis->points;
    //const fullMatrix<double> &nodes;// = m.getNodes(); //I think this is ref cords of nodes
    scalar xiFace_A[2*D];
    scalar xiFace_B[2*D];
    xiFace_A[0*D + 0] = RefPoints(markerA[0], 0);
    xiFace_A[0*D + 1] = RefPoints(markerA[0], 1);
    xiFace_A[1*D + 0] = RefPoints(markerA[1], 0);
    xiFace_A[1*D + 1] = RefPoints(markerA[1], 1);

    xiFace_B[0*D + 0] = RefPoints(markerB[0], 0);
    xiFace_B[0*D + 1] = RefPoints(markerB[0], 1);
    xiFace_B[1*D + 0] = RefPoints(markerB[1], 0);
    xiFace_B[1*D + 1] = RefPoints(markerB[1], 1);

    if (verbose > 0)
      {
	printf("Reference coordinates for the interface vertices:\n");
	printf("omA:\n");
	printf("first node: xi=%f, eta=%f\t || \t second node: xi=%f, eta=%f\n", xiFace_A[0*D+0], xiFace_A[0*D+1], xiFace_A[1*D+0], xiFace_A[1*D+1]);
	printf("omB:\n");
	printf("first node: xi=%f, eta=%f\t || \t second node: xi=%f, eta=%f\n", xiFace_B[0*D+0], xiFace_B[0*D+1], xiFace_B[1*D+0], xiFace_B[1*D+1]);
      }

    //printf("The Gaussian quadrature locations from GaussLegendre Routine\n");
    fullMatrix<double> GQpoints, weight;
    if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(p*2+1, GQpoints, weight);
    //if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(p*2, GQpoints, weight);
    if(inputs.getElemType() == "qua") gaussIntegration::getQuad(p*2+1, GQpoints, weight);
    if (verbose > 0)
      {
	printf("All of the Gaussian quadrature points, ref coordinates:\n");
	for (int g = 0; g < weight.size1(); g++)
	  {
	    printf("g=%d: xi=%f, eta=%f\n", g, GQpoints(g,0), GQpoints(g,1));
	  }
      }
    //Next: Find the N_EP quadrature nodes closest to the interface. These will be
    //the equivalence nodes.
    //Use the distance_to_edge subroutine from simplemesh.cc for this task
    scalar EqDistA[weight.size1()]; 
    scalar EqDistB[weight.size1()];
    for (int g = 0; g < weight.size1(); g++)
      {
	//This is all done wrt reference coordinates:
	//distance to edge(xNode, yNode, x face vertex, y face vertex, x face vertex, y face vertex
	EqDistA[g] = Dist2Edge_2D(GQpoints(g,0), GQpoints(g,1), xiFace_A[0*D+0], xiFace_A[0*D+1], xiFace_A[1*D+0], xiFace_A[1*D+1]);
        EqDistB[g] = Dist2Edge_2D(GQpoints(g,0), GQpoints(g,1), xiFace_B[0*D+0], xiFace_B[0*D+1], xiFace_B[1*D+0], xiFace_B[1*D+1]);
	if (EqDistA[g] < pow(10,-11) || EqDistB[g] < pow(10,-11))
	  {
	    printf("ALERT: In Nodal Eq Cord subroutine: EqDistA[%d]=%f, EqDistB[%d]=%f\n",g, EqDistA[g], g, EqDistB[g]); fflush(stdout);
	    exit(1);
	  }
      }
    /*
    printf("Distances from face to quadrature nodes, reference coordinates:\n");
    printf("omA:\n");
    for (int g = 0; g < weight.size1(); g++){
      printf("g=%d: distance=%f\n", g, EqDistA[g]); }
    printf("omB:\n");
    for (int g = 0; g < weight.size1(); g++){
      printf("g=%d: distance=%f\n", g, EqDistB[g]); }
    */
 
    int sortA[weight.size1()];
    int sortB[weight.size1()];
    for (int j = 0; j < weight.size1(); j++)
      {
	sortA[j] = weight.size1()-1;
	sortB[j] = weight.size1()-1;
      }
    for (int g1 = 0; g1 < weight.size1(); g1++)
      {
	for (int g2 = 0; g2 < weight.size1(); g2++)
	  {
	    if (EqDistA[g1] < EqDistA[g2] && g1 != g2)
	      {
		sortA[g1] -= 1;
	      }
	    if (EqDistB[g1] < EqDistB[g2] && g1 != g2)
	      {
		sortB[g1] -= 1;
	      }
	  }
      }

    /*
    printf("Element A sort, pre-tiebraker:\n");
    for (int g =0; g < weight.size1(); g++)
      {
	printf("sortA[node %d] = place %d\n", g, sortA[g]);
      }
    printf("Element B sort, pre-tiebraker:\n");
    for (int g =0; g < weight.size1(); g++)
      {
	printf("sortB[node %d] = place %d\n", g, sortB[g]);
      }
    */

    //Some sortA entries are identical because some nodes are same distance from interface.
    //So, use an arbitrary tiebracker to get a one-to-one mapping

    //this one worked early for simplex, but not quad:
    /*
    for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 != g2)
		  {
		    sortA[g1] -= 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 != g2)
		  {
		    sortB[g1] -= 1;
		  }
	      }
	  }
    */

    /*
    for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 < g2)
		  {
		    sortA[g1] -= 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 < g2)
		  {
		    sortB[g1] -= 1;
		  }
	      }
	  }
*/

    
    int tiebreak = 1;
    while(tiebreak == 1)
      {
        for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 < g2)
		  {
		    sortA[g1] -= 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 < g2)
		  {
		    sortB[g1] -= 1;
		  }
	      }
	  }
	tiebreak = 0;
	for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
	      }
	  }
      }
    /*
    printf("Element A sort, post-tiebraker:\n");
    for (int g =0; g < weight.size1(); g++)
      {
	printf("sortA[node %d] = place %d\n", g, sortA[g]);
      }
    printf("Element B sort, post-tiebraker:\n");
    for (int g =0; g < weight.size1(); g++)
      {
	printf("sortB[node %d] = place %d\n", g, sortB[g]);
      }
    */
    int FinalMapA[N_EP];
    int FinalMapB[N_EP];
    for (int j = 0; j < N_EP; j++)
      {
	for (int g = 0; g < weight.size1(); g++)
	  {
	    if (sortA[g] == j)
	      {
		FinalMapA[j] = g;
	      }
	    if (sortB[g] == j)
	      {
		FinalMapB[j] = g;
	      }
	  }
      }

    /*
    printf("Final Map A\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("j=%d, corresponding node=%d\n", j, FinalMapA[j]);
      }
    printf("Final Map B\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("j=%d, corresponding node=%d\n", j, FinalMapB[j]);
      }
    */

    //Get the closest N_EP nodes in each element, store the reference coordinates
    //and also get the recovery coordinates
    for (int j = 0; j < N_EP; j++)
      {
	for (int a = 0; a < D; a++)
	  {
	    xirgA[j*D + a] = GQpoints(FinalMapA[j], a);
	    xirgB[j*D + a] = GQpoints(FinalMapB[j], a);
	  }
      }

    /*
    printf("Reference coordinates, equivalence nodes\n");
    printf("omA:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, xirgA[j*D + 0], xirgA[j*D + 1]);
      }
    printf("omB:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, xirgB[j*D + 0], xirgB[j*D + 1]);
      }
    */

    //Okay, I now have xi coordinates of equivalence nodes. Get physical coordinates
    int GQres = weight.size1();
    fullMatrix<scalar> phi (GQres,N_s); 
    fullMatrix<double> phiD (GQres,N_s); 
    //RefBasis->f (RefPoints, phiD);
    RefBasis->f (GQpoints, phiD);
    for(int g = 0; g < GQres; g++){
      for(int i = 0; i < N_s; i++){
	phi(g,i) = (scalar)phiD(g,i);
      }
    }    

    //physical coordinates of equivalence nodes:
    //scalar* XrgA = new scalar[N_EP*D];
    //scalar* XrgB = new scalar[N_EP*D];
    

    //use basis functions to get the physical locations of equivalence nodes
    for (int g = 0; g < N_EP; g++){
      for (int a = 0; a < D; a++){
	XrgA[g*D + a] = 0.0;
	XrgB[g*D + a] = 0.0;
	for (int k = 0; k < N_s; k++)
	  {
	    //Basis functions are multiplied by physical solution point coordinates
	    //to get physical equivalence node coordinates.
	    XrgA[g*D + a] += phi(FinalMapA[g], k) * xGeo_A[k*D + a];
	    XrgB[g*D + a] += phi(FinalMapB[g], k) * xGeo_B[k*D + a];
	  } //end k summation
      }} //end equivalence node and dimension loop

    /*
    printf("Physical coordinates, equivalence nodes\n");
    printf("omA:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, XrgA[j*D + 0], XrgA[j*D + 1]);
      }
    printf("omB:\n");
    for (int j = 0; j < N_EP; j++)
      {
	printf("node %d: (x,y) = (%f,%f)\n", j, XrgB[j*D + 0], XrgB[j*D + 1]);
      }
    */

    scalar* RrgFace = new scalar[D*M_G];
    //Now get the recovery coordinates corresponding to xirgA and xirgB
    //TransformPhysical2Recovery(N_s, GQres, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, RrgA, RrgB, RrgFace);
    TransformPhysical2Recovery(N_s, N_EP, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, FaceNormal, RrgA, RrgB, RrgFace);

    //One more thing: Remove the periodicity correction placed on xGeo_B; these
    //coordinates are an input from main.cc and should not be altered
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, removing periodicity fix on xGeo_B at end of Nodal Eq Cord subroutine\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] -= sign_shift*Ldomain;
	      }
	  }
      }


    if (verbose > 0) {printf("End of NodalEqCord after calling Recovery coordinates\n");}
    //printf("FlagPeri Size is %d\n", int(FlagPeri.size1()));
    if (verbose > 0) {printf("Check 1\n");}
    fflush(stdout);
    delete[] RrgFace;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.7\n");}
    delete[] FlagPeri;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.3\n");}
    delete[] FaceCent;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.4\n");}
    //delete[] XrgA;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.5\n");}
    //delete[] XrgB;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.6\n");}
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.1\n");}
    delete[] sign_shift_relay;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.2\n");}
    delete[] Ldomain_relay;
    fflush(stdout);
    
    
    //return 0;
    /*
    //Execute periodicity fux on xGeo_B; if the interface
    //is interior, then correction will be set to zero
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    printf("FlagPeri[%d] nonzero, commencing periodicity fix\n",a);
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] += sign_shift*Ldomain;
	      }
	  }
      }
    //Find interface node in element A:
    for (int k = 0; k < N_s; k++)
      {
	scalar dist = fabs(FaceCent - xGeo_A[k*D + 0]);
	if (dist < eps)
	  {
	    //the node corresponds to the interface location
	    markerA = k;
	    printf("kA=%d, dist=%f\n", k, dist);
	  }
      }
    //Find interface node in element B:
    for (int k = 0; k < N_s; k++)
      {
	scalar dist = fabs(FaceCent - xGeo_B[k*D + 0]);
	if (dist < eps)
	  {
	    //the node corresponds to the interface location
	    markerB = k;
	    printf("kB=%d, dist=%f\n", k, dist);
	  }
      }
    printf("FaceCent(phys) = %f\n", FaceCent);
  
    printf("markerA=%d, markerB=%d\n",markerA,markerB);

    //if either markerA or markerB is still 0, edge is periodic
    //I think it's always marker B because the interface is attached primarily to omA
  


    //Given the solution node number, I can get its reference coordinates
    const fullMatrix<double> &nodes = m.getNodes(); //I think this is ref cords of nodes
    scalar xiFace_A = nodes(0, markerA); //xi_1 coordinate of interface node
    scalar xiFace_B = nodes(0, markerB);

    printf("xiFace_A=%f, xiFace_B=%f\n", xiFace_A, xiFace_B);

    //Identify the p+1 gauss quadrature location on reference element
    fullMatrix<double> points, weight;
    gaussIntegration::getLine(fmax((p+1)*2,1), points, weight);

    printf("The %d quadrature point collection for choosing nodal equivalence:\n",int(points.size1()));
    for (int g = 0; g < points.size1(); g++)
      {
	printf("xi_gq(%d) = %f\n", g, points(g,0));
      }

    //Get the standard polynomial basis for the reference element:
    int msh_qua;
    int msh_tri;
    int msh_lin;
    COPY_get_element_types(p, msh_qua, msh_tri, msh_lin);
    int elem_type = msh_lin; //DIMENSION_DEPENDENT~ Do not forget.
    const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
  
    //Identify the xi_quadrature location closest to marker node for each element
    xirgA[0] = points(0,0);
    xirgB[0] = points(0,0);
    scalar dist_winner_A = fabs(xiFace_A - xirgA[0]);
    scalar dist_winner_B = fabs(xiFace_B - xirgB[0]);
    int index_winner_A = 0;
    int index_winner_B = 0;
    int GQres = points.size1();
    //For A element:
    for (int g = 0; g < points.size1(); g++)
      {
	scalar dist = fabs(xiFace_A - points(g,0));
	if (dist < dist_winner_A)
	  {
	    xirgA[0] = points(g,0);
	    dist_winner_A = dist;
	    index_winner_A = g;
	  }
      }
    //For B element:
    for (int g = 0; g < points.size1(); g++)
      {
	scalar dist = fabs(xiFace_B - points(g,0));
	if (dist < dist_winner_B)
	  {
	    xirgB[0] = points(g,0);
	    dist_winner_B = dist;
	    index_winner_B = g;
	  }
      }
 
    //Okay, I now have xi coordinates of equivalence nodes. Get physical coordinates
    fullMatrix<scalar> phi (GQres,N_s); 
    fullMatrix<double> phiD (GQres,N_s); 
    basis->f (points, phiD);
    for(int g = 0; g < GQres; g++){
      for(int i = 0; i < N_s; i++){
	phi(g,i) = (scalar)phiD(g,i);
      }
    }    

    //physical coordinates of equivalence nodes:
    scalar* XrgA = new scalar[N_EP*D];
    scalar* XrgB = new scalar[N_EP*D];
    
    //use basis functions to get the physical locations of equivalence nodes
    for (int g = 0; g < N_EP; g++){
      for (int a = 0; a < D; a++){
	XrgA[g*D + a] = 0.0;
	XrgB[g*D + a] = 0.0;
	for (int k = 0; k < N_s; k++)
	  {
	    //index_winner is a specific volume quadrature node associated
	    //with a given equivalence point, so there are N_EP of them per element.
	    //Note, index_winner is not basis function index (that's k).
	    //Basis functions are multiplied by physical solution point coordinates
	    //to get physical equivalence node coordinates.
	    XrgA[g*D + a] += phi(index_winner_A[g], k) * xGeo_A[k*D + a];
	    XrgB[g*D + a] += phi(index_winner_B[g], k) * xGeo_B[k*D + a];
	  } //end k summation
      }} //end equivalence node and dimension loop
    
    scalar* RrgFace = new scalar[D*M_G];
    //get recovery coordinates from physical coordinates
    TransformPhysical2Recovery(N_s, GQres, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, RrgA, RrgB, RrgFace);
    //xirgA, xirgB, RrgA, RrgB have been populated. This subroutine's work is done
    */  
    if (verbose > 0) {printf("Check 2\n");}
  } //end TWOD case 
#endif
  if (verbose > 0) {printf("Check 2.5\n");}
#ifdef THREED
  {
    if (N_N == 6)
      {
	//Do nothing
      }
    else
      {
	printf("CATASTROPHE!!! in NodalEqCord: D=3, but N_N != 6\n");
      }
    //Execute periodicity fix on xGeo_B; if the interface
    //is interior, then correction will be set to zero
    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, commencing periodicity fix on xGeo_B\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] += sign_shift*Ldomain;
	      }
	  }
      }

    //1) for each element, identify the solution point numbers
    //corresponding to the interface
    scalar eps = pow(10,-10);
    //marker indicates the 4 solution nodes of each element corresponding to interface
    int markerA[4];
    int markerB[4];
    for (int j = 0; j < 4; j++){
      markerA[j] = 0;
      markerB[j] = 0; }
    //Need REFERENCE coordinates of the face vertex
    //coordinates in each element. Massive pain.

    //Face centroid:
    scalar* FaceCent = new scalar[D];
    for (int a = 0; a < D; a++){
      FaceCent[a] = 0.0; }
    
    for (int g = 0; g < GQresFace; g++){
      for (int a = 0; a < D; a++){
	FaceCent[a] += 1.0/(GQresFace+0.0) * xFace[g*D+a]; }}
    
    //Now, for a given element, pick any two cube vertices
    // (nodes 0,1,2,3,4,5,6,7 in gmsh ordering). The combination
    //that yields a centroid closest to FaceCent is the pair of
    //solution points corresponding to the face. Also, I think that in 3D
    //it is important to place a check that all j indices are different when setting marker.
    //First, take care of element A
    scalar dist_local[8*8*8*8]; //8 corner vertices crossed with themselves
    for (int j1 = 0; j1 < 8; j1++) 
      {
	for (int j2 = 0; j2 < 8; j2++) 
	  {
	    for (int j3 = 0; j3 < 8; j3++)
	      {
		for (int j4 = 0; j4 < 8; j4++)
		  {
		    //centroid of the four solution points:
		    scalar* midpoint = new scalar[D];
		    for (int a = 0; a < D; a++)
		      {
			midpoint[a] = 0.25*(xGeo_A[j1*D + a] + xGeo_A[j2*D + a] + xGeo_A[j3*D + a] + xGeo_A[j4*D + a]);
		      }
		    dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4] = distance(FaceCent, midpoint);
		    
		    //	    printf("omA: dist_local(%d,%d,%d,%d)=%f\n", j1,j2,j3,j4, dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4]);
		    delete[] midpoint;
		  }
	      }
	  }
      }
    //Whichever dist_local is the smallest, that is the pair of vertices associated with the face
    scalar dist_min = dist_local[0];
    if (verbose > 0) {printf("omA: original dist_min=%f\n", dist_min);}
    for (int j1 = 0; j1 < N_N; j1++)
      {
	for (int j2 = 0; j2 < N_N; j2++)
	  {
	    for (int j3 = 0; j3 < 8; j3++)
	      {
		for (int j4 = 0; j4 < 8; j4++)
		  {
		    //use leq so that markerA gets set properly of zeroth interface is winner
		    if (dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4] <= dist_min)
		      {
			if (j1 != j2 && j1 != j3 && j1 != j4 && j2 != j3 && j2 != j4 && j3 != j4)
			  {
			    dist_min = dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4];
			    markerA[0] = j1;
			    markerA[1] = j2;
			    markerA[2] = j3;
			    markerA[3] = j4;
			    //	printf("Set dist_min=%f, markerA[0]=%d, markerA[1]=%d\n", dist_min, markerA[0], markerA[1]);
			  }
		      }
		  }
	      }
	  }
      }
    if (verbose > 0) {printf("After running through possible face vertices for omA: dist_min = %f\n", dist_min);} 
    //Now, repeat that entire process for omB:
    for (int j1 = 0; j1 < 8; j1++) 
      {
	for (int j2 = 0; j2 < 8; j2++) 
	  {
	    for (int j3 = 0; j3 < 8; j3++)
	      {
		for (int j4 = 0; j4 < 8; j4++)
		  {
		    //centroid of the four solution points:
		    scalar* midpoint = new scalar[D];
		    for (int a = 0; a < D; a++)
		      {
			midpoint[a] = 0.25*(xGeo_B[j1*D + a] + xGeo_B[j2*D + a] + xGeo_B[j3*D + a] + xGeo_B[j4*D + a]);
		      }
		    dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4] = distance(FaceCent, midpoint);
		    //	    printf("omA: dist_local(%d,%d,%d,%d)=%f\n", j1,j2,j3,j4, dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4]);
		    delete[] midpoint;
		  }
	      }
	  }
      }
    //Whichever dist_local is the smallest, that is the pair of vertices associated with the face
    dist_min = dist_local[0];
    if (verbose > 0) {printf("omB: original dist_min=%f\n", dist_min);}
    for (int j1 = 0; j1 < N_N; j1++)
      {
	for (int j2 = 0; j2 < N_N; j2++)
	  {
	    for (int j3 = 0; j3 < 8; j3++)
	      {
		for (int j4 = 0; j4 < 8; j4++)
		  {
		    //use leq so that markerA gets set properly of zeroth interface is winner
		    if (dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4] <= dist_min)
		      {
			if (j1 != j2 && j1 != j3 && j1 != j4 && j2 != j3 && j2 != j4 && j3 != j4)
			  {
			    dist_min = dist_local[j1*(8*8*8) + j2*(8*8) + j3*8 + j4];
			    markerB[0] = j1;
			    markerB[1] = j2;
			    markerB[2] = j3;
			    markerB[3] = j4;
			    //	printf("Set dist_min=%f, markerA[0]=%d, markerA[1]=%d\n", dist_min, markerA[0], markerA[1]);
			  }
		      }
		  }
	      }
	  }
      }
    if (verbose > 0)
      {
	printf("After running through possible face vertices for omB: dist_min = %f\n", dist_min); 
	//Now, markerB holds the indices of the face vertices from omB.
	printf("Idetified vertex nodes in both elements.\n");
	printf("omA: vertex indices are %d, %d, %d, and %d\n", markerA[0], markerA[1], markerA[2], markerA[3]);
	printf("omB: vertex indices are %d, %d, %d, and %d\n", markerB[0], markerB[1], markerB[2], markerB[3]);
      }
    //Given A solution node number, I can get its reference coordinates
    fullMatrix<double> RefPoints = RefBasis->points;
    //const fullMatrix<double> &nodes;// = m.getNodes(); //I think this is ref cords of nodes
    //4 corners per interface:
    scalar xiFace_A[4*D];
    scalar xiFace_B[4*D];
    xiFace_A[0*D + 0] = RefPoints(markerA[0], 0);
    xiFace_A[0*D + 1] = RefPoints(markerA[0], 1);
    xiFace_A[0*D + 2] = RefPoints(markerA[0], 2);////

    xiFace_A[1*D + 0] = RefPoints(markerA[1], 0);
    xiFace_A[1*D + 1] = RefPoints(markerA[1], 1);
    xiFace_A[1*D + 2] = RefPoints(markerA[1], 2);////

    xiFace_A[2*D + 0] = RefPoints(markerA[2], 0);
    xiFace_A[2*D + 1] = RefPoints(markerA[2], 1);
    xiFace_A[2*D + 2] = RefPoints(markerA[2], 2);////

    xiFace_A[3*D + 0] = RefPoints(markerA[3], 0);
    xiFace_A[3*D + 1] = RefPoints(markerA[3], 1);
    xiFace_A[3*D + 2] = RefPoints(markerA[3], 2);////


    xiFace_B[0*D + 0] = RefPoints(markerB[0], 0);
    xiFace_B[0*D + 1] = RefPoints(markerB[0], 1);
    xiFace_B[0*D + 2] = RefPoints(markerB[0], 2);////

    xiFace_B[1*D + 0] = RefPoints(markerB[1], 0);
    xiFace_B[1*D + 1] = RefPoints(markerB[1], 1);
    xiFace_B[1*D + 2] = RefPoints(markerB[1], 2);////

    xiFace_B[2*D + 0] = RefPoints(markerB[2], 0);
    xiFace_B[2*D + 1] = RefPoints(markerB[2], 1);
    xiFace_B[2*D + 2] = RefPoints(markerB[2], 2);////

    xiFace_B[3*D + 0] = RefPoints(markerB[3], 0);
    xiFace_B[3*D + 1] = RefPoints(markerB[3], 1);
    xiFace_B[3*D + 2] = RefPoints(markerB[3], 2);////
    
    fullMatrix<double> GQpoints, weight;
    if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(p*2 + 1, GQpoints, weight);

    //Next: Find the N_EP quadrature nodes closest to the interface. These will be
    //the equivalence nodes.
    //Use the distance_to_face subroutine from simplemesh.cc for this task
    scalar EqDistA[weight.size1()]; 
    scalar EqDistB[weight.size1()];
    for (int g = 0; g < weight.size1(); g++)
      {
	//This is all done wrt reference coordinates:
	//distance to edge(out-of-plane point, 3*In-plane points)
	EqDistA[g] = Dist2Plane_3D(GQpoints(g,0), GQpoints(g,1), GQpoints(g,2), 
     				   xiFace_A[0*D + 0], xiFace_A[0*D + 1], xiFace_A[0*D + 2],
				   xiFace_A[1*D + 0], xiFace_A[1*D + 1], xiFace_A[1*D + 2],
				   xiFace_A[2*D + 0], xiFace_A[2*D + 1], xiFace_A[2*D + 2]);
	//printf("g=%d: distance=%f\n", g, EqDistA[g]); 
        EqDistB[g] = Dist2Plane_3D(GQpoints(g,0), GQpoints(g,1), GQpoints(g,2), 
     				   xiFace_B[0*D + 0], xiFace_B[0*D + 1], xiFace_B[0*D + 2],
				   xiFace_B[1*D + 0], xiFace_B[1*D + 1], xiFace_B[1*D + 2],
				   xiFace_B[2*D + 0], xiFace_B[2*D + 1], xiFace_B[2*D + 2]);
      }
    if (verbose == 1)
      {
	printf("Minimal Quadrature points in REFERENCE coordinates\n");
	for (int g = 0; g < weight.size1(); g++)
	  {
	    printf("node %d: (xi,eta,zeta) = (%f, %f, %f)\n", g, GQpoints(g,0), GQpoints(g,1), GQpoints(g,2));
	  }
	printf("omA: The face vertices, reference coordinates\n");
	for (int j = 0; j < 4; j++)
	  {
	    printf("vertex %d: (xi,eta,zeta) = (%f, %f, %f)\n",j, xiFace_A[j*D+0], xiFace_A[j*D+1], xiFace_A[j*D+2]);
	  }
	printf("omB: The face vertices, reference coordinates\n");
	for (int j = 0; j < 4; j++)
	  {
	    printf("vertex %d: (xi,eta,zeta) = (%f, %f, %f)\n",j, xiFace_B[j*D+0], xiFace_B[j*D+1], xiFace_B[j*D+2]);
	  }
	    
	printf("Distances from face to quadrature nodes, ELEMENT REFERENCE coordinates:\n");
	printf("omA:\n");
	for (int g = 0; g < weight.size1(); g++){
	  printf("g=%d: distance=%f\n", g, EqDistA[g]); }
	printf("omB:\n");
	for (int g = 0; g < weight.size1(); g++){
	  printf("g=%d: distance=%f\n", g, EqDistB[g]); }
      }
    //Now, need to sort the quadrature nodes according to who is closest to the interface.
    //This is a two-step process, requires a tiebreaking step.
    int sortA[weight.size1()];
    int sortB[weight.size1()];
    for (int j = 0; j < weight.size1(); j++)
      {
	sortA[j] = weight.size1()-1;
	sortB[j] = weight.size1()-1;
      }
    for (int g1 = 0; g1 < weight.size1(); g1++)
      {
	for (int g2 = 0; g2 < weight.size1(); g2++)
	  {
	    if (EqDistA[g1] < EqDistA[g2] && g1 != g2)
	      {
		sortA[g1] -= 1;
	      }
	    if (EqDistB[g1] < EqDistB[g2] && g1 != g2)
	      {
		sortB[g1] -= 1;
	      }
	  }
      }
    if (verbose == 1)
      {
	printf("Element A sort, pre-tiebraker:\n");
	for (int g =0; g < weight.size1(); g++)
	  {
	    printf("sortA[node %d] = place %d\n", g, sortA[g]);
	  }
	printf("Element B sort, pre-tiebraker:\n");
	for (int g =0; g < weight.size1(); g++)
	  {
	    printf("sortB[node %d] = place %d\n", g, sortB[g]);
	  }
      }

    //This next step assigns each node a unique spot in sorting order
    int tiebreak = 1;
    while(tiebreak == 1)
      {
        for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 < g2)
		  {
		    sortA[g1] -= 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 < g2)
		  {
		    sortB[g1] -= 1;
		  }
	      }
	  }
	tiebreak = 0;
	for (int g1 = 0; g1 < weight.size1(); g1++)
	  {
	    for (int g2 = 0; g2 < weight.size1(); g2++)
	      {
		if (sortA[g1] == sortA[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
		if (sortB[g1] == sortB[g2] && g1 != g2)
		  {
		    tiebreak = 1;
		  }
	      }
	  }
      }
    if (verbose == 1)
      {
	printf("Element A sort, post-tiebraker:\n");
	for (int g =0; g < weight.size1(); g++)
	  {
	    printf("sortA[node %d] = place %d\n", g, sortA[g]);
	  }
	printf("Element B sort, post-tiebraker:\n");
	for (int g =0; g < weight.size1(); g++)
	  {
	    printf("sortB[node %d] = place %d\n", g, sortB[g]);
	  }
      }
    
    int FinalMapA[N_EP];
    int FinalMapB[N_EP];
    for (int j = 0; j < N_EP; j++)
      {
	for (int g = 0; g < weight.size1(); g++)
	  {
	    if (sortA[g] == j)
	      {
		FinalMapA[j] = g;
	      }
	    if (sortB[g] == j)
	      {
		FinalMapB[j] = g;
	      }
	  }
      }

    if (verbose == 1)
      {
	printf("Final Map A\n");
	for (int j = 0; j < N_EP; j++)
	  {
	    printf("j=%d, corresponding quadrature node=%d\n", j, FinalMapA[j]);
	  }
	printf("Final Map B\n");
	for (int j = 0; j < N_EP; j++)
	  {
	    printf("j=%d, corresponding quadrature node=%d\n", j, FinalMapB[j]);
	  }
      }

    //Get the closest N_EP nodes in each element, store the reference coordinates
    //and also get the recovery coordinates
    for (int j = 0; j < N_EP; j++)
      {
	for (int a = 0; a < D; a++)
	  {
	    xirgA[j*D + a] = GQpoints(FinalMapA[j], a);
	    xirgB[j*D + a] = GQpoints(FinalMapB[j], a);
	  }
      }

    if (verbose == 1)
      {
	printf("Reference coordinates, equivalence nodes\n");
	printf("omA:\n");
	for (int j = 0; j < N_EP; j++)
	  {
	    printf("node %d: (x,y,z) = (%f,%f,%f)\n", j, xirgA[j*D+0], xirgA[j*D+1], xirgA[j*D+2]);
	  }
	printf("omB:\n");
	for (int j = 0; j < N_EP; j++)
	  {
	    printf("node %d: (x,y,z) = (%f,%f,%f)\n", j, xirgB[j*D+0], xirgB[j*D+1], xirgB[j*D+2]);
	  }
      }

    //Okay, I now have xi coordinates of equivalence nodes. Get physical coordinates using nodal basis
    int GQres = weight.size1();
    fullMatrix<scalar> phi (GQres,N_s); 
    fullMatrix<double> phiD (GQres,N_s); 
    //RefBasis->f (RefPoints, phiD);
    RefBasis->f (GQpoints, phiD);
    for(int g = 0; g < GQres; g++){
      for(int i = 0; i < N_s; i++){
	phi(g,i) = (scalar)phiD(g,i);
      }
    }    

    //physical coordinates of equivalence nodes:
    //scalar* XrgA = new scalar[N_EP*D];
    //scalar* XrgB = new scalar[N_EP*D];
    
    //use basis functions to get the physical locations of equivalence nodes
    for (int g = 0; g < N_EP; g++){
      for (int a = 0; a < D; a++){
	XrgA[g*D + a] = 0.0;
	XrgB[g*D + a] = 0.0;
	for (int k = 0; k < N_s; k++)
	  {
	    //Basis functions are multiplied by physical solution point coordinates
	    //to get physical equivalence node coordinates.
	    XrgA[g*D + a] += phi(FinalMapA[g], k) * xGeo_A[k*D + a];
	    XrgB[g*D + a] += phi(FinalMapB[g], k) * xGeo_B[k*D + a];
	  } //end k summation
      }} //end equivalence node and dimension loop

    if (verbose == 1)
      {
	if (D == 2)
	  {
	    printf("Physical coordinates, equivalence nodes\n");
	    printf("omA:\n");
	    for (int j = 0; j < N_EP; j++)
	      {
		printf("node %d: (x,y) = (%f,%f)\n", j, XrgA[j*D + 0], XrgA[j*D + 1]);
	      }
	      printf("omB:\n");
	      for (int j = 0; j < N_EP; j++)
		{
		printf("node %d: (x,y) = (%f,%f)\n", j, XrgB[j*D + 0], XrgB[j*D + 1]);
	      }
	      }
		if (D == 3)
		  {
		printf("Physical coordinates, equivalence nodes\n");
		printf("omA:\n");
		for (int j = 0; j < N_EP; j++)
		  {
		printf("node %d: (x,y,z) = (%f,%f,%f)\n", j, XrgA[j*D + 0], XrgA[j*D + 1], XrgA[j*D + 2]);
	      }
		  printf("omB:\n");
		for (int j = 0; j < N_EP; j++)
		  {
		printf("node %d: (x,y,z) = (%f,%f,%f)\n", j, XrgB[j*D + 0], XrgB[j*D + 1], XrgB[j*D + 2]);
	      }
	  }
      }

    scalar* RrgFace = new scalar[D*M_G];
    //Now get the recovery coordinates corresponding to xirgA and xirgB
    //TransformPhysical2Recovery(N_s, GQres, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, RrgA, RrgB, RrgFace);
    TransformPhysical2Recovery(N_s, N_EP, M_G, xGeo_A, xGeo_B, XrgA, XrgB, xFace, FaceNormal, RrgA, RrgB, RrgFace);

    for (int a = 0; a < D; a++)
      {
	if (FlagPeri[a] > 0)
	  {
	    if (verbose > 0) {printf("FlagPeri[%d] nonzero, removing periodicity fix on xGeo_B in Nodal Eq Cord\n",a);}
	    for (int k = 0; k < N_s; k++)
	      {
		xGeo_B[k*D + a] -= sign_shift*Ldomain;
	      }
	  }
      }

    if (verbose > 0) {printf("End of NodalEqCord after calling Recovery coordinates\n");}
    //printf("FlagPeri Size is %d\n", int(FlagPeri.size1()));
    if (verbose > 0) {printf("Check 1\n");}
    fflush(stdout);
    delete[] RrgFace;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.7\n");}
    delete[] FlagPeri;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.3\n");}
    delete[] FaceCent;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.4\n");}
    //delete[] XrgA;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.5\n");}
    //delete[] XrgB;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.6\n");}
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.1\n");}
    //delete[] sign_shift_relay;
    fflush(stdout);
    if (verbose > 0) {printf("Check 1.2\n");}
    //delete[] Ldomain_relay;
    fflush(stdout);

    //printf("CATASTROPHE! NodalEqNode not ready to move past 1D,2D\n");
  }
#endif
  if (verbose > 0) {printf("Check 3\n");}
}

void ICB_Cords_from_Physical(int N_s, int N_N, int GQresA, int GQresFace, int N_EP, scalar* xFace, scalar* xGeo_A, scalar* xGeo_B, scalar* xGQ_A, scalar* xGQ_B, scalar* XrgA, scalar* XrgB, /*begin output*/ scalar* rA_omA, scalar* rB_omB, scalar* rFace_omA, scalar* rFace_omB, scalar* r_rgB_omA, scalar* r_rgA_omB)
{
  /*!
    \brief For an interface, get the icb coordinate system of the two adjacent elements
    \param[in] N_s solution points per element
    \param[in] N_N sides per element
    \param[in] GQresA element quadratyre resoluytion
    \param[in] GQresFace interface quadrature resolution
    \param[in] N_EP non-dominant equivalence points per element
    \param[in] xFace physical coordinates of face quadrature points
    \param[in] xGeo_A the solution point nodes, elelement A
    \param[in] xGeo_B the solution point nodes, element B
    \param[in] xGQ_A quadrature point element A physical coordinates
    \param[in] xGQ_B physical coordinates of element B quadrature points
    \param[in] XrgA physical coordinates, element A nodal equivalence points
    \param[in] XrgB physical coordinates, element B nodal equivalence points
    \param[out] rA_omA A-native ICB coordinates for omA quadrature points
    \param[out] rB_omB B-native ICB coordinates for omB quadrature points
    \param[out] rFace_omA A-native ICB coordinates for face quadrature points
    \param[out] rFace_omB B-native ICB coordinates for face quadrature points
    \param[out] r_rgB_omA A-native ICB coordinates for omB equivalence points
    \param[out] r_rgA_omB B-native ICB coordinates for omA equivalence points
   */
  //Get the ICB coordinates for both elements
  int verbose = 0;
  //Step 1: Get the normal from the interface
  scalar normal[D];
  scalar Mrotate[D][D];
  scalar Mrotate_Inv[D][D];
#ifdef ONED
  normal[0] = 1.0; //don't care about +/- x direction here
#endif
#ifdef TWOD
  //https://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment
  normal[0] =  xFace[0*D + 1] - xFace[(GQresFace-1)*D + 1]; //y differential
  normal[1] = -xFace[0*D + 0] + xFace[(GQresFace-1)*D + 0]; //x differential
  scalar tangent[D];
  tangent[0] = xFace[0*D+0] - xFace[(GQresFace-1)*D + 0]; 
  tangent[1] = xFace[0*D+1] - xFace[(GQresFace-1)*D + 1]; 
  scalar magTan = 0.0;
  for (int a = 0; a < D; a++)
    {
      magTan += tangent[a]*tangent[a];
    }
  magTan = sqrt(magTan);
  for (int a = 0; a < D; a++)
    {
      tangent[a] = tangent[a] / magTan;
    }
#endif
#ifdef THREED
  //printf("WARNING: ICB_Cords_from_Physical not ready for 3D yet.\n");
  scalar* dir_r = new scalar[D];
  scalar* dir_s = new scalar[D];
  scalar* dir_t = new scalar[D];
  //Find the first four points in xGeo_A that lie in the plane.
  scalar* XV1 = new scalar[D];
  scalar* XV2 = new scalar[D];
  scalar* XV3 = new scalar[D];
  scalar* XV4 = new scalar[D];
  GrabPlaneVertices_3D(GQresFace, xGeo_A, xFace, XV1, XV2, XV3, XV4);
  //Now, get the normalized vectors defining the recovery coordinates.
  Reco3DVector_From_Vertices(XV1, XV2, XV3, XV4, dir_r, dir_s, dir_t);
  scalar tangent[D];
  scalar tan2gent[D];
  for (int a = 0; a < D; a++)
    {
      normal[a] = dir_r[a];
      tangent[a] = dir_s[a];
      tan2gent[a] = dir_t[a];
    }
  delete[] dir_r;
  delete[] dir_s;
  delete[] dir_t;
  delete[] XV1;
  delete[] XV2;
  delete[] XV3;
  delete[] XV4;
#endif
  scalar magNorm = 0.0;
  for (int a = 0; a < D; a++)
    {
      magNorm += normal[a]*normal[a];
    }
  magNorm = sqrt(magNorm);
  for (int a = 0; a < D; a++)
    {
      normal[a] = normal[a] / magNorm;
    }
  if (verbose > 0){
    printf("In ICB_Cords_from_Physical: normal = (");
    for (int a = 0; a < D; a++)
      {
	printf("%f, ",normal[a]);}
    printf(")\n");
#ifdef TWOD
    printf("\t\t\t\t\tAdditionally, tangent = (");
    for (int a = 0; a < D; a++){
      printf("%f, ",tangent[a]);}
    printf(")\n");
#endif
#ifdef THREED
    printf("\t\t\t\t\tAdditionally, tangent = (");
    for (int a = 0; a < D; a++){
      printf("%f, ",tangent[a]);}
    printf(")\n");
    printf("\t\t\t\t\tAnd the other tangent = (");
    for (int a = 0; a < D; a++){
      printf("%f, ",tan2gent[a]);}
    printf(")\n");
  }
#endif

  //Build the rotation matrix (for the last step of this subroutine).
  //The rotation matrix gets physical coordinates from ICB-directed coordinates;
  //I need to go opposite direction
#ifdef ONED
  Mrotate[0][0] = normal[0];
  Mrotate_Inv[0][0] = Mrotate[0][0];
#endif
#ifdef TWOD
  Mrotate[0][0] = normal[0];
  Mrotate[1][0] = normal[1];
  Mrotate[0][1] = tangent[0];
  Mrotate[1][1] = tangent[1];
  scalar detRot = Mrotate[0][0]*Mrotate[1][1] - Mrotate[0][1]*Mrotate[1][0];
  Mrotate_Inv[0][0] = 1.0/detRot * Mrotate[1][1];
  Mrotate_Inv[1][1] = 1.0/detRot * Mrotate[0][0];
  Mrotate_Inv[1][0] = -1.0/detRot * Mrotate[1][0];
  Mrotate_Inv[0][1] = -1.0/detRot * Mrotate[0][1];
#endif
#ifdef THREED
  Mrotate[0][0] = normal[0];
  Mrotate[1][0] = normal[1];
  Mrotate[2][0] = normal[2];
  Mrotate[0][1] = tangent[0];
  Mrotate[1][1] = tangent[1];
  Mrotate[2][1] = tangent[2];
  Mrotate[0][2] = tan2gent[0];
  Mrotate[1][2] = tan2gent[1];
  Mrotate[2][2] = tan2gent[2];
  //Use the matrix algebra library to deal with this inversion
  fullMatrix<scalar> Mrotate_Mat(3,3);
  for (int a1 = 0; a1 < D; a1++){
    for (int a2 = 0; a2 < D; a2++){
      Mrotate_Mat(a1,a2) = Mrotate[a1][a2]; }}
  //Invert it:
  fullMatrix<scalar> MRinv;
  Mrotate_Mat.invert(MRinv);
  for (int a1 = 0; a1 < D; a1++){
    for (int a2 = 0; a2 < D; a2++){
      Mrotate_Inv[a1][a2] = MRinv(a1,a2); }}
    
#endif

  //Step 2: Calculate centroids of both elements
  //Get centroids of the two elements
  scalar* Ac = new scalar[D]; //centroid of element A
  scalar* Bc = new scalar[D]; //centroid of element B
  //This isn't a strict centroid, I'm just taking arithmetic average of the nodes
  for (int a = 0; a < D; a++)
    {
      Ac[a] = 0.0;
      Bc[a] = 0.0;
    }
  for (int j = 0; j < N_s; j++)
    {
      for (int a = 0; a < D; a++)
	{
	  Ac[a] += xGeo_A[j*D + a];
	  Bc[a] += xGeo_B[j*D + a];
	}
    }
  for (int a = 0; a < D; a++)
    {
      Ac[a] = Ac[a] / (N_s + 0.0);
      Bc[a] = Bc[a] / (N_s + 0.0);
    }
  if (verbose > 0)
    {
      printf("Element Centroids in ICB_cords_from_Physical, pre-periodicity correction\n");
      for (int a = 0; a < D; a++)
	{
	  printf("Ac[%d]=%f\n",a,Ac[a]);
	}
       for (int a = 0; a < D; a++)
	{
	  printf("Bc[%d]=%f\n",a,Bc[a]);
	}
    }
  
  //Get the element scaling:
  scalar ICB_scaleA[D];
  scalar ICB_scaleB[D];
  for (int a = 0; a < D; a++){
    ICB_scaleA[a] = 0.0;
    ICB_scaleB[a] = 0.0; }
  int VertCount = N_N;
  if (D == 3)
    {
      VertCount = 8;
    }
  for (int k = 0; k < VertCount; k++) {
    for (int a = 0; a < D; a++) {
      ICB_scaleA[a] = fmax(ICB_scaleA[a], fabs(xGeo_A[k*D+a] - Ac[a]));
      ICB_scaleB[a] = fmax(ICB_scaleB[a], fabs(xGeo_B[k*D+a] - Bc[a])); }}
  if (verbose == 1)
    {
      printf("Element scaling:\n");
      for (int a = 0; a < D; a++)
	{
	  printf("ICB_scaleA[%d] = %f, ICB_scalB[%d] = %f\n", a, ICB_scaleA[a], a, ICB_scaleB[a]);
	}
    }

  //Deal with periodicity problem:
  int* FlagPeri = new int[D];
  int* sign_shift_relay = new int[1];
  scalar* Ldomain_relay = new scalar[1];
  PeriodicityFix(N_s, GQresFace, xGeo_A, xGeo_B, xFace, FlagPeri, sign_shift_relay, Ldomain_relay);
  int sign_shift = sign_shift_relay[0];
  scalar Ldomain = Ldomain_relay[0];
  if (verbose > 0){
    printf("In ICB_Cords_from_Physical. Exited Periodicity fix. sign_shift=%d, Ldomain=%f\n",sign_shift, Ldomain);
    for (int a = 0; a < D; a++){
      printf("FlagPeri[%d]=%d\n",a,FlagPeri[a]);}}

  //Get centroid-relevant coordinates
  //physical coordiate differences, element-centroid to element points:
  scalar* XrelA_fromA = new scalar[GQresA*D];
  scalar* XrelB_fromB = new scalar[GQresA*D];
  //physical coordiate differences, element-centroid to interface points
  scalar* XrelFace_fromA = new scalar[GQresFace*D];
  scalar* XrelFace_fromB = new scalar[GQresFace*D];

  //physical coordinate difference to non-dominant equivalence coordinates:
  scalar* Xrel_rgA_fromB = new scalar[N_EP*D];
  scalar* Xrel_rgB_fromA = new scalar[N_EP*D];
  //element relative coordinates:
  for (int g = 0; g < GQresA; g++){
    for (int a = 0; a < D; a++)      {
      XrelA_fromA[g*D+a] = xGQ_A[g*D + a] - Ac[a];
      XrelB_fromB[g*D+a] = xGQ_B[g*D + a] - Bc[a]; }}
    
  //face, relative coordinates:
  for (int g = 0; g < GQresFace; g++){
    for (int a = 0; a < D; a++){
      XrelFace_fromA[g*D+a] = xFace[g*D + a] - Ac[a];
      XrelFace_fromB[g*D+a] = xFace[g*D + a] - Bc[a]; }}
    
  //equivalence points, relative physical coordinates:
  for (int g = 0; g < N_EP; g++){
    for (int a = 0; a < D; a++){
      Xrel_rgB_fromA[g*D+a] = XrgB[g*D + a] - Ac[a];  
      Xrel_rgA_fromB[g*D+a] = XrgA[g*D + a] - Bc[a]; }}

  //Use information from periodicity fix to adjust element-B related coordinates
  for (int a = 0; a < D; a++)
    {
      if (FlagPeri[a] > 0)
	{
	  Bc[a] += sign_shift*Ldomain;
	  for (int g = 0; g < GQresFace; g++) {
	    XrelFace_fromB[g*D+a] -= sign_shift*Ldomain; }
	  for (int g = 0; g < N_EP; g++) {
	    //Xrel_rgB_fromA[g*D+a] -= sign_shift*Ldomain; 
	    /*Important issue: Upon import to this subroutine,
	      XrgB has already experienced periodicity correction,
	      so there is not need to repeat it here
	    */
	    Xrel_rgA_fromB[g*D+a] -= sign_shift*Ldomain; }
	}
    }

  if (verbose > 0)
    {
      printf("In ICB_Cords_from_Physical: After periodicity correction, before scaling:\n");
      printf("Element-Introspective:\n");
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelA_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelB_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-face:\n");
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelFace_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelFace_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-equivalence points\n");
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omA (points in B), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrel_rgB_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omB (points in A), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrel_rgA_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
    }
    

  //Scale the relative coordiates according to distance from centroid
  //to element vertices.
  for (int g = 0; g < GQresA; g++){
    for (int a = 0; a < D; a++)      {
      XrelA_fromA[g*D + a] /=  ICB_scaleA[a];
      XrelB_fromB[g*D + a] /=  ICB_scaleB[a];}}
  for (int g = 0; g < GQresFace; g++){
    for (int a = 0; a < D; a++){
      XrelFace_fromA[g*D + a] /= ICB_scaleA[a];
      XrelFace_fromB[g*D + a] /= ICB_scaleB[a]; }}
  for (int g = 0; g < N_EP; g++){
    for (int a = 0; a < D; a++){
      Xrel_rgB_fromA[g*D + a] /= ICB_scaleA[a];  
      Xrel_rgA_fromB[g*D + a] /= ICB_scaleB[a]; }}
      

  if (verbose > 0)
    {
      printf("After periodicity correction and scaling correction in ICB_Cords_from_Physical:\n");
      printf("Element-Introspective:\n");
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelA_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelB_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-face:\n");
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelFace_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrelFace_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-equivalence points\n");
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omA (points in B), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrel_rgB_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omB (points in A), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrel_rgA_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
    }

  //Transform scaled, physical-direction coordinates to ICB coordinate direction
  scalar* XrotA_fromA = new scalar[GQresA*D];
  scalar* XrotB_fromB = new scalar[GQresA*D];
  scalar* XrotFace_fromA = new scalar[GQresFace*D];
  scalar* XrotFace_fromB = new scalar[GQresFace*D];
  scalar* Xrot_rgA_fromB = new scalar[N_EP*D];
  scalar* Xrot_rgB_fromA = new scalar[N_EP*D];

  for (int g = 0; g < GQresA; g++) {
    for (int a = 0; a < D; a++) {
      XrotA_fromA[g*D+a] = 0.0;
      XrotB_fromB[g*D+a] = 0.0;
      for (int alpha = 0; alpha < D; alpha++) {
	//XrotA_fromA[g*D+a] += XrelA_fromA[g*D + alpha] * normal[a];
	//XrotB_fromB[g*D+a] += XrelB_fromB[g*D + alpha] * normal[a]; 
	XrotA_fromA[g*D+a] += XrelA_fromA[g*D + alpha] * Mrotate_Inv[a][alpha];
	XrotB_fromB[g*D+a] += XrelB_fromB[g*D + alpha] * Mrotate_Inv[a][alpha];}}}

  for (int g = 0; g < GQresFace; g++) {
    for (int a = 0; a < D; a++) {
      XrotFace_fromA[g*D+a] = 0.0;
      XrotFace_fromB[g*D+a] = 0.0;
      for (int alpha = 0; alpha < D; alpha++) {
	//XrotFace_fromA[g*D+a] += XrelFace_fromA[g*D + alpha] * normal[a];
	//XrotFace_fromB[g*D+a] += XrelFace_fromB[g*D + alpha] * normal[a]; 
	XrotFace_fromA[g*D+a] += XrelFace_fromA[g*D + alpha] * Mrotate_Inv[a][alpha];
	XrotFace_fromB[g*D+a] += XrelFace_fromB[g*D + alpha] * Mrotate_Inv[a][alpha]; }}}
  for (int g = 0; g < N_EP; g++) {
    for (int a = 0; a < D; a++) {
      Xrot_rgB_fromA[g*D+a] = 0.0;
      Xrot_rgA_fromB[g*D+a] = 0.0;
      for (int alpha = 0; alpha < D; alpha++) {
	//Xrot_rgB_fromA[g*D+a] += Xrel_rgB_fromA[g*D + alpha] * normal[a];
	//Xrot_rgA_fromB[g*D+a] += Xrel_rgA_fromB[g*D + alpha] * normal[a]; 
	Xrot_rgB_fromA[g*D+a] += Xrel_rgB_fromA[g*D + alpha] * Mrotate_Inv[a][alpha];
	Xrot_rgA_fromB[g*D+a] += Xrel_rgA_fromB[g*D + alpha] * Mrotate_Inv[a][alpha];  }}}

  if (verbose > 0) {
    printf("After periodicity correction, scaling correction, and rotation in ICB_Cords_from_Physical:\n");
      printf("Element-Introspective:\n");
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrotA_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresA; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrotB_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-face:\n");
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omA, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrotFace_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < GQresFace; g++)
	{
	  printf("omB, g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,XrotFace_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
      printf("Element-to-equivalence points\n");
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omA (points in B), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrot_rgB_fromA[g*D+a]);
	    }
	  printf(")\n");
	}
      for (int g = 0; g < N_EP; g++)
	{
	  printf("omB (points in A), g=%d:  (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("Xrel[a=%d]=%f, ",a,Xrot_rgA_fromB[g*D+a]);
	    }
	  printf(")\n");
	}
    }

  //Anything starting with Xrot is an ICB recovery coordinate.
  for (int j = 0; j < GQresA*D; j++) {
      rA_omA[j] = XrotA_fromA[j];
      rB_omB[j] = XrotB_fromB[j];}
  for (int j = 0; j < GQresFace*D; j++){
    rFace_omA[j] = XrotFace_fromA[j];
    rFace_omB[j] = XrotFace_fromB[j];}
  for (int j = 0; j < N_EP*D; j++)
    {
      r_rgB_omA[j] = Xrot_rgB_fromA[j];
      r_rgA_omB[j] = Xrot_rgA_fromB[j];}
  /*
  rA_omA = XrotA_fromA;
  rB_omB = XrotB_fromB;
  rFace_omA = XrotFace_fromA;
  rFace_omB = XrotFace_fromB;
  r_rgB_omA = Xrot_rgB_fromA;
  r_rgA_omB = Xrot_rgA_fromB;
  */
  //Subroutine concluded, whoever called it now has the proper ICB coordinates.
  delete[] Ac;
  delete[] Bc;
  delete[] FlagPeri;
  delete[] sign_shift_relay;
  delete[] Ldomain_relay;
  delete[] XrelA_fromA;
  delete[] XrelB_fromB;
  delete[] XrelFace_fromA;
  delete[] XrelFace_fromB;
  delete[] Xrel_rgA_fromB;
  delete[] Xrel_rgB_fromA;
  delete[] XrotA_fromA;
  delete[] XrotB_fromB;
  delete[] XrotFace_fromA;
  delete[] XrotFace_fromB;
  delete[] Xrot_rgA_fromB;
  delete[] Xrot_rgB_fromA;
  
}

void GetPSIxR_Biased_oneface_Nov2017(deck inputs, int N_N, TIMERS &timers, int N_s, int p, /*simpleMesh &m,*/ int N_basis, int M_G, const polynomialBasis *basis, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, scalar* PSIxR_biA, scalar* PSIxR_biB, scalar* PSIxR_biA_A, scalar* PSIxR_biA_B, scalar* PSIxR_biB_A, scalar* PSIxR_biB_B)
{
/*!
  \breif Build the Biased Recovery operators for a given interface, non-boundary. This routine uses interface coordinates to build recovery basis, as opposed to my newer approach, where each element's PSI basis is built to nearly replicate DG basis
  \param[in] N_s = DOF per element
  \param[in] p = polynomial order
  \param[in] N_basis = dimension of recovery basis
  \param[in] M_G = quadrature points per interface
  \param[in] xGeo_A = gmsh nodes for element A
  \param[in] xGeo_B = gmsh nodes for element B
  \param[in] xFace = quadrature node physical locations on the interface
  \param[in] FaceNormal = the normal vector OUT of the interface
  \param[out] PSIxR_biA: the interface's A-biased recovery operator
  \param[out] PSIxR_biB: the interface's B-biased recovery operator
  \param[out] PSIxR_biA_A: the half of A-baised PSIxR corresponding to element A
  \param[out] PSIxR_biA_B: the half of A-biased PSIxR corresponding to element B
  \param[out] PSIxR_biB_A: the half of B-baised PSIxR corresponding to element A
  \param[out] PSIxR_biB_B: the half of B-biased PSIxR corresponding to element B
*/

  //To get physical node locations: multiply nodal DG shape functions
  //by nodal coordinate values. First step is to populate
  //the DG basis on reference element
    int verbose = 0;
    if (verbose > 0) {printf("Entered  GetPSIxR_Biased_oneface_Nov2017. N_basis=%d\n",N_basis);}
    timers.start_timer(62);
    int Rescale_Mats = 0; //if 1, I adjust recovery system for better condition number.
    int CordStyle = 1; //0 for usual recovery coordinates, 1 for element-centered r/s 
    //1) Get reference basis on reference element at all N_superG quadrature points,
    //Quite a bit of this code is copied from early segments of main.cc because
    //I need to populate the DG basis functions
    // Get the method order
    int order = inputs.getOrder();
    bool order0 = false; if (order==0) {order0 = true; order = 1;}
    /*
    int elem_type;
    int face_type;
    int msh_hex;
    int msh_qua;
    int msh_tri;
    int msh_lin;
    int nsides; // this offsets j in buildInterfaces function
    int N_N;    // number of neighbors to an element
    scalar refArea; // area of reference element
    COPY_get_element_types(order, msh_hex, msh_qua, msh_tri, msh_lin);
    if     (inputs.getElemType() == "lin"){face_type = MSH_PNT, elem_type = msh_lin; nsides = 0; N_N = 2;}
    else if(inputs.getElemType() == "tri"){face_type = msh_lin, elem_type = msh_tri; nsides = 3; N_N = 3; refArea = 0.5;}
    else if(inputs.getElemType() == "qua"){face_type = msh_lin, elem_type = msh_qua; nsides = 4; N_N = 4; refArea = 4;}
    else if(inputs.getElemType() == "hex"){face_type = msh_qua, elem_type = msh_hex; nsides = 6; N_N = 6; refArea = 8;}
    else printf("Invalid element type in deck");
    //const polynomialBasis *basis  = polynomialBases::find(elem_type);  // for the element
    
    fullMatrix<double> RefPoints = basis->points;
    //  int arg0 = RefPoints.size0();
    int arg1 = RefPoints.size1(); //number of solution points (N_s)
    int arg2 = RefPoints.size2(); //spatial dimension (D)
    
    //  printf("test=%f\n",test);
    const std::vector<std::vector<int> > &closures = basis->closures;
    fullMatrix<double> points, weight;
     //Using very high resolution quadrature because I can
    if     (inputs.getElemType() == "lin") gaussIntegration::getLine(20, points, weight);
    else if(inputs.getElemType() == "tri") gaussIntegration::getTriangle(order*4+1, points, weight);
    else if(inputs.getElemType() == "qua") gaussIntegration::getQuad(order*4+1, points, weight);
    else if(inputs.getElemType() == "hex") gaussIntegration::getHexahedron(3*order + 1, points, weight);
    
    //Now, I have super-resolution gaussian quadratures and weights
    int N_superG = points.size1();
    if (verbose > 1)
    {
      printf("Quadrature nodes, ref element:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",points(g,a));
	    }
	  printf(")\n");
	}
    }

    //Get the DG basis over reference element
    fullMatrix<scalar> phiRef (N_superG , N_s); 
    fullMatrix<double> phiD (N_superG , N_s); 
    fullMatrix<scalar> dphi(N_superG*D,N_s);
    basis->f (points, phiD);
    for(int g = 0; g < N_superG; g++){
      for(int i = 0; i < N_s; i++){
	phiRef(g,i) = (scalar)phiD(g,i);
	//    printf("phiRef(%d,%d)=%f\n",g,i,phiRef(g,i));
      }
    }   

    double grads[N_s][3];  
    for(int g = 0; g < N_superG; g++){
      basis->df(points(g,0),points(g,1),points(g,2),grads);
      for(int alpha = 0; alpha < D; alpha ++){
	for(int i = 0; i < N_s; i++){
	  dphi(g*D+alpha,i) = (scalar)grads[i][alpha];  //see paper for indexing p.6
	}	  
      }    
    }
    */

    fullMatrix<double> RefPoints = basis->points;
    //  int arg0 = RefPoints.size0();
    int arg1 = RefPoints.size1(); //number of solution points (N_s)
    int arg2 = RefPoints.size2(); //spatial dimension (D)

    //Get the detJxW distribution for each element
  scalar* detJxW_A = new scalar[N_superG];
  scalar* detJxW_B = new scalar[N_superG];
  dg_detJ_OneElement(N_superG, xGeo_A, dphi, detJxW_A);
  dg_detJ_OneElement(N_superG, xGeo_B, dphi, detJxW_B);
  //That was just the detJ distribution; multiply by quadrature weights for det_JxW
  for (int g = 0; g < N_superG; g++){
    detJxW_A[g] = detJxW_A[g] * weight(g,0);
    detJxW_B[g] = detJxW_B[g] * weight(g,0); }
  //Also, this is an appropriate spot to scale the det_JxW values; they can
  //be quite small, and I don't want an overflow error on the LHS inversion
  
  scalar detNorm = 0.0;
  for (int g = 0; g < N_superG; g++)
    {
      detNorm += detJxW_A[g] + detJxW_B[g];}
  detNorm = detNorm / (2.0*N_superG);
  for (int g = 0; g < N_superG; g++){
    detJxW_A[g] = detJxW_A[g] / detNorm;
    detJxW_B[g] = detJxW_B[g] / detNorm; }


    //  printf("No Segfault yet A\n");
    //The solution basis has been populated on reference element: relay it to phiA, phiB storage
    fullMatrix<scalar> phi_A = phiRef;
    fullMatrix<scalar> phi_B = phiRef;
    scalar* xGQ_A = new scalar[N_superG*D];
    scalar* xGQ_B = new scalar[N_superG*D];
   
    //Get physical GQ locations in each element.
    //This procedure is particluar to a nodal basis; xGeo must be the solution node locations
    for (int g = 0; g < N_superG; g++)
      {
	for (int a = 0; a < D; a++)
	  {
	    xGQ_A[g*D + a] = 0.0;
	    xGQ_B[g*D + a] = 0.0;
	    for (int k = 0; k < N_s; k++)
	      {
		xGQ_A[g*D + a] += phi_A(g,k) * xGeo_A[k*D + a];
		xGQ_B[g*D + a] += phi_B(g,k) * xGeo_B[k*D + a];
	      }
	  }
      }
    
    //  printf("No Segfault yet C\n");
    if (verbose > 1)
      {
	printf("Qudarature node locations in GetPSIxR_Biased_oneface_Nov2017\n");
	printf("xGQ, element A:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",xGQ_A[g*D + a]);
	    }
	  printf(")\n");
	}
      printf("xGQ, element B:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",xGQ_B[g*D + a]);
	    }
	  printf(")\n");
	}
    }
    /*
    int N_EP; //nodal equivalence points per element
    if     (inputs.getElemType() == "lin") {N_EP=1;}
    else if(inputs.getElemType() == "tri") {N_EP=p+1;}
    else if(inputs.getElemType() == "qua") {N_EP=p+1;}
    else if(inputs.getElemType() == "hex") {N_EP=pow(p+1,2);}
    */

    timers.stop_timer(62);
    timers.start_timer(63);

    int N_icb = N_basis;//N_s + N_EP;
    int N_EP = N_basis - N_s; //nodal equivalence points in non-dominant element
    
    //Get the physical, recovery, and reference coordinates of the 
    //non-dominant equivalence points
    scalar* xirgA = new scalar[N_EP*D];
    scalar* xirgB = new scalar[N_EP*D];
    scalar* RrgA = new scalar[N_EP*D];
    scalar* RrgB = new scalar[N_EP*D];
    scalar* XrgA = new scalar[N_EP*D];
    scalar* XrgB = new scalar[N_EP*D];
      
    if (verbose > 0)
      {
	printf("xGeo_B before NodalEqCord call:\n");
	for (int k = 0; k < N_s; k++)
	  {
	    printf("node %d: (x,y,z) = (%f, %f, %f)\n",k,xGeo_B[k*D+0],xGeo_B[k*D+1],0.0);
	  }    
      }
    //Call NodalEqCord to get the physical and recovery coordinates of equivalence points.
    //I will not use the recovery coordinates, but instead use the ICB approach.
    if (verbose > 0){printf("Calling NodalEqCord from Get_PSIxR_biased_OneFace_Nov2017\n");}
    NodalEqCord(inputs, N_s, p, M_G, N_EP, N_N, basis, xGeo_A, xGeo_B, xFace, FaceNormal, xirgA, xirgB, RrgA, RrgB, XrgA, XrgB);
    if (verbose > 0) {printf("Back in PSIxR_biased_OneFace_Oct2017 after calling NodalEqCord\n");}
    if (verbose > 0)
      {
	printf("xGeo_B After NodalEqCord call:\n");
	for (int k = 0; k < N_s; k++)
	  {
	    printf("node %d: (x,y,z) = (%f, %f, %f)\n",k,xGeo_B[k*D+0],xGeo_B[k*D+1],0.0);
	  }    
	printf("Also in OneFace_Oct2017, here are the equivalence coordinates:\n");
	printf("Reference coordinates of equivalence nodes in omA:\n");
	for (int j = 0; j < N_EP; j++){
	  printf("node %d: Xi = (",j);
	  for (int a = 0; a < D; a++){
	    printf("%f, ",xirgA[j*D+a]);}
	  printf(")\n");}
	printf("Reference coordinates of equivalence nodes in omB:\n");
	for (int j = 0; j < N_EP; j++){
	  printf("node %d: Xi = (",j);
	  for (int a = 0; a < D; a++){
	    printf("%f, ",xirgB[j*D+a]);}
	  printf(")\n");}
	printf("Physical coordinates of equivalence nodes in omA:\n");
	for (int j = 0; j < N_EP; j++){
	  printf("node %d: Xi = (",j);
	  for (int a = 0; a < D; a++){
	    printf("%f, ",XrgA[j*D+a]);}
	  printf(")\n");}
	printf("Physical coordinates of equivalence nodes in omB:\n");
	for (int j = 0; j < N_EP; j++){
	  printf("node %d: Xi = (",j);
	  for (int a = 0; a < D; a++){
	    printf("%f, ",XrgB[j*D+a]);}
	  printf(")\n");}
      
	
      }
    timers.stop_timer(63);
    timers.start_timer(64);
    //Get element ICB coordinate system. r#_GQ_* = #-centered coordinate system populated in element *
    scalar* rA_GQ_A = new scalar[N_superG*D];
    scalar* rB_GQ_B = new scalar[N_superG*D];
    scalar* rA_Face = new scalar[M_G*D];
    scalar* rB_Face = new scalar[M_G*D];
    scalar* RB_rg_A = new scalar[N_EP*D];//B-dominant recovery coordinates for nodal equivalence in omA
    scalar* RA_rg_B = new scalar[N_EP*D];//A-dominant recovery coordinates for nodal equivalence in omB
    //these coordinates are centered at element centroid, and r coordinate is directed to closest face intersection
    if (CordStyle == 0)
      {
	//face-centered ICB coordinates, same as standard recovery operation
	TransformPhysical2Recovery(N_s, N_superG, M_G, xGeo_A, xGeo_B, xGQ_A, xGQ_B, xFace, FaceNormal, rA_GQ_A, rB_GQ_B, rA_Face);
	//Need to copy rA_Face in to rB_Face (both A and B use same coordinate system)
	for (int g = 0; g < M_G; g++){
	  for (int a = 0; a < D; a++){
	    rB_Face[g*D + a] = rA_Face[g*D+a];}}
	//What remains: the recovery coordinates of the non-dominant nodal equivalence nodes.
        //For this, copy in RrgA and RrgB
	for (int g = 0; g < N_EP; g++){
	  for (int a = 0; a < D; a++){
	    RB_rg_A[g*D+a] = RrgA[g*D+a];//recovery coordinates in element A
	    RA_rg_B[g*D+a] = RrgB[g*D+a];//recovery coordinates in element B 
	  }}
      }
    if (CordStyle == 1)
      {
	//element-centered ICB coordinates.
	if (verbose > 0){printf("Calling ICB_Cords_from_Physical from  Get_PSIxR_biased_OneFace_Oct2017\n");}
	ICB_Cords_from_Physical(N_s, N_N, N_superG, M_G, N_EP, xFace, xGeo_A, xGeo_B, xGQ_A, xGQ_B, XrgA, XrgB, rA_GQ_A, rB_GQ_B, rA_Face, rB_Face, RA_rg_B, RB_rg_A);
	if (verbose > 0){printf("Finished ICB_Cords_from_Physical call from  Get_PSIxR_biased_OneFace_Oct2017\n");}
      }
    if (verbose > 1)
    {
      printf("Got rA_GQ_A, rB_GQ_B, and rFace, check return:");
      printf("The ICB-recovery coordinates:\n");
      printf("rA, element A:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rA_GQ_A[g*D + a]);
	    }
	  printf(")\n");
	}
      printf("rB, element B:\n");
      for (int g = 0; g < N_superG; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rB_GQ_B[g*D + a]);
	    }
	  printf(")\n");
	}
    
      printf("rA, interface:\n");
      for (int g = 0; g < M_G; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rA_Face[g*D + a]);
	    }
	  printf(")\n");
	}

      printf("rB, interface:\n");
      for (int g = 0; g < M_G; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",rB_Face[g*D + a]);
	    }
	  printf(")\n");
	}

      printf("rA, equivalence nodes in omB:\n");
      for (int g = 0; g < N_EP; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",RA_rg_B[g*D + a]);
	    }
	  printf(")\n");
	}

      printf("rB, equivalence nodes in omA:\n");
      for (int g = 0; g < N_EP; g++)
	{
	  printf("node %d: (r1,r2,r3) = (",g);
	  for (int a = 0; a < D; a++)
	    {
	      printf("%f, ",RB_rg_A[g*D + a]);
	    }
	  printf(")\n");
	}
    }
    
     //Populate the ICB bases at: element quadrautre points, element equivalence points, face quadrature points
    scalar* PsiA_A = new scalar[N_icb*N_superG]; //A-dominant Psi in element A
    scalar* PsiB_B = new scalar[N_icb*N_superG]; //B-dominant Psi in element B
    scalar* PsiA_Face = new scalar[M_G*N_basis]; //A-dominant basis along interface quadrature points
    scalar* PsiB_Face = new scalar[M_G*N_basis]; //B-dominant basis along interface quadrature points
    scalar* PsiB_rgA = new scalar[N_basis*N_EP]; //B-dominant basis in element A
    scalar* PsiA_rgB = new scalar[N_basis*N_EP]; //A-dominant basis in element B

    PsiAtPoints(N_s, p, N_icb, N_N, N_superG, rA_GQ_A, PsiA_A);
    PsiAtPoints(N_s, p, N_icb, N_N, N_superG, rB_GQ_B, PsiB_B);
    
    PsiAtPoints(N_s, p, N_icb, N_N, M_G, rA_Face, PsiA_Face);
    PsiAtPoints(N_s, p, N_icb, N_N, M_G, rB_Face, PsiB_Face);

    PsiAtPoints(N_s, p, N_icb, N_N, N_EP, RA_rg_B, PsiA_rgB);
    PsiAtPoints(N_s, p, N_icb, N_N, N_EP, RB_rg_A, PsiB_rgA);

    if (verbose > 0)
      {
	printf("Basis Contents in Get_PSIxR_OneFace_Oct2017:\n");
	printf("PsiB_B:\n");
	for (int g = 0; g < N_superG; g++)
	  {
	    printf("g=%d, R=(",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("%f,",rB_GQ_B[g*D+a]);}
	    printf(")\n");
	    for (int k = 0; k < N_icb; k++)
	      {
		printf("\t\tPsiB_B[deg=%d] = %f\n",k, PsiB_B[k*N_superG + g]);}
	  }
	printf("PsiB_rgA:\n");
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d, R(in omA)=(",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("%f,",RB_rg_A[g*D+a]);}
	    printf(")\n");
	    for (int k = 0; k < N_icb; k++)
	      {
		printf("\t\tPsiB_rgA[deg=%d] = %f\n",k, PsiB_rgA[k*N_EP + g]);}
	  }
	printf("PsiB_Face:\n");
	for (int g = 0; g < M_G; g++)
	  {
	    printf("g=%d, rFace(from omB)=(",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("%f,",rB_Face[g*D+a]);}
	    printf(")\n");
	    for (int k = 0; k < N_icb; k++)
	      {
		printf("\t\tPsiB_Face[deg=%d] = %f\n",k, PsiB_Face[k*M_G + g]);}
	  }
      }
    timers.stop_timer(64);
    timers.start_timer(65);
    //Populate DG basis functions using xirg and Rrg for both elements
    //xi coordinates must be in matrif form for basis call
    fullMatrix<double> dummyNodesA (N_EP, D);
    fullMatrix<double> dummyNodesB (N_EP, D);
    for (int g = 0; g < N_EP; g++)
      {
	for (int a = 0; a < D; a++)
	  {
	    dummyNodesA(g,a) = xirgA[g*D + a];
	    dummyNodesB(g,a) = xirgB[g*D + a];
	  }
      }
    
    
    //DG basis at the nodal equivalence points:
    fullMatrix<scalar> phi_rgA (N_EP , N_s); 
    fullMatrix<double> phiD_rgA (N_EP , N_s); 
    
    basis->f (dummyNodesA, phiD_rgA);
    for(int g = 0; g < N_EP; g++){
      for(int i = 0; i < N_s; i++){
	phi_rgA(g,i) = (scalar)phiD_rgA(g,i);
	//    printf("phiRef(%d,%d)=%f\n",g,i,phiRef(g,i));
      }
    }   
    //Rinse and repeat for phi_rgB:
    fullMatrix<scalar> phi_rgB (N_EP , N_s); 
    fullMatrix<double> phiD_rgB (N_EP , N_s); 
    basis->f (dummyNodesB, phiD_rgB);
    for(int g = 0; g < N_EP; g++){
      for(int i = 0; i < N_s; i++){
	phi_rgB(g,i) = (scalar)phiD_rgB(g,i);
	//    printf("phiRef(%d,%d)=%f\n",g,i,phiRef(g,i));
      }
    }   
    timers.stop_timer(65);
    timers.start_timer(66);
    //Next, build the biased recovery system. Start with A-dominant system
    fullMatrix<scalar> LHS (N_basis, N_basis);
    fullMatrix<scalar> RHS (N_basis, 2*N_s);
    fullMatrix<scalar> PHIxdetJxW_A(N_s, N_superG);
    fullMatrix<scalar> PHIxdetJxW_B(N_s, N_superG);
    for (int k = 0; k < N_s; k++){
      for (int g = 0; g < N_superG; g++){
	PHIxdetJxW_A(k,g) = phi_A(g,k)*detJxW_A[g];
	PHIxdetJxW_B(k,g) = phi_B(g,k)*detJxW_B[g];
      }}
    //A-dominant case:
    //First N_s rows: equivalence over element A:
    for (int row = 0; row < N_s; row++)
      {
	//LHS: DG test basis vs icb solution basis:
	for (int col = 0; col < N_basis; col++)
	  {
	    LHS(row,col) = 0.0;
	    for (int g = 0; g < N_superG; g++)
	      {
		//phi call is quadrature node, then index
		//psi call is different, because it is.
		int slot = col*N_superG + g;
		//LHS(row,col) += phi_A(g, row) * PsiA_A[slot] * weight(g,0);
		//LHS(row,col) += phi_A(g, row) * PsiA_A[slot] * detJxW_A[g];
		LHS(row,col) += PsiA_A[slot] * PHIxdetJxW_A(row,g);
	      }
	  }
	//RHS: DG test basis vs DG solution basis:
	for (int col = 0; col < N_s; col++)
	  {
	    RHS(row,col) = 0.0;
	    for (int g = 0; g < N_superG; g++)
	      {
		//phi call is quadrature node, then index,
		//RHS(row,col) += phi_A(g, row) * phi_A(g, col) * weight(g,0);
		//RHS(row,col) += phi_A(g, row) * phi_A(g, col) * detJxW_A[g];
		RHS(row,col) += phi_A(g, col) * PHIxdetJxW_A(row,g);
	      }
	  }
      }
    //Next (N_basis - Nsolut) rows: nodal equivalence in element B
    for (int row = N_s; row < N_basis; row++)
      {
	int gIndex = row - N_s; //testing node in non-dominant element
	//Left side: the nodal value of ICB solution
	for (int col = 0; col < N_basis; col++) //ICB solution mode
	  {
	    int slot = col*N_EP + gIndex;
	    LHS(row,col) = PsiA_rgB[slot];
	  }
	//right side: the nodal value of DG solution
	for (int col = 0; col < N_s; col++)
	  {
	    RHS(row, col + N_s) = phi_rgB(gIndex, col);
	  }
      }
    
    //PEj 10/22/2017: Ading an option to adjust matrix scaling;
    //the nodal equivalence rows generally carry much larger values
    //than the galerkin equivalence rows
    if (Rescale_Mats == 1)
      {
	scalar sum_galerkin = 0.0;
	scalar sum_nodal = 0.0;
	for (int row = 0; row < N_s; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		sum_galerkin += fabs(LHS(row,col));
	      }
	  }
	for (int row = N_s; row < N_icb; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		sum_nodal += fabs(LHS(row,col));
	      }
	  }
	scalar ratio = sum_nodal/sum_galerkin;
	//Now, multiply the Galerkin equivalence rows by ratio
	for (int row = 0; row < N_s; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		LHS(row,col) = LHS(row,col) * ratio;
	      }
	    for (int col = 0; col < 2*N_s; col++)
	      {
		RHS(row,col) = RHS(row,col) * ratio;
	      }
	  }
      }

    //Invert the system and multiply by basis functions to get the Biased recovery operator
    
    //found this inversion command in dg/dg_functions.cc/dg_inverse_mass_matrix.
    //Brought in french as well to match style of code
    // Inverser la matrice de masse.
    if (verbose > 0) {printf("About to invert LHS, A-biased recovery\n");
      printf("The LHS:\n");
      for (int row = 0; row < N_basis; row++)
	{
	  printf("row %d:  (",row);
	  for (int col = 0; col < N_basis; col++)
	    {
	      printf("%6.4f, ",LHS(row,col));
	    }
	  printf("\n");
	}
    }
    fullMatrix<scalar> LHS_inv(N_basis,N_basis);
    LHS.invert(LHS_inv);
    if (verbose > 0) {printf("Inverted LHS, A-biased recovery\n");}
    
    fullMatrix<scalar> Hope_Iden(N_basis, N_basis);
    if (verbose > 0){
    Hope_Iden.gemm(LHS,LHS_inv);
    }
    //Recovery matrix is LHS^-1 * RHS
    //found the matrix multiplication command in main.cc, subroutine  LagMono2DTransformsCartesian, see last line
    fullMatrix<scalar> MatReco(N_basis, 2*N_s);
    MatReco.gemm(LHS_inv, RHS);
    //Now, MatReco is the interface's recovery matrix
    //Now, multiply the recovery matrix by the recovery basis (at interface) for the PSIxR structure.
    //How to use PSIxR: f(r=0) = PSIxR times (DOF_A; DOF_B)
    
    fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
    for (int g = 0; g < M_G; g++)
      {
	for (int col = 0; col < 2*N_s; col++)
	  {
	    scalar sum = 0;
	    for (int jj = 0; jj < N_basis; jj++)
	      {
		sum += PsiA_Face[jj*M_G + g] * MatReco(jj,col);
	      }
	    PSIxR_matrix(g,col) = sum;
	  }
      }

    //Store PSIxR_matrix as PSIxR_biA:
    //Serialize the PSIxR matrix
    //The organization here must pair properly with use of PSIxR in src/dg/kernels_phil.cu
    for (int g = 0; g < M_G; g++)
      {
	for (int col = 0; col < 2*N_s; col++)
	  {
	    int slot = g*2*N_s + col;
	    PSIxR_biA[slot] = PSIxR_matrix(g,col);
	  }
      }
  //Serialize the split matrix components
    for (int g = 0; g < M_G; g++)
      {
	for (int col = 0; col < N_s; col++)
	  {
	    int slot = g*N_s + col;
	    PSIxR_biA_A[slot] = PSIxR_matrix(g,col + 0);
	    PSIxR_biA_B[slot] = PSIxR_matrix(g,col + N_s);
	  }
      }
    if (verbose > 0)
      {
	printf("Basis contents (created by subroutine find(int tag) in polynomialBasis.cc:\n");
	printf("Size of RefPoints matrix: arg1=%d, arg2=%d\n", arg1, arg2);
	printf("The Basis reference coordinates:\n");
	for (int k = 0; k < N_s; k++)
	  {
	    printf("solution node %d: ",k);
	    for (int a = 0; a < D; a++)
	      {
		printf("xi_%d = %f,\t\t",a,RefPoints(k,a));
	      }
	    printf("\n");
	  }
	printf("\n");
	printf("In GetPSIxR_Biased_oneface. number of non-dominant equivalence points N_EP=%d\n", N_EP);
	printf("In GetPSIxR_Biased_oneface. The Nodal Equivalence coordinates:\n");
	printf("El A:\n");
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Reco(g=%d,a=%d)=%f, ",g,a,RA_rg_B[g*D+a]);
	      }
	    printf("\n");
	  }
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Xi(g=%d,a=%d)=%f, ",g,a,xirgA[g*D+a]);
	      }
	    printf("\n");
	  }
	printf("El B:\n");
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Reco(g=%d,a=%d)=%f, ",g,a,RB_rg_A[g*D+a]);
	      }
	    printf("\n");
	  }
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Xi(g=%d,a=%d)=%f, ",g,a,xirgB[g*D+a]);
	      }
	    printf("\n");
	  }
	if (verbose > 1)
	  {
	    printf("PsiA_A from ICB_PsiAtPoints:\n");
  
	    printf("---Element A:---\n");
	    for (int k = 0; k < N_basis; k++)
	      {
		printf("index %d:\n",k);
		for (int g = 0; g < N_superG; g++)
		  {
#ifdef ONED
		    printf("\t\tPsi[%d][g=%d](r1=%f) = %f\n", k, g, rA_GQ_A[g*D+0],PsiA_A[k*N_superG + g]);
#endif
#ifdef TWOD
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rA_GQ_A[g*D+0], rA_GQ_A[g*D+1],PsiA_A[k*N_superG + g]);
#endif
#ifdef THREED
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f,r3=%f) = %f\n", k, g, rA_GQ_A[g*D+0], rA_GQ_A[g*D+1], rA_GQ_A[g*D+2], PsiA_A[k*N_superG + g]);
#endif
		  }
	      }
  
	    printf("\n---The Interface:---\n");
	    for (int k = 0; k < N_basis; k++)
	      {
		printf("index %d:\n",k);
		for (int g = 0; g < M_G; g++)
		  {
#ifdef ONED
		    printf("\t\tPsi[%d][g=%d](r1=%f) = %f\n", k, g, rA_Face[g*D+0] ,PsiA_Face[k*M_G + g]);
#endif
#ifdef TWOD
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rA_Face[g*D+0], rA_Face[g*D+1],PsiA_Face[k*M_G + g]);
#endif
#ifdef THREED
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f,r3=%f) = %f\n", k, g, rA_Face[g*D+0], rA_Face[g*D+1], rA_Face[g*D+2], PsiA_Face[k*M_G + g]);
#endif
		  }
	      }
	  }

	printf("Results of A-dominant setup:\n");
	printf("The LHS\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",LHS(row,col));
	      }
	    printf("\n");
	  }
	printf("The RHS\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%6.4f, ",RHS(row,col));
	      }
	    printf("\n");
	  }
	printf("The LHS_inverse\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",LHS_inv(row,col));
	      }
	    printf("\n");
	  }
	printf("LHS x LHS_inverse\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",Hope_Iden(row,col));
	      }
	    printf("\n");
	  }
	printf("The Recovery matrix:\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d: ",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%f, ",MatReco(row,col));
	      }
	    printf("\n");
	  }
	printf("The PSixR matrix:\n");
	for (int row = 0; row < M_G; row++)
	  {
	    printf("row [g=%d]: ",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%f, ",PSIxR_matrix(row,col));
	      }
	    printf("\n");
	  }
      }


    //====================================================================
    //===================================================================
    //===================================================================
    //Redo the operation for B-biased recovery operation
    //Zero a few things
    
    for (int row = 0; row < N_basis; row++){
      for (int col = 0; col < N_basis; col++) {
	LHS(row,col) = 0.0; 
	LHS_inv(row,col) = 0.0;
	Hope_Iden(row,col) = 0.0;}
      for (int col = 0; col < 2*N_s; col++) {
	RHS(row,col) = 0.0;
	MatReco(row,col) = 0.0; }}
    
    for (int g = 0; g < M_G; g++) {
      for (int col = 0; col < 2*N_s; col++) {
	PSIxR_matrix(g,col) = 0.0; }}

    //B-dominant case:
    //first (N_basis - Nsolut) rows: nodal equivalence in element A
    for (int row = 0; row < N_basis - N_s; row++)
      {
	int gIndex = row; //testing node in non-dominant element
	//Left side: the nodal value of B-dominant ICB solution
	for (int col = 0; col < N_basis; col++)
	  {
	    int slot = col*N_EP + gIndex;
	    LHS(row,col) = PsiB_rgA[slot];
	  }
	//right side: the nodal value of DG solution, element A
	for (int col = 0; col < N_s; col++)
	  {
	    RHS(row, col + 0) = phi_rgA(gIndex, col);
	  }
      }
    
    //Last N_s rows: equivalence over element B:
    for (int row = N_basis - N_s; row < N_basis; row++)
      {
	int testIndex = row - (N_basis - N_s);
	//LHS: DG test basis vs icb solution basis:
	for (int col = 0; col < N_basis; col++)
	  {
	    LHS(row,col) = 0.0;
	    for (int g = 0; g < N_superG; g++)
	      {
		//phi call is quadrature node, then index
		//psi call is different, because it is.
		int slot = col*N_superG + g;
		//LHS(row,col) += phi_B(g, testIndex) * PsiB_B[slot] * weight(g,0);
		//LHS(row,col) += phi_B(g, testIndex) * PsiB_B[slot] * detJxW_B[g];
		LHS(row,col) += PsiB_B[slot] * PHIxdetJxW_B(testIndex,g);
	    }
	  }
	//RHS: DG test basis vs DG solution basis:
	for (int col = N_s; col < 2*N_s; col++)
	  {
	    int solIndex = col - N_s;
	    RHS(row,col) = 0.0;
	    for (int g = 0; g < N_superG; g++)
	      {
		//phi call is quadrature node, then index,
		//RHS(row,col) += phi_B(g, testIndex) * phi_B(g, solIndex) * weight(g,0);
		//RHS(row,col) += phi_B(g, testIndex) * phi_B(g, solIndex) * detJxW_B[g];
		RHS(row,col) += phi_B(g, solIndex) * PHIxdetJxW_B(testIndex,g);
	      }
	  }
      }
    
     //PEj 10/22/2017: Ading an option to adjust matrix scaling;
    //the nodal equivalence rows generally carry much larger values
    //than the galerkin equivalence rows
    
    if (Rescale_Mats == 1)
      {
	scalar sum_galerkin = 0.0;
	scalar sum_nodal = 0.0;
	for (int row = N_EP; row < N_icb; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		sum_galerkin += fabs(LHS(row,col));
	      }
	  }
	for (int row = 0; row < N_EP; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		sum_nodal += fabs(LHS(row,col));
	      }
	  }
	scalar ratio = sum_nodal/sum_galerkin;
	//Now, multiply the Galerkin equivalence rows by ratio
	for (int row = N_EP; row < N_icb; row++)
	  {
	    for (int col = 0; col < N_icb; col++)
	      {
		LHS(row,col) = LHS(row,col) * ratio;
	      }
	    for (int col = 0; col < 2*N_s; col++)
	      {
		RHS(row,col) = RHS(row,col) * ratio;
	      }
	  }
      }
    

    //Invert the system and multiply by basis functions to get the Biased recovery operator
    
    //found this inversion command in dg/dg_functions.cc/dg_inverse_mass_matrix.
    //Brought in french as well to match style of code
    // Inverser la matrice de masse.
    //fullMatrix<scalar> LHS_inv(N_basis,N_basis);
    if (verbose > 0) {printf("About to invert LHS, B-biased recovery\n");}
    LHS.invert(LHS_inv);
    if (verbose > 0) {printf("Inverted LHS, B-biased recovery\n");}
    //fullMatrix<scalar> Hope_Iden(N_basis, N_basis);
    if (verbose > 0){
    Hope_Iden.gemm(LHS,LHS_inv);
    }
    //Recovery matrix is LHS^-1 * RHS
    //found the matrix multiplication command in main.cc, subroutine  LagMono2DTransformsCartesian, see last line
    //fullMatrix<scalar> MatReco(N_basis, 2*N_s);
    MatReco.gemm(LHS_inv, RHS);
    //Now, MatReco is the interface's recovery matrix
    //Now, multiply the recovery matrix by the recovery basis (at interface) for the PSIxR structure.
    //How to use PSIxR: f(r=0) = PSIxR times (DOF_A; DOF_B)
    
    //fullMatrix<scalar> PSIxR_matrix (M_G, 2*N_s);
    for (int g = 0; g < M_G; g++)
      {
	for (int col = 0; col < 2*N_s; col++)
	  {
	    scalar sum = 0;
	    for (int jj = 0; jj < N_basis; jj++)
	      {
	      sum += PsiB_Face[jj*M_G + g] * MatReco(jj,col);
	      }
	    PSIxR_matrix(g,col) = sum;
	  }
      }

  //Store PSIxR_matrix as PSIxR_biB:
  //Serialize the PSIxR matrix
  //The organization here must pair properly with use of PSIxR in src/dg/kernels_phil.cu
  for (int g = 0; g < M_G; g++)
    {
      for (int col = 0; col < 2*N_s; col++)
	{
	  int slot = g*2*N_s + col;
	  PSIxR_biB[slot] = PSIxR_matrix(g,col);
	}
    }
  //Serialize the split matrix components
  for (int g = 0; g < M_G; g++)
    {
      for (int col = 0; col < N_s; col++)
	{
	  int slot = g*N_s + col;
	  PSIxR_biB_A[slot] = PSIxR_matrix(g,col + 0);
	  PSIxR_biB_B[slot] = PSIxR_matrix(g,col + N_s);
	}
    }

  if (verbose > 0)
      {
	printf("Basis contents (created by subroutine find(int tag) in polynomialBasis.cc:\n");
	printf("Size of RefPoints matrix: arg1=%d, arg2=%d\n", arg1, arg2);
	printf("The Basis reference coordinates:\n");
	for (int k = 0; k < N_s; k++)
	  {
	    printf("solution node %d: ",k);
	    for (int a = 0; a < D; a++)
	      {
		printf("xi_%d = %f,\t\t",a,RefPoints(k,a));
	      }
	    printf("\n");
	  }
	printf("\n");
	printf("In GetPSIxR_Biased_oneface. number of non-dominant equivalence points N_EP=%d\n", N_EP);
	printf("In GetPSIxR_Biased_oneface. The Nodal Equivalence coordinates:\n");
	printf("El A:\n");
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Reco(g=%d,a=%d)=%f, ",g,a,RA_rg_B[g*D+a]);
	      }
	    printf("\n");
	  }
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Xi(g=%d,a=%d)=%f, ",g,a,xirgA[g*D+a]);
	      }
	    printf("\n");
	  }
	printf("El B:\n");
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Reco(g=%d,a=%d)=%f, ",g,a,RB_rg_A[g*D+a]);
	      }
	    printf("\n");
	  }
	for (int g = 0; g < N_EP; g++)
	  {
	    printf("g=%d: ",g);
	    for (int a = 0; a < D; a++)
	      {
		printf("Xi(g=%d,a=%d)=%f, ",g,a,xirgB[g*D+a]);
	      }
	    printf("\n");
	  }
	if (verbose > 1)
	  {
	    printf("PsiB_B from ICB_PsiAtPoints:\n");
  
	    printf("---Element B:---\n");
	    for (int k = 0; k < N_basis; k++)
	      {
		printf("index %d:\n",k);
		for (int g = 0; g < N_superG; g++)
		  {
#ifdef ONED
		    printf("\t\tPsi[%d][g=%d](r1=%f) = %f\n", k, g, rB_GQ_B[g*D+0],PsiB_B[k*N_superG + g]);
#endif
#ifdef TWOD
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rB_GQ_B[g*D+0], rB_GQ_B[g*D+1],PsiB_B[k*N_superG + g]);
#endif
#ifdef THREED
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f,r3=%f) = %f\n", k, g, rB_GQ_B[g*D+0], rB_GQ_B[g*D+1], rB_GQ_B[g*D+2], PsiB_B[k*N_superG + g]);
#endif
		  }
	      }
  
	    printf("\n---The Interface:---\n");
	    for (int k = 0; k < N_basis; k++)
	      {
		printf("index %d:\n",k);
		for (int g = 0; g < M_G; g++)
		  {
#ifdef ONED
		    printf("\t\tPsi[%d][g=%d](r1=%f) = %f\n", k, g, rB_Face[g*D+0] ,PsiB_Face[k*M_G + g]);
#endif
#ifdef TWOD
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f) = %f\n", k, g, rB_Face[g*D+0], rB_Face[g*D+1],PsiB_Face[k*M_G + g]);
#endif
#ifdef THREED
		    printf("\t\tPsi[%d][g=%d](r1=%f,r2=%f,r3=%f) = %f\n", k, g, rB_Face[g*D+0], rB_Face[g*D+1], rB_Face[g*D+2], PsiB_Face[k*M_G + g]);
#endif
		  }
	      }
	  }

	printf("Results of B-dominant setup:\n");
	printf("The LHS\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",LHS(row,col));
	      }
	    printf("\n");
	  }
	printf("The RHS\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%6.4f, ",RHS(row,col));
	      }
	    printf("\n");
	  }
	printf("The LHS_inverse\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",LHS_inv(row,col));
	      }
	    printf("\n");
	  }
	printf("LHS x LHS_inverse\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d:  (",row);
	    for (int col = 0; col < N_basis; col++)
	      {
		printf("%6.4f, ",Hope_Iden(row,col));
	      }
	    printf("\n");
	  }
	printf("The Recovery matrix:\n");
	for (int row = 0; row < N_basis; row++)
	  {
	    printf("row %d: ",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%f, ",MatReco(row,col));
	      }
	    printf("\n");
	  }
	printf("The PSixR matrix:\n");
	for (int row = 0; row < M_G; row++)
	  {
	    printf("row [g=%d]: ",row);
	    for (int col = 0; col < 2*N_s; col++)
	      {
		printf("%f, ",PSIxR_matrix(row,col));
	      }
	    printf("\n");
	  }
      }
  timers.stop_timer(66);
  if (verbose > 0) {printf("GetPSIxR_Biased_oneface_NOv2017: Preparing for pointer deletion\n");}
  if (verbose > 0) {fflush(stdout); printf("check A\n");}
  delete[] xGQ_A;
  if (verbose > 0) {fflush(stdout); printf("check B\n");}
  delete[] xGQ_B;
  if (verbose > 0) {fflush(stdout); printf("check C\n");}
  delete[] xirgA;
  if (verbose > 0) {fflush(stdout); printf("check J\n");}
  delete[] xirgB;
  if (verbose > 0) {fflush(stdout); printf("check K\n");}
  delete[] RrgA; delete[] RrgB;
  delete[] XrgA; delete[] XrgB;
  delete[] rA_GQ_A; delete[] rB_GQ_B;
  delete[] rA_Face; delete[] rB_Face;
  delete[] RB_rg_A; delete[] RA_rg_B;
  delete[] PsiA_A; delete[] PsiB_B;
  delete[] PsiA_Face; delete[] PsiB_Face;
  delete[] PsiB_rgA; delete[] PsiA_rgB;
  delete[] detJxW_A; delete[] detJxW_B;

  if (verbose > 0) {printf("GetPSIxR_Biased_oneface_Nov2017: Finished with pointer deletion\n");}

} //End of the subroutine
