/*!
  \file kernels_phil.cu
  \brief Kernels to be used for DG operations, specifically for mixed-form approach
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "kernels.h"
#include "kernels_phil.h"
//#include<cuda_runtime.h>
//#include "cuPrintf.cu"
#include<stdio.h> //need this for printf, assuming sufficient nvcc compiler
//==========================================================================

#ifdef USE_GPU
__global__ void print_kernel() {
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}
#endif


#ifdef USE_GPU
__device__ double ATOMICADD(double* address, double val)
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


//Some stuff that only gets used when I run on GPU (Phil doesn't like mixing cpu/gpu subroutines)
#ifdef USE_GPU
__global__ void Uhat_to_UhCommonHalf_GPU(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf)
{
    /*!
    \brief Send each element's UhCommon contribution to the interface. Another routine must add the contributions.
    \This is the pure GPU version because I want to get better at compiling
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

  */
  //Big picture: See Uhat_to_UhCommonHalf. Multiplying PSIxR_elemwise by Uhat to get interface contributions
  scalar sum = 0;

  //Configuration 1:
  //int e = blockIdx.x*blockDim.x + threadIdx.z;
  int e = blockIdx.x*blkE + threadIdx.z;
  int s = threadIdx.y;
  //End configuratin 1.

  /*
  //Configuration 2:
  int e = blockIdx.z*blockDim.z + threadIdx.z;
  int s = blockIdx.y*blockDim.y + threadIdx.y;
  //End configuration 2.
  */

  int tGlo = RecoPair[e*(D+N_N) + s]; 
  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
  //  printf("e=%d,s=%d,tGlo=%d,ArgSide=%d\n", e, s, tGlo, ArgSide);
  //Configuration 1:
  int g = threadIdx.x;
  //end configuration 1.

  /*
  //Configuration 2:
  int g = blockIdx.x*blockDim.x + threadIdx.x;
  //End configuration 2.
  */
  /*
    The problem: I need a += approach for setting UhCommonHalf because MANY elements send contributions
    to side 0 of interface 0. I could not determine another approach to dealing with the duplicate
    periodic interfaces. += is tricky for GPU. Must use atomicAdd from kernels.cu

    Don't know why, but calling atomicAdd results in compile error. So, I copied it
    in to this .cu as ATOMICADD
  */

  if (e >= N_E || s >= (D+N_N) || g >= M_G) return;
  for(int fc = 0; fc < N_F; fc++)
    {
      for (int k = 0; k < N_s; k++)
	{
	  sum += PSIxR_elemwise[e*((D+N_N)*M_G*N_s) + s*(M_G*N_s) + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	}
      //have to use a += because interface 0 gets many zeros sent to it. 
      /*
      sum += UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide];
      UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide] =  sum; 
      */
      ATOMICADD(&UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide], sum);
      //  printf("\tprintf(e=%d,s=%d,argside=%d); sum to UhCommonHalf(%d,%d,%d,side=%d) = %f; Actual UhCommonHalf = %f\n",e,s,ArgSide,tGlo,g,fc,ArgSide,sum,UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide]);
    } 
}

__global__ void UhCommonHalf_to_UhCommon_GPU(int M_T, int M_G, scalar* UhCommonHalf, scalar* UhCommon)
{
  /*!
    \brief combine the UhCommonHalf contributions from each side of each interface to get UhCommon on the face
    
    \param[in] M_T number of interfaces
    \param[out] UhCommon: interface solution, specifically at quadrature points
   */
  int t = blockIdx.x*blkT + threadIdx.z;
  int g = threadIdx.y;
  int fc = threadIdx.x;
  if (t >= M_T || g >= M_G || fc >= N_F) return;
  //  printf("Left-hand value = %f, right-hand = %f, ",UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
  UhCommon[t*M_G*N_F + g*N_F + fc] = UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0] + UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1];
  // printf("Left-hand value = %f, right-hand = %f, UhCommon(%d,%d,%d) = %f\n", UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1], t,g,fc,UhCommon[t*M_G*N_F + g*N_F + fc]);
}

__global__ void Uhat_to_GradCommonHalf_GPU(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf)
{
  scalar sum = 0;
  int e = blockIdx.x*blkE + threadIdx.z;
  int s = threadIdx.y;
  int tGlo = RecoPair[e*(D+N_N) + s]; 
  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
  //  printf("e=%d,s=%d,tGlo=%d,ArgSide=%d\n", e, s, tGlo, ArgSide);
  int g = threadIdx.x;
  
  if (e >= N_E || s >= (D+N_N) || g >= M_G) return;
  for(int fc = 0; fc < N_F; fc++)
    {
      for (int a = 0; a < D; a++)
	{
     
	  //printf("fc=%d,a=%d\n",fc,a);
	  for (int k = 0; k < N_s; k++)
	    {
	      //	      gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	      //	      gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	      sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	      //sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	    }
	  //gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide] = sum;
	  ATOMICADD(&gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide], sum);
	  sum = 0;
	  //printf("sum=%f\n",sum);
	  //gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide]  =  sum; 
	  //	  printf("\tprintf(e=%d,s=%d,g=%d,a=%d,argside=%d); sum to gradCommonHalf(%d,%d,%d,a=%d,side=%d) = %f\n",e,s,g,a,ArgSide,tGlo,g,fc,a,ArgSide,gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide]);
	}
    }
}

__global__ void GradCommonHalf_to_GradCommon_GPU(int M_T, int M_G, scalar* GradCommonHalf, scalar* GradCommon)
{
  int t = blockIdx.x*blkT + threadIdx.z;
  int g = threadIdx.y;
  int fc = threadIdx.x;
  if (t >= M_T || g >= M_G || fc >= N_F) return;
  for (int a = 0; a < D; a++)
    {
      GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] = GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 0] + GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 1];
      //    printf("GradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n", t,g,fc,a,GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
    }
}

__global__ void Uhat_to_UicbHalf_GPU(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_biased_elemwise, scalar* Uhat, scalar* UicbHalf)
{
  printf("ENTERING Uhat_to_UicbHalf_GPU\n"); //fflush(stdout);
   //Big picture: See Uhat_to_UicbHalf. Multiplying PSIxR_biased_elemwise by Uhat to get interface contributions for biased ICBN solutions 
  scalar sum = 0;

  //Configuration 1:

  int e = blockIdx.x*blkE + threadIdx.z;
  if (e < N_E)
    {
      int s = threadIdx.y; //goes to D+N_N
      int dom = threadIdx.x; //goes to 2, for number of elements per interface (the element passes information for both A-biased and B-biased reconstructions)
      
      int tGlo = RecoPair[e*(D+N_N) + s]; //finds theglobal interface for local (e,s) pair
      int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
      
      
      //int g = threadIdx.x; //goes to M_G (not N_G, because I'm working on an interface of element)
      
      /*
	The problem: I need a += approach for setting UicbHalf because MANY elements send contributions
	to side 0 of interface 0. I could not determine another approach to dealing with the duplicate
	periodic interfaces. += is tricky for GPU. Must use atomicAdd from kernels.cu

	Don't know why, but calling atomicAdd results in compile error. So, I copied it
	in to this .cu as ATOMICADD
      */

      // if (e >= N_E || s >= (D+N_N) || dom >= 2) return;
      for (int g = 0; g < M_G; g++)
	{
	  for(int fc = 0; fc < N_F; fc++)
	    {
	      sum = 0;
	      for (int k = 0; k < N_s; k++)
		{
		  sum += 
		    PSIxR_biased_elemwise[e*((D+N_N)*2*M_G*N_s) + s*(2*M_G*N_s) + dom*M_G*N_s + g*N_s + k] * 
		    Uhat[(e*N_F+fc)*N_s+k];
		}
	      printf("e=%d,s=%d, dom=%d, g=%d: Uhat summation to UicbHalf(t=%d, g=%d, dom=%d, fc=%d, side=%d) = % f\n",e,s,dom,g,tGlo,g,dom,fc,ArgSide,sum);
	      ATOMICADD(&UicbHalf[tGlo*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + ArgSide], sum);
	    }
	  //have to use a += because interface 0 gets many zeros sent to it. 
	} 
    }
}

__global__ void UicbHalf_to_Uicb_GPU(int M_T, int M_G, scalar* UicbHalf, scalar* Uicb)
{
  printf("ENTERING UicbHalf_to_Uicb_GPU\n");
  
  //REMARK: Indexing of UicbHalf is poor form because fc is at start of Uicb call
  //but at end of UicbHalf call.
  
  int t = blockIdx.x*blkT + threadIdx.z;
  if (t < M_T) 
    {
      int g = threadIdx.y; //quadrature node on interface (M_G per interface)
      int dom = threadIdx.x; //dom side  (2 dominant sides per interface)
      scalar sum;
       printf("t=%d, g=%d, dom=%d\n", t, g, dom);
      //scalar sum = 0;
      //if (t >= M_T || g >= M_G || dom >= 2) return;
      for (int fc = 0; fc < N_F; fc++)
	{
	  /*
	    sum = UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0] + 
	    UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1];
	    ATOMICADD(&Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g], sum);
	  */
      
	  Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = 
	    UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0] + 
	    UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1];
      
	  printf("t=%d, g=%d, dom=%d, fc=%d: sumA=%f, sumB=%f, Uicb=%f\n",t,g,dom,fc,UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0],UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1],Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g]);
      
	}
    }
}

#endif

arch_global void Uhat_to_UhCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* Alt_FaceFromElem, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommon)
{
    /*!
    \brief Send each element's UhCommon contribution to the interface and store complete UhCommon.
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

  */
#ifdef USE_CPU
  {
  //Set UhCommon to zero
  for (int j = 0; j < M_T*M_G*N_F; j++)
    {
      UhCommon[j] = 0.0;
    }

  scalar* PushReco2Face = new scalar[N_N*M_G*N_F]; 
  for (int e = 0; e < N_E; e++)
    {
      //Matrix-multiply PSIxR_A_local's N_N superblocks by Uhat_A
      //to get the interface contributions
      for (int j = 0; j < N_N*M_G*N_F; j++)
	{
	  PushReco2Face[j] = 0.0;
	}
      for (int s = 0; s < (D+N_N); s++) //face of the present element, with room for D extra faces because sometimes faces are duplicated
	{
	  //int tGlo = Alt_FaceFromElem[e*N_N + s]; //this may not be the proper interface address to call
	  int tGlo = RecoPair[e*(D+N_N) + s]; 
	  //printf("Sending element %d, face %d contribution to t=%d\n",e, s, t);
	  for (int g = 0; g < M_G; g++) //quadrature node on the face
	    {
	      for (int fc = 0; fc < N_F; fc++) //field variable
		{
		  scalar sum = 0;
		  for (int k = 0; k < N_s; k++) //matrix-multiply PSIxR entries with element's Uhat
		    {
		      sum += PSIxR_elemwise[e*((D+N_N)*M_G*N_s) + s*(M_G*N_s) + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
		    }
		  //PushReco2Face[/*e*(N_N*M_G*N_F) + */s*(M_G*N_F) + g*N_F + fc] = sum;
		  UhCommon[tGlo*M_G*N_F + g*N_F + fc] += sum; //this += gets half of the common solution on the interface
		}
	    }
	}
    }
  /*
  //Now, communicate PushReco2Face to the appropriate interface addresses
  for (int s = 0; s < N_N; s++)
    {
      int t = Alt_FaceFromElem[e*N_N + s]; //this may not be the proper interface address to call
      for (int g = 0; g < M_G; g++)
	{
	  for (int fc = 0; fc < N_F; fc++)
	    {
	      UhCommon[t*M_G*N_F + g*N_F + fc] += PushReco2Face[]
	    }
	}
    }
  */
  delete[] PushReco2Face;
  }
#endif
#ifdef USE_GPU 
  {
    //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    /*
# if __CUDA_ARCH__>=200
    printf("%d \n", 5);
    
#endif
    */  
    /*
    int *jk;
    int value = 4;
    (*jk) = value*2;
    //printf("value = %d, address = %p\n",*jk,(void *)jk);
    */
    //printf("displayed from inside the kernel :\nvalue of jk = %d\nvalue of *jk = %p\n",jk,*jk);
    //doing this like kernels.cu/redistribute_sf
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int s = threadIdx.x; //side of local element
	int g = threadIdx.y; //quadrature node along the interface
	printf("Inside Uhat_to_UhCommon:e=%d, s=%d, g=%d\n",e,s,g);
	
	//printf("Not done with GPU version of Uhat_to_UhCommon_v2 yet\n");
	scalar sum = 0;
	int tGlo = RecoPair[e*(D+N_N) + s]; 
	/*
# if __CUDA_ARCH__>=200
	printf("Inside Uhat_to_UhCommon_v2: e=%d \n", e);
    
#endif  
	*/
	for (int fc = 0; fc < N_F; fc++)
	  {
	    for (int k = 0; k < N_s; k++)
	      {
		sum += PSIxR_elemwise[e*((D+N_N)*M_G*N_s) + s*(M_G*N_s) + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	      }
	    //	    UhCommon[tGlo*M_G*N_F + g*N_F + fc] += sum; //this += gets half of the common solution on the interface
	    UhCommon[tGlo*M_G*N_F + g*N_F + fc] = UhCommon[tGlo*M_G*N_F + g*N_F + fc] + sum; //this += gets half of the common solution on the interface
	    printf("\tsum to UhCommon(%d,%d,%d) = %f, corresponding UhCommon=%f\n",tGlo,g,fc,sum,UhCommon[tGlo*M_G*N_F + g*N_F + fc]); 
	    sum = 0; //re-zero the sum
	  }
      }
    
  }
#endif

  for (int t = 0; t < M_T; t++)
      {
	for (int g = 0; g < M_G; g++)
	  {
	    for (int fc = 0; fc < N_F; fc++)
	      {
		printf("UhCommon(%d,%d,%d) = %f\n",t,g,fc,UhCommon[t*M_G*N_F + g*N_F + fc]);
	      }
	  }
      }
}

//Send element contributions to interfaces for recovered solution storage (could be BR2 implementation if Averaging operator used for UhCommon)
arch_global void Uhat_to_UhCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf)
{
      /*!
    \brief Send each element's UhCommon contribution to the interface. Another routine must add the contributions
    \Uses element-wise storage of the PSIxR operator
    \param[in] N_E number of flesh elements on partition
    \param[in] Ne_AUG N_E+ghost element count
    \param[in] N_N faces per element
    \param[in] M_T interface count
    \param[in] M_G quadrature points per interface side
    \param[in] N_s solution points per element
    \param[in] RecoPair tells each element the indices of its perimeter interfaces 
    \param[in] BinarySideAddress tells each element whether it is A or B on an interface
    \param[in] PSIxR_elemwise the element-wise recovery(or averaging) operator
    \param[in] Uhat: element solution
    \param[out] UhCommonHalf: interface solution contribution, specifically at quadrature points

  */
  //printf("This function called only for CPU case and not ready\n");
#ifdef GPU
  printf("DO NOT CALL Uhat_to_UhCommonHalf on GPU architechture\n");
#endif
  //zero the UhCommonHalf array
  for (int j = 0; j < M_T*M_G*N_F*2; j++) {
  UhCommonHalf[j] = 0.0; }
  scalar sum = 0.0;
  for (int e = 0; e < Ne_AUG; e++) //element loop, over flesh and ghost elements
    {
      for (int s = 0; s < D+N_N; s++) //local face loop
	{
	  int tGlo = RecoPair[e*(D+N_N) + s];  //global face address corresponding to the local face
	  
	  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	  //printf("Uhat2UhCommonHalf: Sending contrbution from om=%d, s=%d to global interface %d, from side %d\n",e,s,tGlo,ArgSide);
	  for (int g = 0; g < M_G; g++) //quadrature node on the target global interface
	    {
	      for(int fc = 0; fc < N_F; fc++) //field variable
		{
		  sum = 0;
		  for (int k = 0; k < N_s; k++) //sum over the local element's DOF
		    {
		      sum += PSIxR_elemwise[e*((D+N_N)*M_G*N_s) + s*(M_G*N_s) + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
		    }
		  //gotta use += here for the sake of interface 0, which recieves zeros from a bunch of different elements. I'm still bitter about this.
		  UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide] += sum;
		  //		  printf("\tprintf(e=%d,s=%d,argside=%d); sum to UhCommonHalf(%d,%d,%d,side=%d) = %f; Actual UhCommonHalf = %f\n",e,s,ArgSide,tGlo,g,fc,ArgSide,sum,UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide]);
		}
	    } //end quadrature node loop
	} //end local face loop
    } //end element loop
}

//Sum components at each interface to get the common interface solution
arch_global void UhCommonHalf_to_UhCommon(int M_T, int M_G, scalar* UhCommonHalf, scalar* UhCommon)
{
   /*!
    \brief combine the UhCommonHalf contributions from each side of each interface to get UhCommon on the face
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_T number of interfaces
    \param[out] UhCommon: interface solution, specifically at quadrature points
   */
#ifdef GPU
  printf("DO NOT CALL UhCommonHalf_to_UhCommon on GPU architechture\n");
#endif
  //printf("This function called only for CPU and not ready\n");
  for (int t = 0; t < M_T; t++) //global interface address
    {
      for (int g = 0; g < M_G; g++) //quadrature node on the interface
	{
	  for (int fc = 0; fc < N_F; fc++) //field variable
	    {
	      UhCommon[t*M_G*N_F + g*N_F + fc] = UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0] + UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1];
	      //printf("Left-hand value = %f, right-hand = %f, UhCommon(%d,%d,%d) = %f\n", UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1], t,g,fc,UhCommon[t*M_G*N_F + g*N_F + fc]);
	    }
	} //end g loop
    } //end global interface loop
}

arch_global void Uhat_to_UhCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf)
{
      /*!
    \brief Send each element's UhCommon contribution to the interface. Another routine must add the contributions
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

  */
  //printf("This function called only for CPU case and not ready\n");
#ifdef GPU
  printf("DO NOT CALL Uhat_to_UhCommonHalf on GPU architechture\n");
#endif
  //zero the UhCommonHalf array
  for (int j = 0; j < M_T*M_G*N_F*2; j++) {
  UhCommonHalf[j] = 0.0; }
  scalar sum = 0.0;
  //for (int e = 0; e < N_E; e++) //element loop
  for (int e = 0; e < Ne_AUG; e++) //element loop, flesh and ghost elements
    {
      if (elemCmax[e] > Cthresh/*eps_gen*/ || sensor[e] > 0)
	{
      for (int s = 0; s < D+N_N; s++) //local face loop
	{
	  int tGlo = RecoPair[e*(D+N_N) + s];  //global face address corresponding to the local face
	  
	  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	  //printf("Uhat2UhCommonHalf: Sending contrbution from om=%d, s=%d to global interface %d, from side %d\n",e,s,tGlo,ArgSide);
	  for (int g = 0; g < M_G; g++) //quadrature node on the target global interface
	    {
	      for(int fc = 0; fc < N_F; fc++) //field variable
		{
		  sum = 0;
		  for (int k = 0; k < N_s; k++) //sum over the local element's DOF
		    {
		      sum += PSIxR_elemwise[e*((D+N_N)*M_G*N_s) + s*(M_G*N_s) + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
		    }
		  //gotta use += here for the sake of interface 0, which recieves zeros from a bunch of different elements. I'm still bitter about this.
		  UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide] += sum;
		  //		  printf("\tprintf(e=%d,s=%d,argside=%d); sum to UhCommonHalf(%d,%d,%d,side=%d) = %f; Actual UhCommonHalf = %f\n",e,s,ArgSide,tGlo,g,fc,ArgSide,sum,UhCommonHalf[tGlo*M_G*N_F*2 + g*N_F*2 + fc*2 + ArgSide]);
		}
	    } //end quadrature node loop
	} //end local face loop
	} //end if "need AD evaluation"
    } //end element loop
}

arch_global void UhCommonHalf_to_UhCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* UhCommonHalf, scalar* UhCommon)
{
   /*!
    \brief combine the UhCommonHalf contributions from each side of each interface to get UhCommon on the face
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_T number of interfaces
    \param[out] UhCommon: interface solution, specifically at quadrature points
   */
#ifdef GPU
  printf("DO NOT CALL UhCommonHalf_to_UhCommon on GPU architechture\n");
#endif
  //printf("This function called only for CPU and not ready\n");
  for (int t = 0; t < M_T; t++) //global interface address
    {
      //this routine needs to execute ONLY when a neighboring
      //element needs artificial dissipation
      int omA = BR2_Map[t*4 + 2*0 + 0];
      int omB = BR2_Map[t*4 + 2*1 + 0];
      if (elemCmax[omA] > Cthresh/*eps_gen*/ || sensor[omA] > 0 || elemCmax[omB] > Cthresh/*eps_gen*/ || sensor[omB] > 0)
	{
      for (int g = 0; g < M_G; g++) //quadrature node on the interface
	{
	  for (int fc = 0; fc < N_F; fc++) //field variable
	    {
	      UhCommon[t*M_G*N_F + g*N_F + fc] = UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0] + UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1];
	      //printf("Left-hand value = %f, right-hand = %f, UhCommon(%d,%d,%d) = %f\n", UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1], t,g,fc,UhCommon[t*M_G*N_F + g*N_F + fc]);
	    }
	} //end g loop
	} //end if "need AD evaluation" case
    } //end global interface loop
}

//Send element contributions to interfaces for ICB solution storage
arch_global void Uhat_to_UicbHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_biased_elemwise, scalar* Uhat, scalar* UicbHalf)
{
  /*!
    \brief Send each element's Uicb contribution to the interface. Another routine must reconcile the contributions
    \Uses element-wise storage of the biased PSIxR operator
    \param[in] N_E flush elements
    \param[in] N_E+Number of ghost elements
    \param[in] N_N faces per element
    \param[in] M_T total interface count
    \param[in] M_G quadrature points per side of interface
    \param[in] N_s solution points per element
    \param[in] RecoPair tells each element who its perimeter interfaces are
    \param[in] BinarySideAddress tells each element its position from interface perspective (A or B)
    \param[in] PSIxR_biased_elemwise: biased recovery operator
    \param[in] Uhat: element solution
    \param[out] UicbHalf: partial interface solution, specifically at quadrature points

  */
#ifdef GPU
  printf("DO NOT CALL Uhat_to_UicbHalf on GPU architechture\n");
#endif
  //zero the UicbHalf array
  //printf("From the top: imported M_T = %d\n", M_T); fflush(stdout);
  for (int j = 0; j < M_T*M_G*2*N_F*2; j++) {
  UicbHalf[j] = 0.0; }
  scalar sum = 0.0;
  for (int e = 0; e < Ne_AUG; e++) //element loop, over flesh and ghost elements
    {
      for (int s = 0; s < D+N_N; s++) //local face loop
	{
	  //printf("e=%d, s=%d:\n",e,s);
	  //printf("About to access RecoPair, input = %d\n",e*(D+N_N) + s); fflush(stdout);
	  int tGlo = RecoPair[e*(D+N_N) + s];  //global face address corresponding to the local face
	  //printf("RecoPair[%d] = %d\n",e*(D+N_N) + s, tGlo);
	  //printf("About to access BinarySideAddress, input = %d\n",e*(D+N_N) + s); fflush(stdout);
	  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	  //printf("BinarySideAddress[%d] = %d\n", e*(D+N_N) + s, ArgSide); fflush(stdout);
	  for (int g = 0; g < M_G; g++) //quadrature node on the target global interface
	    {
	      for (int dom = 0; dom < 2; dom++) //which icb solution is being populated (a-dominate or b-dominated)
		{
		  for(int fc = 0; fc < N_F; fc++) //field variable
		    {
		      sum = 0;
		      for (int k = 0; k < N_s; k++) //sum over the local element's DOF
			{
			  //sum += 1.0*Uhat[(e*N_F + fc)*N_s+k];
			  sum += PSIxR_biased_elemwise[e*((D+N_N)*2*M_G*N_s) + s*(2*M_G*N_s) + dom*M_G*N_s + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
			}
		      //gotta use += here for the sake of interface 0, which recieves zeros from a bunch of different elements. I'm still bitter about this.
		      UicbHalf[tGlo*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + ArgSide] += sum;
		      
		      //printf("\tprintf(e=%d,s=%d,dom=%d,argside=%d); sum to UicbHalf(%d,%d,dom=%d,fc=%d,side=%d) = %f; Actual UicbHalf = %f\n",e,s,dom,ArgSide,tGlo,g,dom,fc,ArgSide,sum,UicbHalf[tGlo*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + ArgSide]); fflush(stdout);
		    } //end field variable loop
		} //end dom loop
	    } //end quadrature node loop
	} //end local face loop
    } //end element loop
}

//Sum components at each interface to get the ICB solutions (2 per interface quadrature point, to be sent to Riemann solver)
arch_global void UicbHalf_to_Uicb(int M_T, int M_G, scalar* UicbHalf, scalar* Uicb)
{
   /*!
    \brief combine the UhCommonHalf contributions from each side of each interface to get UhCommon on the face
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_T number of interfaces
    \param[in] M_G quadrature points per interface
    \param[in] UicbHalf element contributions to icb solutions at each interface quadrature point
    \param[out] Uicb: interface solution, specifically at quadrature points
   */
  /*
    VERY IMPORTANT!! Organization of Uicb at end of this routine needs to match Marc's UintegF organization 
    used for calculation of interface fluxes (see evaluate_q in physics/physics.cu).
    UintegF[t*N_F*2*M_G + fc*2*M_G + side*M_G + g]
   */
#ifdef GPU
  printf("DO NOT CALL UicbHalf_to_Uicb on GPU architechture\n");
#endif
  //printf("This function called only for CPU and not ready\n");
  for (int t = 0; t < M_T; t++) //global interface address
    {
      for (int fc = 0; fc < N_F; fc++) //field
	{
	  for (int dom = 0; dom < 2; dom++) //the side of interface whose ICB solution is being fully summed
	    {
	      for (int g = 0; g < M_G; g++) //quadrature node on the interface
		{
		  //zero and one are for the contributing elements
		  //Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + 0*N_F*2 + fc*2 + dom] + UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + 1*N_F*2 + fc*2 + dom];
		  Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0] + UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1];
		  //printf("t=%d,f=%d,dom=%d,g=%d: left sum=%f, right sum=%f, Uicb=%f\n", t, fc, dom, g, UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 0], UicbHalf[t*M_G*2*N_F*2 + g*2*N_F*2 + dom*N_F*2 + fc*2 + 1], Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g]);
		} //end g loop
	    } //end dom loop
	} //end field variable loop
    } //end global interface loop
}

//Direct relay from DOF to interface ICB solutions
arch_global void Uhat_to_UicbDirect(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_biased_Global, scalar* Uhat, scalar* Uicb)
{
  /*!
    \brief Get the ICB solution along each interface
    \param[in] M_G number of quadrature nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_biased_Global: the biased recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] Uicb: interface ICB solution, 1 per side, specifically at quadrature points

  */

  //On a given interface:each row of UhCommon corresponds
  //to a quadrature point on the interface, column is field

#ifdef USE_CPU
  {
//Some variables that will be rewritten for each interface:
  scalar sum;
  int omA;
  int omB;
  scalar* RecoDOF = new scalar[N_F*2*N_s];
  
  for (int t = 0; t < M_T; t++){
    //Get the element addresses
    omA = BR2_Map[t*4];
    omB = BR2_Map[t*4 + 2];
    //Now, for this interface, load the DOF vectors
    //Element A DOF:
    
    for (int fc = 0; fc < N_F; fc++){
      for (int k = 0; k < N_s; k++){
	//U is in column-major form
	RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k]; }}
    
    //Element B DOF:
    for (int fc = 0; fc < N_F; fc++){
      for (int k = 0; k < N_s; k++){
	RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k]; }}
    
    //RecoDOF holds the DG DOG of the bi-element union
    //Now, perform the matrix multiplication (summation) to get UhCommon along the interface
    //UhCommon is never sent through BlasGemm, so I can organize it however I wish.
    
    for (int dom = 0; dom < 2; dom++){
      for (int fc = 0; fc < N_F; fc++){
 	for (int g = 0; g < M_G; g++){ 
	  //Initialize sum to zero:
	  sum = 0.0;
	  for (int k = 0; k < 2*N_s; k++){
	    //here, I am summing over the 2K DOF contained in the element union for the specific field
	    sum += PSIxR_biased_Global[t*(2*M_G*2*N_s) + dom*M_G*2*N_s + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	  } //end the k summation loop
	  //Relay the sum to global UhCommon storage:
	  Uicb[t*N_F*2*M_G + fc*2*M_G + dom*M_G + g] = sum;
	} //end the quadrature node (g) loop
      } //end the field variable (f) loop
    } //end "sominant side" loop
    
  } //end interface loop (t)
  }  //end (USE_CPU) section
#endif

#ifdef USE_GPU
#endif

}//End of subroutine

//Send element contributions to interfaces for common interface gradient.
arch_global void Uhat_to_GradCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf)
{
  /*!
    \brief send the element contribution to the common interface gradient
    \param[in] N_E flesh elements on partition
    \param[in] Ne_AUG N_E + number of ghost elements
    \param[in] N_N faces per element
    \param[in] M_T total interface count
    \param[in] M_G quadrature points per side of interface
    \param[in] N_s solution points per element
    \param[in] RecoPair tells each element who its perimeter interfaces are
    \param[in] BinarySideAddress tells each interface its orientation (0 or 1) relative to interface
    \param[in] SigFace_from_Uhat_elemwise the element-wise Uhat->interface gradient operator
    \param[in] Uhat the DOF
    \param[out] gradCommonHalf the element contributions to common interface gradient at each interface
   */
#ifdef USE_GPU
  printf("Don't call GradCommonHalf_to_GradCommon on GPU\n");
#endif

  //zero the gradCommonHalf array
  for (int j = 0; j < M_T*M_G*N_F*D*2; j++){
    gradCommonHalf[j] = 0.0; }

  
  for (int e = 0; e < Ne_AUG; e++) //loop over flesh and ghost elements
    {
      scalar sum = 0.0;
      //scalar FaceRelay[M_G*N_F*D];
      for (int s = 0; s < D+N_N; s++)
	{
	  
	  int tGlo = RecoPair[e*(D+N_N) + s]; 
	  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	  int fetch_op = e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s); //neighborhood in the SigFace operator
	  int fetch_U = e*N_F*N_s; //neighborhood in the U vector
	  int fetch_out = tGlo*(M_G*N_F*D*2) + ArgSide; //neighborhood in the gradCommonHalf vector
	  //printf("e=%d,s=%d,tGlo=%d,ArgSide=%d\n", e, s, tGlo, ArgSide);
	  int index = 0;
	  for (int g = 0; g < M_G; g++)
	    {
	      for(int fc = 0; fc < N_F; fc++)
		{
		  for (int a = 0; a < D; a++)
		    {
		      for (int k = 0; k < N_s; k++)
			{
			  //gradCommonHalf[fetch_out + g*(N_F*D*2) + fc*D*2 + a*2] += SigFace_from_Uhat_elemwise[fetch_op + g*(D*N_s) + a*N_s + k] * Uhat[fetch_U + fc*N_s + k];
			  sum += SigFace_from_Uhat_elemwise[fetch_op + g*(D*N_s) + a*N_s + k] * Uhat[fetch_U + fc*N_s + k];
			  //sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
			}
		      //gotta use += here because many elements send zeros to the zero interface
		      gradCommonHalf[fetch_out + g*(N_F*D*2) + fc*D*2 + a*2] += sum;
		      //gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide] += sum;
		      sum = 0;
		      //  printf("\tprintf(e=%d,s=%d,g=%d,a=%f,argside=%d); sum to gradCommonHalf(%d,%d,%d,a=%d,side=%d) = %f\n",e,s,g,a,ArgSide,tGlo,g,fc,a,ArgSide,gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide]);
		    } //end derivative component looop
		} //end field variable loop
	    } //end quadrature node loop
	  /*
	  //The relay to global storage location
	  for (int j = 0; j < M_G*N_F*D; j++){
	    gradCommonHalf[fetch_out + 2*j] = FaceRelay[j];}
	  */
	} //end local face loop
    } //end element loop
}

//Sum components for common gradient at each interface quadrature point
arch_global void GradCommonHalf_to_GradCommon(int M_T, int M_G, scalar* GradCommonHalf, scalar* GradCommon)
{
#ifdef USE_GPU
  printf("Don't call GradCommonHalf_to_GradCommon on GPU\n");
#endif
  for (int t = 0; t < M_T; t++)
    {
  for (int g = 0; g < M_G; g++)
    {
  for (int fc = 0; fc < N_F; fc++)
    {
  for( int a = 0; a < D; a++)
    {
  GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] = GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 0] + GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 1];
  //printf("GradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n", t,g,fc,a,GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
}
} //end field variable loop
} //end quadrature node loop
} //end interface loop
}

//For artificial dissipatio: repeat those subroutines, but with sensor-based considerations
arch_global void Uhat_to_GradCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf)
{
#ifdef USE_GPU
  printf("Don't call GradCommonHalf_to_GradCommon on GPU\n");
#endif

  //zero the gradCommonHalf array
  for (int j = 0; j < M_T*M_G*N_F*D*2; j++){
    gradCommonHalf[j] = 0.0; }

  
  //for (int e = 0; e < N_E; e++)
  for (int e = 0; e < Ne_AUG; e++) //loop over all flesh and ghost elements
    {
      if (elemCmax[e] > Cthresh/*eps_gen*/ || sensor[e] > 0)
	{
	  //	  printf("Executing Uhat_to_GradCommonHalf_AD for element %d; elemCmax=%f, sensor=%d\n", e, elemCmax[e], sensor[e]);
      scalar sum = 0.0;
      //scalar FaceRelay[M_G*N_F*D];
      for (int s = 0; s < D+N_N; s++)
	{
	  
	  int tGlo = RecoPair[e*(D+N_N) + s]; 
	  int ArgSide = BinarySideAddress[e*(D+N_N) + s]; //either zero or 1, tells the interface which side e is on
	  int fetch_op = e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s); //neighborhood in the SigFace operator
	  int fetch_U = e*N_F*N_s; //neighborhood in the U vector
	  int fetch_out = tGlo*(M_G*N_F*D*2) + ArgSide; //neighborhood in the gradCommonHalf vector
	  //printf("e=%d,s=%d,tGlo=%d,ArgSide=%d\n", e, s, tGlo, ArgSide);
	  int index = 0;
	  for (int g = 0; g < M_G; g++)
	    {
	      for(int fc = 0; fc < N_F; fc++)
		{
		  for (int a = 0; a < D; a++)
		    {
		      for (int k = 0; k < N_s; k++)
			{
			  //gradCommonHalf[fetch_out + g*(N_F*D*2) + fc*D*2 + a*2] += SigFace_from_Uhat_elemwise[fetch_op + g*(D*N_s) + a*N_s + k] * Uhat[fetch_U + fc*N_s + k];
			  sum += SigFace_from_Uhat_elemwise[fetch_op + g*(D*N_s) + a*N_s + k] * Uhat[fetch_U + fc*N_s + k];
			  //sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
			}
		      //gotta use += here because many elements send zeros to the zero interface
		      gradCommonHalf[fetch_out + g*(N_F*D*2) + fc*D*2 + a*2] += sum;
		      //gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide] += sum;
		      sum = 0;
		      //  printf("\tprintf(e=%d,s=%d,g=%d,a=%f,argside=%d); sum to gradCommonHalf(%d,%d,%d,a=%d,side=%d) = %f\n",e,s,g,a,ArgSide,tGlo,g,fc,a,ArgSide,gradCommonHalf[tGlo*(M_G*N_F*D*2) + g*(N_F*D*2) + fc*D*2 + a*2 + ArgSide]);
		    } //end derivative component looop
		} //end field variable loop
	    } //end quadrature node loop
	  /*
	  //The relay to global storage location
	  for (int j = 0; j < M_G*N_F*D; j++){
	    gradCommonHalf[fetch_out + 2*j] = FaceRelay[j];}
	  */
	} //end local face loop
	} //end if case for "The element needs AD"
    } //end element loop
}

//Sum components for common gradient at each interface quadrature point
arch_global void GradCommonHalf_to_GradCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* GradCommonHalf, scalar* GradCommon)
{
#ifdef USE_GPU
  printf("Don't call GradCommonHalf_to_GradCommon on GPU\n");
#endif
  for (int t = 0; t < M_T; t++)
    {
      //this routine needs to execute ONLY when a neighboring
      //element needs artificial dissipation
      int omA = BR2_Map[t*4 + 2*0 + 0];
      int omB = BR2_Map[t*4 + 2*1 + 0];
      if (elemCmax[omA] > Cthresh/*eps_gen*/ || sensor[omA] > 0 || elemCmax[omB] > Cthresh/*eps_gen*/ || sensor[omB] > 0)
	{
	  //	  printf("Executing GradCommonHalf_to_GradCommon_AD for t=%d: Cmax[%d]=%f, Cmax[%d]=%f, sensor[%d]=%d, sensor[%d]=%d\n",t, omA, elemCmax[omA], omB, elemCmax[omB], omA, sensor[omA], omB, sensor[omB]);
  for (int g = 0; g < M_G; g++)
    {
  for (int fc = 0; fc < N_F; fc++)
    {
  for( int a = 0; a < D; a++)
    {
  GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] = GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 0] + GradCommonHalf[t*M_G*N_F*D*2 + g*N_F*D*2 + fc*D*2 + a*2 + 1];
  //printf("GradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n", t,g,fc,a,GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
}
} //end field variable loop
} //end quadrature node loop
	} //end if loop for "the neighborhood needs AD"
} //end interface loop
}

arch_global void Uhat_to_GradCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, scalar*SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommon)
{
  /*!
    \brief Send each element's common gradient contribution to the interface and store complete UhCommon.
    \Uses element-wise storage of the PSIxR operator
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

  */
  #ifdef USE_CPU
  {
  //Set common gradient to zero
  for (int j = 0; j < M_T*M_G*N_F*D; j++) {
    gradCommon[j] = 0.0; }
  
  //The indexing used for SigFace_from_Uhat_elemwise needs to match the mapping applied in main.cc
  for (int e = 0; e < N_E; e++)
    {
      for (int s = 0; s < (D+N_N); s++) //face of the present element, with room for D extra faces because sometimes faces are duplicated
	{
	  //int tGlo = Alt_FaceFromElem[e*N_N + s]; //this may not be the proper interface address to call
	  int tGlo = RecoPair[e*(D+N_N) + s]; 
	  //printf("Sending element %d, face %d contribution to t=%d\n",e, s, t);
	  /*
	  for (int g = 0; g < M_G; g++) //quadrature node on the face
	    {
	      for (int fc = 0; fc < N_F; fc++) //field variable
		{
		  for (int a = 0; a < D; a++) //derivative direction
		    {
		      for (int k = 0; k < N_s; k++) //DOF of the current element
			{
			  gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += SigFace_from_Uhat_elemwise[e*((D+N_N)*D*M_G*N_s) + s*(D*M_G*N_s) + a*M_G*N_s + g*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
			}
		    }
		}
	    }
	  */
	  for (int g = 0; g < M_G; g++) //quadrature node on interface
	    {
	      for (int fc = 0; fc < N_F; fc++) //field variable
		{
		  for (int a = 0; a < D; a++) //derivative direction
		    {
		      for (int k = 0; k < N_s; k++)
			{
			  gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
			}
		    }
		}
	    }
	}
    }
  /*
  //Now, communicate PushReco2Face to the appropriate interface addresses
  for (int s = 0; s < N_N; s++)
    {
      int t = Alt_FaceFromElem[e*N_N + s]; //this may not be the proper interface address to call
      for (int g = 0; g < M_G; g++)
	{
	  for (int fc = 0; fc < N_F; fc++)
	    {
	      UhCommon[t*M_G*N_F + g*N_F + fc] += PushReco2Face[]
	    }
	}
    }
  */
  }
#endif

#ifdef USE_GPU 
  {
    for (int j = 0; j < M_T*M_G*N_F*D; j++){
      gradCommon[j] = 0.0; }

    scalar sum = 0;
    int e = blockIdx.x*blkE + threadIdx.z;
    int s = threadIdx.y;
    int g = threadIdx.x;
    int tGlo = RecoPair[e*(D+N_N) + s]; 
    for (int fc = 0; fc < N_F; fc++) {
      for (int a = 0; a < D; a++) {
	for (int k = 0; k < N_s; k++) {
	  sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k]; }
	gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += sum; //this += gets half of the common solution on the interface
	printf("\tprintf(e=%d,s=%d,g=%d); sum to gradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n",e,s,g,tGlo,g,fc,a,sum);
	printf("\tPresent gradCommon(t=%d,g=%d,f=%d,a=%d) = %f\n\n",tGlo,g,fc,a,gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a]);
      } //end a loop
    } //end fc loop
  /*
    //doing this like kernels.cu/redistribute_sf
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int s = threadIdx.x; //side of local element
	int g = threadIdx.y; //quadrature node along the interface
	//printf("Not done with GPU version of Uhat_to_UhCommon_v2 yet\n");
	scalar sum = 0;
	int tGlo = RecoPair[e*(D+N_N) + s]; 
	for (int fc = 0; fc < N_F; fc++) {
             for (int a = 0; a < D; a++) {
                 for (int k = 0; k < N_s; k++) {
  //gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += 
  //                                                                    SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * 
  //                                                                    Uhat[(e*N_F+fc)*N_s+k];
		sum += SigFace_from_Uhat_elemwise[e*((D+N_N)*M_G*D*N_s) + s*(M_G*D*N_s) + g*(D*N_s) + a*N_s + k] * Uhat[(e*N_F+fc)*N_s+k];
	      }
	    gradCommon[tGlo*(M_G*N_F*D) + g*(N_F*D) + fc*D + a] += sum; //this += gets half of the common solution on the interface
	    sum = 0; //re-zero the sum
}
      }
}
  */
  } //end GPU section
#endif
  


} 

arch_global void Uhat_to_CUGUF(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* SigFace_from_Uhat_ROG, scalar* Uhat, scalar* UhCommon, scalar* GradCommon){
  /*!
    \brief Get the interface UhCommon AND common gradient at each interface from neighboring DOF
    \param[in] M_G quadrature points per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] BR2_Map tells each interface who its neighboring elements are
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] SigFace_from_Uhat_ROG each interface's operator to get gradient from neighboring DOF (order is t-g-a-k)
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points
    \param]out] GradCommon: interface gradient, specifically at quadrature points
  */

  //On a given interface:each row of UhCommon corresponds
  //to a quadrature point on the interface, column is field

#ifdef USE_CPU
  {
    //Some variables that will be rewritten for each interface:
    scalar sum;
    int omA;
    int omB;
    scalar* RecoDOF = new scalar[N_F*2*N_s];
    
    for (int t = 0; t < M_T; t++){
      //Get the element addresses
      omA = BR2_Map[t*4 + 0];
      omB = BR2_Map[t*4 + 2];
      //Now, for this interface, load the DOF vectors
      //Element A DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k]; }}
      
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k]; }}
      
      //RecoDOF holds the DG DOG of the bi-element union
      //Now, perform the matrix multiplication (summation) to get UhCommon along the interface
      //UhCommon is never sent through BlasGemm, so I can organize it however I wish.
      
      for (int g = 0; g < M_G; g++){ 
	int row = t*M_G + g;
	for (int fc = 0; fc < N_F; fc++){
	  //Initialize sum to zero:
	  sum = 0.0;
	  
	  for (int k = 0; k < 2*N_s; k++){
	    //here, I am summing over the 2K DOF contained in the element union for the specific field
	    sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	  } //end the k summation loop
	  //Relay the sum to global UhCommon storage:
	  UhCommon[row*N_F + fc] = sum;
	  /*
	  if (sum > 3.1)
	    {
	      printf("Excessive UhCommon detected: t=%d, omA=%d, omB=%d, g=%d\n", t, omA, omB, g);
	    }
	  */
	} //end the field variable (f) loop
      } //end the quadrature node (g) loop
      //UhCommon has now been populated for the entire interface


      //While RecoDOF are available, get the common gradient as well
      for (int g = 0; g < M_G; g++){
	for (int fc = 0; fc < N_F; fc++){
	  for (int a = 0; a < D; a++){
	    sum = 0.0;
	    for (int k = 0; k < 2*N_s; k++){
	      sum += SigFace_from_Uhat_ROG[t*(M_G*D*2*N_s) + g*(D*2*N_s) + a*(2*N_s) + k] * RecoDOF[fc*2*N_s + k]; //need to build this reorganized SigFace array
	      //   printf("g=%d, fc=%d, a=%d, k=%d: SigFace_from_Uhat=%f, RecoDOF=%f\n", g, fc, a, k, SigFace_from_Uhat_ROG[t*(M_G*D*2*N_s) + g*(D*2*N_s) + a*(2*N_s) + k], RecoDOF[fc*2*N_s + k]);
	      //sum += SigFace_from_Uhat[t*(D*M_G*2*N_s) + a*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	    }
	    GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] = sum;
	    //    printf("Inside kernel: Grad(t=%d, g=%d, fc=%d, a=%d) = %f\n", t,g,fc,a,GradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a]);
	  }}}
      
      //GradCommon has now been populated along the interface
				     
    } //end interface loop (t)
  }  //end (USE_CPU) section
#endif
  
#ifdef USE_GPU
  {
    //Perform the above operation on GPU architecture instead.
    //This structuring taken from redistribuite_q in physics.cu, since that
    //function must run over interface nodes
    int omA;
    int omB;
    scalar* RecoDOF = new scalar[2*N_s*N_F];
    int t = blockIdx.x*blkT+threadIdx.z;
    if ( t < M_T){
      //Get the element addresses
      omA = BR2_Map[t*4 + 2*0 + 0];
      omB = BR2_Map[t*4 + 2*1 + 0];
      //Load the local DG degrees of freedom
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k];}}
     
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k];}}
      
      //threads must correspond to the dimBlock call in LUhat_to_UhCommon
      int g = threadIdx.x;
      int fc= threadIdx.y;
    
      //Marc does a summation in GPU architecture in redistribute_sf in physics.cu, so I guess it works
      scalar sum = 0.0;
      for (int k = 0; k < 2*N_s; k++){
	//here, I am summing over the 2K DOF contained in the element union for the specific field
	//UhCommon[row*N_F + f] += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[f*2*N_s + k];
	sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
      } //end the k summation loop
      //Relay the sum to global UhCommon storage:
      UhCommon[t*M_G*N_F + g*N_F + fc] = sum;
  
  
    } //End interface loop, t
  } //End (USE_GPU) section
#endif

}//End of subroutine

arch_global void Uhat_to_UhCommon(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* Uhat, scalar* UhCommon){
  /*!
    \brief Map solution to the faces
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

  */

  //On a given interface:each row of UhCommon corresponds
  //to a quadrature point on the interface, column is field

#ifdef USE_CPU
  {
//Some variables that will be rewritten for each interface:
  scalar sum;
  int omA;
  int omB;
  scalar* RecoDOF = new scalar[N_F*2*N_s];

  for (int t = 0; t < M_T; t++){
 //Get the element addresses
  omA = BR2_Map[t*4 + 2*0 + 0];
  omB = BR2_Map[t*4 + 2*1 + 0];
  //Now, for this interface, load the DOF vectors
  //Element A DOF:

for (int fc = 0; fc < N_F; fc++){
for (int k = 0; k < N_s; k++){
//U is in column-major form
RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k]; }}

//Element B DOF:
for (int fc = 0; fc < N_F; fc++){
for (int k = 0; k < N_s; k++){
RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k]; }}

//RecoDOF holds the DG DOG of the bi-element union
//Now, perform the matrix multiplication (summation) to get UhCommon along the interface
//UhCommon is never sent through BlasGemm, so I can organize it however I wish.

for (int g = 0; g < M_G; g++){ 
  int row = t*M_G + g;
  for (int fc = 0; fc < N_F; fc++){
  //Initialize sum to zero:
  sum = 0.0;

  for (int k = 0; k < 2*N_s; k++){
  //here, I am summing over the 2K DOF contained in the element union for the specific field
  sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
  } //end the k summation loop
  //Relay the sum to global UhCommon storage:
 UhCommon[row*N_F + fc] = sum;

} //end the field variable (f) loop
} //end the quadrature node (g) loop
//UhCommon has now been populated for the entire interface
} //end interface loop (t)
  }  //end (USE_CPU) section
#endif

#ifdef USE_GPU
  {
    //Perform the above operation on GPU architecture instead.
    //This structuring taken from redistribuite_q in physics.cu, since that
    //function must run over interface nodes
    int omA;
    int omB;
    scalar* RecoDOF = new scalar[2*N_s*N_F];
    int t = blockIdx.x*blkT+threadIdx.z;
    if ( t < M_T){
      //Get the element addresses
      omA = BR2_Map[t*4 + 2*0 + 0];
      omB = BR2_Map[t*4 + 2*1 + 0];
      //Load the local DG degrees of freedom
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k];}}
     
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k];}}
      
      //threads must correspond to the dimBlock call in LUhat_to_UhCommon
      int g = threadIdx.x;
      int fc= threadIdx.y;
    
      //Marc does a summation in GPU architecture in redistribute_sf in physics.cu, so I guess it works
      scalar sum = 0.0;
      for (int k = 0; k < 2*N_s; k++){
	//here, I am summing over the 2K DOF contained in the element union for the specific field
	//UhCommon[row*N_F + f] += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[f*2*N_s + k];
	sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
      } //end the k summation loop
      //Relay the sum to global UhCommon storage:
      UhCommon[t*M_G*N_F + g*N_F + fc] = sum;
  
  
    } //End interface loop, t
  } //End (USE_GPU) section
#endif

}//End of subroutine

arch_global void Uhat_to_CUGUF_OLD(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* serial_SigFace_from_DOF, scalar* Uhat, scalar* UhCommon, scalar* gradCommon){
  /*!
    \brief Map solution to the faces
    \param[in] M_s number of nodes per interface
    \param[in] M_T number of interfaces
    \param[in] N_s number of nodes per element
    \param[in] PSIxR_Global: the recovery operator for each interface
    \param[in] Uhat: element solution
    \param[out] UhCommon: interface solution, specifically at quadrature points

    \brief Take into account the geometry by multiplying with Jacobians
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[out] qJ q multiplied by Jacobian
    \param[in] q interface flux array
    \param[in] JF Jacobian of the faces
  */

  //On a given interface:each row of UhCommon corresponds
  //to a quadrature point on the interface, column is field

#ifdef USE_CPU
  {
    //Some variables that will be rewritten for each interface:
    scalar sum;
    int omA;
    int omB;
    scalar* RecoDOF = new scalar[N_F*2*N_s];
    
    for (int t = 0; t < M_T; t++){
      //Get the element addresses
      omA = BR2_Map[t*4 + 2*0 + 0];
      omB = BR2_Map[t*4 + 2*1 + 0];
      //Now, for this interface, load the DOF vectors
      //Element A DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k]; }}
      
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k]; }}
      
      //RecoDOF holds the DG DOG of the bi-element union
      //Now, perform the matrix multiplication (summation) to get UhCommon along the interface
      //UhCommon is never sent through BlasGemm, so I can organize it however I wish.
      
      for (int g = 0; g < M_G; g++){ 
	int row = t*M_G + g;
	for (int fc = 0; fc < N_F; fc++){
	  //Initialize sum to zero:
	  sum = 0.0;
	  
	  for (int k = 0; k < 2*N_s; k++){
	    //here, I am summing over the 2K DOF contained in the element union for the specific field
	    sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
	  } //end the k summation loop
	  //Relay the sum to global UhCommon storage:
	  UhCommon[row*N_F + fc] = sum;
	  
	} //end the field variable (f) loop
      } //end the quadrature node (g) loop
      //UhCommon has now been populated for the entire interface
      
      //Now, for the same interface, populate the gradient from DG DOF.
      //For now, don't want to do simple summation because of how I built sigFace_From_DOF.
      //May change in future, depending on the efficiency of this subroutine
      for (int g = 0; g < M_G; g++) {
	for (int f = 0; f < N_F; f++) {
	  for (int a = 0; a < D; a++){
	    gradCommon[t*M_G*N_F*D + g*N_F*D + f*D + a] = 0.0; }}}
      
      //With gradCommon zeroed, begin summing
      for (int g = 0; g < M_G; g++) {
	for (int fc = 0; fc < N_F; fc++) {
	  for (int a = 0; a < D; a++) {
	    for (int k = 0; k < 2*N_s; k++) {
	      gradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] += serial_SigFace_from_DOF[t*D*M_G*2*N_s + 
										     a*M_G*2*N_s + 
										     g*2*N_s +
										     k] * RecoDOF[fc*2*N_s + k];
	    } //end k loop
	  } //end a loop
	} //end f loop
      } //end g loop
      //gradCommon has now  been populated for the interface
      
      /*
      for (int a = 0; a < D; a++) {
	for (int g = 0; g < M_G; g++) {
	  sum = 0.0;
	  for (int k = 0; k < 2*N_s; k++)
	    {
	      sum += serial_SigFace_from_DOF[t*D*M_G*2*N_s + 
					       a*M_G*2*N_s + 
					           g*2*N_s +
					                 k] * RecoDOF[fc*2*N_s + k];
	    }
	  gradCommon[t*M_G*N_F*D + g*N_F*D + f*D + a]
      */
    } //end interface loop (t)
  }  //end (USE_CPU) section
#endif

#ifdef USE_GPU
  {
    //Perform the above operation on GPU architecture instead.
    //This structuring taken from redistribuite_q in physics.cu, since that
    //function must run over interface nodes
    int omA;
    int omB;
    scalar* RecoDOF = new scalar[2*N_s*N_F];
    int t = blockIdx.x*blkT+threadIdx.z;
    if ( t < M_T){
      //Get the element addresses
      omA = BR2_Map[t*4 + 2*0 + 0];
      omB = BR2_Map[t*4 + 2*1 + 0];
      //Load the local DG degrees of freedom
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  //U is in column-major form
	  RecoDOF[fc*2*N_s + k] = Uhat[(omA*N_F+fc)*N_s+k];}}
     
      //Element B DOF:
      for (int fc = 0; fc < N_F; fc++){
	for (int k = 0; k < N_s; k++){
	  RecoDOF[fc*2*N_s + N_s + k] = Uhat[(omB*N_F+fc)*N_s + k];}}
      
      //zero the gradU array before threading g, fc
      for (int g = 0; g < M_G; g++) {
	for (int f = 0; f < N_F; f++) {
	  for (int a = 0; a < D; a++){
	    gradCommon[t*M_G*N_F*D + g*N_F*D + f*D + a] = 0.0; }}}

      //threads must correspond to the dimBlock call in LUhat_to_UhCommon
      int g = threadIdx.x;
      int fc= threadIdx.y;
      /*
      //Marc does a summation in GPU architecture in redistribute_sf in physics.cu, so I guess it works
      scalar sum = 0.0;
      for (int k = 0; k < 2*N_s; k++){
	//here, I am summing over the 2K DOF contained in the element union for the specific field
	//UhCommon[row*N_F + f] += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[f*2*N_s + k];
	sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * RecoDOF[fc*2*N_s + k];
      } //end the k summation loop
      //Relay the sum to global UhCommon storage:
      UhCommon[t*M_G*N_F + g*N_F + fc] = sum;
  

      //Repeat for gradCommon
      for (int a = 0; a < D; a++) {
	for (int k = 0; k < 2*N_s; k++) {
	  gradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] += serial_SigFace_from_DOF[t*D*M_G*2*N_s + 
										 a*M_G*2*N_s + 
										 g*2*N_s +
										 k] * RecoDOF[fc*2*N_s + k];
	} //gradCommon k loop
										 
      */
      //trying the above without using RecoDOF, because I fear a race condition.
      //This did not help the situation
      scalar sum = 0;
      for (int k = 0; k < N_s; k++)
	{
	  sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + k] * Uhat[(omA*N_F+fc)*N_s + k];
	}
      for (int k = 0; k < N_s; k++)
	{
	  sum += PSIxR_Global[t*(M_G*2*N_s) + g*(2*N_s) + N_s + k] * Uhat[(omB*N_F+fc)*N_s + k];
	} //end k summation
       //Relay the sum to global UhCommon storage:
      UhCommon[t*M_G*N_F + g*N_F + fc] = sum;
      //Repeat for gradCommon
      for (int a = 0; a < D; a++) {
	for (int k = 0; k < N_s; k++) {
	  gradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] += serial_SigFace_from_DOF[t*D*M_G*2*N_s + 
										 a*M_G*2*N_s + 
										 g*2*N_s +
										 k] * Uhat[(omA*N_F+fc)*N_s + k];
	}
	for (int k = 0; k < N_s; k++) {
	  gradCommon[t*M_G*N_F*D + g*N_F*D + fc*D + a] += serial_SigFace_from_DOF[t*D*M_G*2*N_s + 
										 a*M_G*2*N_s + 
										 g*2*N_s +
										 N_s + k] * Uhat[(omB*N_F+fc)*N_s + k];
	} //gradcommon k loop

      } //gradCommon a loop
  
      //g,fc are threaded

    } //End interface loop, t
  } //End (USE_GPU) section
#endif

}//End of subroutine

//Get element contribution to its own sigma approximation (sigma=gradient of U)
arch_global void Uhat_to_Sigma(int N_E, int N_G, int N_s, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el)
{
   /*!
    \bried Added 05/27/2017 by PEJ. Partially populates auxiliary variable from element DOF
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G quadrature nodes per element
    \param[in] serial_MSigVol matrix that grabs auxiliary variabl contribution from DOF
    \param[in] Uhat DG DOF
    \param[out] gradU_el auxiliary variable over element interiors. Not yet completely populated
  */
#ifdef USE_CPU
  {
    for (int e = 0; e < N_E; e++) //element in which we are summing gradient
      {
	scalar sum;
	for (int a = 0; a < D; a++) //direction of gradient component
	  {
	    for (int g = 0; g < N_G; g++) //quadrature point in the element
	      {
		//some indexing parameters:
		int MegaSlot_M =    e*D*N_G*N_s + a*N_G*N_s + g*N_s;
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		
		for (int f = 0; f < N_F; f++) //field variable whose gradient we are calculating
		  {
		    int MegaSlot_U =    e*N_F*N_s   + f*N_s;
		    sum = 0;
		    for (int jj = 0; jj < N_s; jj++) //sum over the element's DOF
		      {
			//sum += serial_MSigVol[MegaSlot_M + jj] * Uhat[MegaSlot_U + jj];
			sum += serial_MSigVol[e*D*N_G*N_s + a*N_G*N_s + g*N_s + jj] * Uhat[e*N_F*N_s   + f*N_s + jj];
		      }
		    //relay the sum to global storage:
		    //gradU_el[MegaSlot_grad + f] = sum;
		    gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] = sum;
		    //		    printf("Partially populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f, sum=%f\n", e,a,g,f,gradU_el[MegaSlot_grad + f],sum);
		  } //end f<N_F loop
	      } //end g<N_G loop
	  } //end a<D loop
      } //end e < N_E loop
  } //end of (USE_CPU) section
#endif
#ifdef USE_GPU
  {
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int a = threadIdx.x;
	int g = threadIdx.y;
	scalar sum = 0;
	//int MegaSlot_M =    e*D*N_G*N_s + a*N_G*N_s + g*N_s;
	//int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int f = 0; f < N_F; f++)
	  {
	    //int MegaSlot_U =    e*N_F*N_s   + f*N_s;
	    //sum = 0;
	    for (int jj = 0; jj < N_s; jj++)
	      {
		sum += serial_MSigVol[e*D*N_G*N_s + a*N_G*N_s + g*N_s + jj] * Uhat[e*N_F*N_s   + f*N_s + jj];
		//		sum += serial_MSigVol[MegaSlot_M + jj] * Uhat[MegaSlot_U + jj];
	      }
	    gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] = sum;
	    //gradU_el[MegaSlot_grad + f] = sum;
	    //ATOMICADD(&gradU_el[MegaSlot_grad + f] , sum);
	    //	    printf("\t\tpresent sum=%f\n",sum);
	    sum = 0;
	    //    printf("Partially populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f, sum=%f\n", e,a,g,f,gradU_el[MegaSlot_grad + f],sum);
	  }
      }
  } //end of (USE_GPU) section
#endif
} //end of suroutine

arch_global void Uhat_to_Sigma_AD(int N_E, int N_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el)
{
   /*!
    \bried Added 05/27/2017 by PEJ. Partially populates auxiliary variable from element DOF
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] N_G quadrature nodes per element
    \param[in] serial_MSigVol matrix that grabs auxiliary variabl contribution from DOF
    \param[in] Uhat DG DOF
    \param[out] gradU_el auxiliary variable over element interiors. Not yet completely populated
  */
#ifdef USE_CPU
  {
    for (int e = 0; e < N_E; e++) //element in which we are summing gradient
      {
	if (elemCmax[e] > Cthresh/*eps_gen*/ || sensor[e] > 0)
	{
	  //	  printf("Executing Uhat_to_GradCommonHalf_AD for element %d; elemCmax=%f, sensor=%d\n", e, elemCmax[e], sensor[e]);
	scalar sum;
	for (int a = 0; a < D; a++) //direction of gradient component
	  {
	    for (int g = 0; g < N_G; g++) //quadrature point in the element
	      {
		//some indexing parameters:
		int MegaSlot_M =    e*D*N_G*N_s + a*N_G*N_s + g*N_s;
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		
		for (int f = 0; f < N_F; f++) //field variable whose gradient we are calculating
		  {
		    int MegaSlot_U =    e*N_F*N_s   + f*N_s;
		    sum = 0;
		    for (int jj = 0; jj < N_s; jj++) //sum over the element's DOF
		      {
			//sum += serial_MSigVol[MegaSlot_M + jj] * Uhat[MegaSlot_U + jj];
			sum += serial_MSigVol[e*D*N_G*N_s + a*N_G*N_s + g*N_s + jj] * Uhat[e*N_F*N_s   + f*N_s + jj];
		      }
		    //relay the sum to global storage:
		    //gradU_el[MegaSlot_grad + f] = sum;
		    gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] = sum;
		    //		    printf("Partially populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f, sum=%f\n", e,a,g,f,gradU_el[MegaSlot_grad + f],sum);
		  } //end f<N_F loop
	      } //end g<N_G loop
	  } //end a<D loop
	} //end if "need AD evalueation" case
      } //end e < N_E loop
  } //end of (USE_CPU) section
#endif
#ifdef USE_GPU
  {
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int a = threadIdx.x;
	int g = threadIdx.y;
	scalar sum = 0;
	//int MegaSlot_M =    e*D*N_G*N_s + a*N_G*N_s + g*N_s;
	//int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int f = 0; f < N_F; f++)
	  {
	    //int MegaSlot_U =    e*N_F*N_s   + f*N_s;
	    //sum = 0;
	    for (int jj = 0; jj < N_s; jj++)
	      {
		sum += serial_MSigVol[e*D*N_G*N_s + a*N_G*N_s + g*N_s + jj] * Uhat[e*N_F*N_s   + f*N_s + jj];
		//		sum += serial_MSigVol[MegaSlot_M + jj] * Uhat[MegaSlot_U + jj];
	      }
	    gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] = sum;
	    //gradU_el[MegaSlot_grad + f] = sum;
	    //ATOMICADD(&gradU_el[MegaSlot_grad + f] , sum);
	    //	    printf("\t\tpresent sum=%f\n",sum);
	    sum = 0;
	    //    printf("Partially populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f, sum=%f\n", e,a,g,f,gradU_el[MegaSlot_grad + f],sum);
	  }
      }
  } //end of (USE_GPU) section
#endif
} //end of suroutine

//Send interface contirbutions to each element's sigma approximation
arch_global void UhCommon_to_Sigma(int N_E, int N_G, int N_N, int M_G, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el)
{
   /*!
    \bried Added 05/27/2017 by PEJ. Partially populates auxiliary variable from UhCommon on surrounding faces
    \param[in] N_N faces per element
    \param[in] N_E number of elements
    \param[in] N_G quadrature nodes per element
    \param[in] M_G quadrature nodes per interface
    \param[in] serial_MSigSurf matrix that grabs auxiliary variabl contribution from UhCommon on interfaces
    \param[in] UhCommon common value distribution on interfaces
    \param[in] gradU_el auxiliary variable over element interiors. Not yet completely populated
    \param[out] gradU_el auxiliary variable over element interiors. Fully populated
  */

  //VERY IMPORTANT: gradU_el has already been partially populated by Uhat_to_Sigma,
  //so this function must add to existing gradU_el values. That's why I don't
  //zero before adding to the summation
  
#ifdef USE_CPU
  {
    //scalar sum;
    for (int e = 0; e < N_E; e++)
      {
	for (int a = 0; a < D; a++) //gradient component
	  {
	    for (int g = 0; g < N_G; g++) //element quadrature point
	      {
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		//Because of how UhCommon is organized (last index is f)
		//I want the f loop to be inside the H loop; this is tricky
		//because the inner f loop is going to be altering
		//different values of gradU_el. Consequently,
		//can't populate sum locally before relaying, as I prefer to do.
		/*
		for (int f = 0; f < Nfield; f++) //field variable whose gradient we are calculating
		  {
		*/
		for (int H = 0; H < N_N; H++) //the contributing interface
		  {
		    int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
		    int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G;
		    int MegaSlot_U = tGlo * (M_G*N_F);
		    for (int jj = 0; jj < M_G; jj++) //summation over the global interface's quadrature nodes
		      {
			for (int f = 0; f < N_F; f++)
			  {
			    gradU_el[MegaSlot_grad + f] += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
			    /*
			    sum += serial_MSigSurf[e*D*N_G*N_N*M_G + 
						   a*N_G*N_N*M_G + 
						   g*N_N*M_G + 
						   H*M_G + 
						   jj] * 
			      UhCommon[tGlo*M_G*N_F + jj*N_F + f];
			    */
			  }
		      } 
		    //printf("Fully populated gradU_el(e=%d,g=%d,f=%d) = %f\n", e,g,f,gradU_el[MegaSlot_grad + f]);
		    //gotta use += here because volume contribution is already stored in gradU_el
		    //gradU_el[MegaSlot_grad + f] += sum;
		  }
	      }//end g < N_G quadature node loop
	  } //end a < D component loop
      } //end e< N_E element loop


    /*
    for (int e = 0; e < N_E; e++)
      {
	for (int a = 0; a < D; a++)
	  {
	    for (int g = 0; g < N_G; g++)
	      {
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		for (int f = 0; f < N_F; f++)
		  {
		    printf("Fully populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f\n", e,a,g,f,gradU_el[MegaSlot_grad + f]);
		  }
	      }
	  }
      }
*/

  } //end CPU routine
#endif
#ifdef USE_GPU
  {
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int a = threadIdx.x;
	int g = threadIdx.y;
	scalar sum = 0.0;
	//int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int f = 0; f < N_F; f++)
	  {
	    for (int H = 0; H < N_N; H++) //the contributing interface
	      {
		int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
		//int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G; //where to look in the serial_SigSurf matrix
		//int MegaSlot_U = tGlo * (M_G*N_F); //where to look in UhCommon
		
		for (int jj = 0; jj < M_G; jj++) //summation over the face's values at gaussian quadrature points
		  {
		    sum += serial_MSigSurf[e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G + jj] * UhCommon[tGlo * (M_G*N_F) + jj*N_F + f];
		    //sum += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
		  }
	      }
	    //gradU_el[MegaSlot_grad + f] += sum;
	    //ATOMICADD(&gradU_el[MegaSlot_grad + f] , sum);
	    ATOMICADD(&gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] , sum);
	    sum = 0;
	    //	    printf("Fully populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f\n", e,a,g,f,gradU_el[MegaSlot_grad + f]);
	  } //end field variable loop
	/*
	int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int H = 0; H < N_N; H++) //the contributing interface
	  {
	    scalar sum = 0;
	    int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
	    int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G;
	    int MegaSlot_U = tGlo * (M_G*N_F);
	    
	    for (int jj = 0; jj < M_G; jj++) //summation over the global interface's quadrature nodes
	      {
		for (int f = 0; f < N_F; f++)
		  {
		    gradU_el[MegaSlot_grad + f] += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
		  }
	      } //end face summation loop
	    
	  } 
	*/
	
      } //end e< N_E element loop
  } //end GPU routine
#endif
}

arch_global void UhCommon_to_Sigma_AD(int N_E, int N_G, int N_N, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el)
{
   /*!
    \bried Added 05/27/2017 by PEJ. Partially populates auxiliary variable from UhCommon on surrounding faces
    \param[in] N_N faces per element
    \param[in] N_E number of elements
    \param[in] N_G quadrature nodes per element
    \param[in] M_G quadrature nodes per interface
    \param[in] serial_MSigSurf matrix that grabs auxiliary variabl contribution from UhCommon on interfaces
    \param[in] UhCommon common value distribution on interfaces
    \param[in] gradU_el auxiliary variable over element interiors. Not yet completely populated
    \param[out] gradU_el auxiliary variable over element interiors. Fully populated
  */

  //VERY IMPORTANT: gradU_el has already been partially populated by Uhat_to_Sigma,
  //so this function must add to existing gradU_el values. That's why I don't
  //zero before adding to the summation
  
#ifdef USE_CPU
  {
    //scalar sum;
    for (int e = 0; e < N_E; e++)
      {
	if (elemCmax[e] > Cthresh/*eps_gen*/ || sensor[e] > 0)
	{
	  //	  printf("Executing Uhat_to_GradCommonHalf_AD for element %d; elemCmax=%f, sensor=%d\n", e, elemCmax[e], sensor[e]);
	for (int a = 0; a < D; a++) //gradient component
	  {
	    for (int g = 0; g < N_G; g++) //element quadrature point
	      {
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		//Because of how UhCommon is organized (last index is f)
		//I want the f loop to be inside the H loop; this is tricky
		//because the inner f loop is going to be altering
		//different values of gradU_el. Consequently,
		//can't populate sum locally before relaying, as I prefer to do.
		/*
		for (int f = 0; f < Nfield; f++) //field variable whose gradient we are calculating
		  {
		*/
		for (int H = 0; H < N_N; H++) //the contributing interface
		  {
		    int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
		    int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G;
		    int MegaSlot_U = tGlo * (M_G*N_F);
		    for (int jj = 0; jj < M_G; jj++) //summation over the global interface's quadrature nodes
		      {
			for (int f = 0; f < N_F; f++)
			  {
			    gradU_el[MegaSlot_grad + f] += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
			    /*
			    sum += serial_MSigSurf[e*D*N_G*N_N*M_G + 
						   a*N_G*N_N*M_G + 
						   g*N_N*M_G + 
						   H*M_G + 
						   jj] * 
			      UhCommon[tGlo*M_G*N_F + jj*N_F + f];
			    */
			  }
		      } 
		    //printf("Fully populated gradU_el(e=%d,g=%d,f=%d) = %f\n", e,g,f,gradU_el[MegaSlot_grad + f]);
		    //gotta use += here because volume contribution is already stored in gradU_el
		    //gradU_el[MegaSlot_grad + f] += sum;
		  }
	      }//end g < N_G quadature node loop
	  } //end a < D component loop
	} //end if "need AD evaluation" case
      } //end e< N_E element loop


    /*
    for (int e = 0; e < N_E; e++)
      {
	for (int a = 0; a < D; a++)
	  {
	    for (int g = 0; g < N_G; g++)
	      {
		int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
		for (int f = 0; f < N_F; f++)
		  {
		    printf("Fully populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f\n", e,a,g,f,gradU_el[MegaSlot_grad + f]);
		  }
	      }
	  }
      }
*/

  } //end CPU routine
#endif
#ifdef USE_GPU
  {
    int e = blockIdx.x*blkE + threadIdx.z;
    if (e < N_E)
      {
	int a = threadIdx.x;
	int g = threadIdx.y;
	scalar sum = 0.0;
	//int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int f = 0; f < N_F; f++)
	  {
	    for (int H = 0; H < N_N; H++) //the contributing interface
	      {
		int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
		//int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G; //where to look in the serial_SigSurf matrix
		//int MegaSlot_U = tGlo * (M_G*N_F); //where to look in UhCommon
		
		for (int jj = 0; jj < M_G; jj++) //summation over the face's values at gaussian quadrature points
		  {
		    sum += serial_MSigSurf[e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G + jj] * UhCommon[tGlo * (M_G*N_F) + jj*N_F + f];
		    //sum += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
		  }
	      }
	    //gradU_el[MegaSlot_grad + f] += sum;
	    //ATOMICADD(&gradU_el[MegaSlot_grad + f] , sum);
	    ATOMICADD(&gradU_el[e*D*N_G*N_F + a*N_G*N_F + g*N_F + f] , sum);
	    sum = 0;
	    //	    printf("Fully populated gradU_el(e=%d,a=%d,g=%d,f=%d) = %f\n", e,a,g,f,gradU_el[MegaSlot_grad + f]);
	  } //end field variable loop
	/*
	int MegaSlot_grad = e*D*N_G*N_F + a*N_G*N_F + g*N_F;
	for (int H = 0; H < N_N; H++) //the contributing interface
	  {
	    scalar sum = 0;
	    int tGlo = Alt_FaceFromElem[e*N_N + H]; //gets global interface address of element's face H
	    int MegaSlot_M =    e * (D*N_G*N_N*M_G) + a * (N_G*N_N*M_G) + g * (N_N*M_G) + H*M_G;
	    int MegaSlot_U = tGlo * (M_G*N_F);
	    
	    for (int jj = 0; jj < M_G; jj++) //summation over the global interface's quadrature nodes
	      {
		for (int f = 0; f < N_F; f++)
		  {
		    gradU_el[MegaSlot_grad + f] += serial_MSigSurf[MegaSlot_M + jj] * UhCommon[MegaSlot_U + jj*N_F + f];
		  }
	      } //end face summation loop
	    
	  } 
	*/
	
      } //end e< N_E element loop
  } //end GPU routine
#endif
}


//Correct UhCommon for boundary interfaces
arch_global void CorrectUhCommonBC(int M_T, int M_G, int M_B, int* boundaryMap, scalar* UintegF, scalar* UhCommon)
{
  for (int tB = 0; tB < M_B; tB++)
    {
      int tG = boundaryMap[tB]; //global interface address
      //printf("In CorrectUhCommonBC: Modifying UhCommon distrubution along boundary interface %d, actual interface %d\n",tB,tG);
      for (int g = 0; g < M_G; g++)
	{
	  for (int fc = 0; fc < N_F; fc++)
	    {
	      //Configuration 1: Set UhCommon to the exterior (dirichlet) value.
	      //The 1*MG means we take this from the exterior side (that's what it means on boundary interfaces)
	      /*
	      printf("tG=%d,g=%d,fc=%d: interior UintegF = %f, current UhCommon = %f, Dirichlet = %f\n",tG,g,fc,
		     UintegF[tG*N_F*2*M_G + fc*2*M_G + 0*M_G + g],
		     UhCommon[tG*M_G*N_F + g*N_F + fc],
		     UintegF[tG*N_F*2*M_G + fc*2*M_G + 1*M_G + g]);
	      */
	      
	      // UhCommon[tG*M_G*N_F + g*N_F + fc] = UintegF[tG*N_F*2*M_G + fc*2*M_G + 1*M_G + g];
	      
	      //Configuration 2: Set UhCommon to average of interior and exterior traces.
	      //This is for when the exterior value is enginnered to achieve a special effect
	      //at the interface, like when velocity is reversed to get a wall BC. In that
	      //case, it makes more sense to take the average
	      UhCommon[tG*M_G*N_F + g*N_F + fc] = 0.5*(UintegF[tG*N_F*2*M_G + fc*2*M_G + 0*M_G + g] + UintegF[tG*N_F*2*M_G + fc*2*M_G + 1*M_G + g]);
	    }
	}
    }
		
}

//Send Uhat contrubution to gradient calculation for boundary interfaces
arch_global void Uhat_to_GradCommonBC(int M_B, int M_G, int N_s, int N_N, int* boundaryMap, int* BR2_Map, scalar* serial_Uhat2GradBC, scalar* Uhat, int* RecoPair, int* Alt_FaceFromElem, scalar* UhCommon, scalar* serial_UhCommon2GradBC, scalar* GradCommon)
{
  /*
    Very nasty subroutine. Looks at the UhCommon distribtion around each boundary element,
    in addition to the DOF of the element, and uses that information
    to build the gradient approximation at the boundary.
    
    GradCommon entries that are touched in this routine should
    not match any interior interface entries.
   */
  for (int tB = 0; tB < M_B; tB++) //the boundary interface
    {
      int tG = boundaryMap[tB]; //global interface address
      int omA = BR2_Map[tG*4 + 2*0 + 0]; //the element we are grabbing Uhat from
      for (int g = 0; g < M_G; g++) //quadrature point on the boundary interface
	{
	  for (int fc = 0; fc < N_F; fc++) //field
	    {
	      for (int a = 0; a < D; a++) //gradient direction
		{
		  scalar sum = 0;
		  for (int k = 0; k < N_s; k++) //summation over the element's DOF
		    {
		      sum += serial_Uhat2GradBC[tB*D*M_G*N_s + a*M_G*N_s + g*N_s  + k] * Uhat[(omA*N_F+fc)*N_s+k];
		    }
		  GradCommon[tG*M_G*N_F*D + g*N_F*D + fc*D + a] = sum;
		}
	    }
	} //end g loop

      //Since knowledge of omA is necessary, I'm also doing the interface contribution inside
      //this element
      scalar UhCommon_Element[N_N*M_G*N_F];
      for (int s = 0; s < N_N; s++) //I think that when true boundary is involved, only need N_N s entries
	{
	  //int tSide = RecoPair[omA*(D+N_N) + s]; //interface we are loading information from
	  int tSide = Alt_FaceFromElem[omA*N_N + s]; //interface we are loading information from
	  for (int j = 0; j < M_G*N_F; j++)
	    {
	      UhCommon_Element[s*M_G*N_F + j] = UhCommon[tSide*M_G*N_F + j];
	    }
	}

      //Now, we know the UhCommon distribution around the entire element. Use it to 
      //finish populating GradCommon
      for (int s = 0; s < N_N; s++) //contributing side of the boundary element
	{
	  for (int g = 0; g < M_G; g++) //quadrature node on boundary interface
	    {
	      for (int fc = 0; fc < N_F; fc++) //field variable
		{
		  for (int a = 0; a < D; a++) //gradient component
		    {
		      scalar sum = 0;
		      for (int j = 0; j < M_G; j++) //contributing quadrature point from face s of omA
			{
			  sum += serial_UhCommon2GradBC[tB*D*M_G*(N_N*M_G) + a*M_G*(N_N*M_G) + g*(N_N*M_G) + s*M_G + j] * UhCommon_Element[s*M_G*N_F + j*N_F + fc];
			}
		      GradCommon[tG*M_G*N_F*D + g*N_F*D + fc*D + a] += sum;
		    }
		}
	    } //end g loop
	} //end s loop on faces of boundary element
      /*
      printf("Boundary Interface %d, Global interface %d:\n",tB,tG);
      for (int g = 0; g < M_G; g++)
	{
	  for (int fc = 0; fc < N_F; fc++)
	    {
	      for (int a = 0; a < D; a++)
		{
		  printf("\t(g=%d,fc=%d,dir=%d) derivative = %f\n",g,fc,a,GradCommon[tG*M_G*N_F*D + g*N_F*D + fc*D + a]);
		}
	    }
	}
      */
    } //end loop on all boundary faces
  
}

arch_global void BuildBC_UintegF(int M_B, int M_G, int M_s, int* boundaryMap, scalar* psi, scalar* UF, scalar* UintegF)
{
  /*!
    \brief routine to get the interface solution along boundary quadrature points when using ICBN procedure
    \param[in] M_B total count of boundary interfaces
    \param[in] M_G face quadrature points per side of interface
    \param[in] M_s supported solution nodes per side of interface
    \param[in] boundaryMap gets the global interface address from the boundary address
    \param[in] psi the face-supported basis functions
    \param[in] UF solution values at face solution nodes
    \param[out] UintegF the solution at the face quadrature points.
  */
  //MapToFace has already been executed, so we have UF along each boundary interface
  for (int tBC = 0; tBC < M_B; tBC++)
    {
      int tGlo = boundaryMap[tBC];
      for (int fc = 0; fc < N_F; fc++) //field variable on the interface
	{
	  for (int d = 0; d < 2; d++) //side of the interface (still need two arguments for Riemann solver)
	    {
	      for (int g = 0; g < M_G; g++) //quadrature node on the interface
		{
		  //Now, sum over the face-supported basis functions
		  scalar sum = 0;
		  for (int m = 0; m < M_s; m++) //sum over supported basis functions
		    {
		      sum += psi[m*M_G + g] * UF[((tGlo*N_F + fc)*2 + d)*M_s + m];
		      //sum += psi(m,g) * UF(tGlo, fc, d, m);
		    }
		  UintegF[((tGlo*N_F + fc)*2 + d)*M_G + g] = sum;
		  //UintegF(tGlo, g, fc, d) = sum;
		}
	    }
	}
    }
}

//========================================================================================================================================================

extern "C" 
void LUhat_to_UhCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf)
{
  //  printf("Preparing to call Uhat_to_UhCommonHalf\n");
#ifdef USE_CPU
  //makeZero(UhCommonHalf, M_T*M_G*N_F*2);
  //no need to zero here, I do it in the subroutine
#endif
#ifdef USE_GPU
  cudaMemset(UhCommonHalf, (scalar)0.0, M_T*M_G*N_F*2*sizeof(scalar));
#endif
  /*
    printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("t=%d,",t);
  for (int g = 0; g < M_G; g++)
    {
      printf("g=%d,",g);
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
  /*
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock((D+N_N),M_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif
  */
#ifdef USE_CPU
  Uhat_to_UhCommonHalf arch_args (N_E, Ne_AUG, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, PSIxR_elemwise, Uhat, UhCommonHalf);
#endif
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;

  //Configuration 1:
  dim3 dimBlock(M_G , D+N_N , blkE); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  //End Configuration 1.


  /*
  //Configuration 2:
  dim3 dimBlock(M_G , D+N_N , N_E);
  dim3 dimGrid(div+mod,1);
  //End Configuration 2.
  */
  
  Uhat_to_UhCommonHalf_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, PSIxR_elemwise, Uhat, UhCommonHalf);
  cudaError_t err  = cudaThreadSynchronize();
  //  printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
  /*  
printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
  for (int g = 0; g < M_G; g++)
    {
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
}

extern "C" 
void LUhCommonHalf_to_UhCommon(int M_T, int M_G, scalar* UhCommonHalf, scalar* UhCommon)
{
  //  printf("Preparing to call UhCommonHalf_to_UhCommon\n");
#ifdef USE_CPU
   UhCommonHalf_to_UhCommon arch_args(M_T, M_G, UhCommonHalf, UhCommon);
#endif
#ifdef USE_GPU
   int div = M_T/blkT;
   int mod = 0;
   if (M_T%blkT != 0) mod = 1;
   dim3 dimBlock(N_F , M_G , blkT); //thread organization within each block
   dim3 dimGrid(div+mod); //block organization
   UhCommonHalf_to_UhCommon_GPU<<<dimGrid, dimBlock>>> (M_T, M_G, UhCommonHalf, UhCommon);
   
#endif
  /*
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T % blkT != 0) mod = 1;
  dim3 dimBlock(M_G, N_F, blkT);
  dim3 dimGrid(div+mod,1);
#endif
  UhCommonHalf_to_UhCommon arch_args(M_T, M_G, UhCommonHalf, UhCommon);
  */
 }

extern "C"
void LUhat_to_UhCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* Alt_FaceFromElem, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommon)
{
    /*!
    \brief Host C function to launch Uhat_to_UhCommon_v2 kernel.
    In GPU mode, launches N_E/blkE blocks of (D+N_N) x M_G x blkE
    threads. blkE controls the number of interfaces to set on each block
  */
  printf("\n\nPreparing to Call LUhat_to_UhCommon\n\n");
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock((D+N_N),M_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Uhat_to_UhCommon_v2 arch_args (N_E, N_N, M_T, M_G, N_s, RecoPair, Alt_FaceFromElem, PSIxR_elemwise, Uhat, UhCommon);
}

extern "C"
void LUhat_to_GradCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommon)
{
    /*!
    \brief Host C function to launch Uhat_to_gradCommon_v2 kernel.
    In GPU mode, launches N_E/blkE blocks of (D+N_N) x M_G x blkE
    threads. blkE controls the number of interfaces to set on each block
  */
#ifdef USE_CPU
  Uhat_to_GradCommon_v2 arch_args (N_E, N_N, M_T, M_G, N_s, RecoPair, SigFace_from_Uhat_elemwise, Uhat, gradCommon);
#endif
#ifdef USE_GPU
  printf("Go try something else\n");
  /*
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(M_G,(D+N_N),blkE);
  dim3 dimGrid(div+mod);
  Uhat_to_GradCommonHalf_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, RecoPair, SigFace_from_Uhat_elemwise, Uhat, gradCommonHalf);
  */

  /*
  //repartition for GradCommonHalf to GradCommon
  div = M_T/blkT;
  mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dimBlock(N_F, M_G, blkT);
  dimGrid(div+mod);
  GradCommonHalf_to_GradCommon_GPU<<<dimGrid, dimBlock>>>
  */
#endif

  
}

extern "C"
void LUhat_to_GradCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf)
{

#ifdef USE_CPU
  //makeZero(gradCommonHalf, M_T*M_G*N_F*D*2);
  Uhat_to_GradCommonHalf arch_args (N_E, Ne_AUG, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, SigFace_from_Uhat_elemwise, Uhat, gradCommonHalf);
#endif
#ifdef USE_GPU
  cudaMemset(gradCommonHalf, (scalar)0.0, M_T*M_G*N_F*D*2*sizeof(scalar));
#endif
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(M_G , D+N_N , blkE); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  Uhat_to_GradCommonHalf_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, SigFace_from_Uhat_elemwise, Uhat, gradCommonHalf);
  cudaError_t err  = cudaThreadSynchronize();
  // printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
}
 
extern "C"
void LGradCommonHalf_to_GradCommon(int M_T, int M_G, scalar* GradCommonHalf, scalar* GradCommon)
{
#ifdef USE_CPU
  GradCommonHalf_to_GradCommon arch_args (M_T, M_G, GradCommonHalf, GradCommon);
#endif
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(N_F , M_G , blkT); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  GradCommonHalf_to_GradCommon_GPU<<<dimGrid, dimBlock>>> (M_T, M_G, GradCommonHalf, GradCommon);
  cudaError_t err  = cudaThreadSynchronize();
  //  printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
}

extern "C"
void LUhat_to_GradCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf)
{

#ifdef USE_CPU
  //makeZero(gradCommonHalf, M_T*M_G*N_F*D*2);
  Uhat_to_GradCommonHalf_AD arch_args (N_E, Ne_AUG, N_N, M_T, M_G, N_s, Cthresh, elemCmax, sensor, RecoPair, BinarySideAddress, SigFace_from_Uhat_elemwise, Uhat, gradCommonHalf);
#endif
#ifdef USE_GPU
  cudaMemset(gradCommonHalf, (scalar)0.0, M_T*M_G*N_F*D*2*sizeof(scalar));
#endif
#ifdef USE_GPU
  printf("LUhat_to_GradCommonHalf_AD is not gpu-ready\n");
  /*
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(M_G , D+N_N , blkE); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  Uhat_to_GradCommonHalf_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, SigFace_from_Uhat_elemwise, Uhat, gradCommonHalf);
  cudaError_t err  = cudaThreadSynchronize();
  */
  // printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
}
 
extern "C"
  void LGradCommonHalf_to_GradCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* GradCommonHalf, scalar* GradCommon)
{
#ifdef USE_CPU
  GradCommonHalf_to_GradCommon_AD arch_args (M_T, M_G, Cthresh, elemCmax, sensor, BR2_Map, GradCommonHalf, GradCommon);
#endif
#ifdef USE_GPU
  printf("LGradCommonHalf_to_GradCommon_AD is not gpu-ready\n");
  /*
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(N_F , M_G , blkT); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  GradCommonHalf_to_GradCommon_GPU<<<dimGrid, dimBlock>>> (M_T, M_G, GradCommonHalf, GradCommon);
  cudaError_t err  = cudaThreadSynchronize();
  */
  //  printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
}
   
extern "C"
void LUhat_to_UhCommon(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* Uhat, scalar* UhCommon){
  /*!
    \brief Host C function to launch Uhat_to_UhCommon kernel.
    \param[in] N_F number of field variables
    \param[in] N_s solution nodes per element
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[in] BR2_Map tells each interface who the neighboring elements are
    \param[out] UhCommon common solution value on interfaces
    \param[in] PSIxR_Global is the recovery operator
    \param[in] Uhat is global DG DOF storage
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_G x N_F x blkT
    threads. blkT controls the number of interfaces to set on each block
  */
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  Uhat_to_UhCommon arch_args (M_T, M_G, N_s, BR2_Map, PSIxR_Global, Uhat, UhCommon);
}

extern "C"
void LUhat_to_Sigma(int N_E, int N_G, int N_s, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el)
{
  //For this routine: copying Marc's implementation of Lsolver, which includes
  //the inverse mass matrix of each element multiplied by global residual.
  //I've noticed that these L(whatever) functions always use 3 arguments in
  //dimBlock. I'm not sure what this argument does, though I've concluded that
  //the thread arguments in the kernel have to match
  //with what is passed in to dimBlock

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(D,N_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Uhat_to_Sigma arch_args (N_E, N_G, N_s, serial_MSigVol, Uhat, gradU_el);
}

extern "C"
void LUhCommon_to_Sigma(int N_E, int N_G, int N_N, int M_G, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el)
{
  
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(D,N_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  UhCommon_to_Sigma arch_args (N_E, N_G, N_N, M_G, serial_MSigSurf, Alt_FaceFromElem, UhCommon, gradU_el);
}

extern "C"
void LUhat_to_CUGUF_OLD(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* serial_SigFace_from_DOF, scalar* Uhat, scalar* UhCommon, scalar* gradCommon){
  /*!
    \brief Host C function to launch Uhat_to_UhCommon kernel.
    \param[in] N_F number of field variables
    \param[in] N_s solution nodes per element
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[in] BR2_Map tells each interface who the neighboring elements are
    \param[out] UhCommon common solution value on interfaces
    \param[in] PSIxR_Global is the recovery operator
    \param[in] serial_SigFace_from_DOF gets interface gradient from DOF
    \param[in] Uhat is global DG DOF storage
    \param[out] gradCommon is the interface gradient
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_G x N_F x blkT
    threads. blkT controls the number of interfaces to set on each block
  */
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T%blkT != 0) mod = 1;
  dim3 dimBlock(M_G,N_F,blkT);
  dim3 dimGrid(div+mod,1);
#endif

  Uhat_to_CUGUF_OLD arch_args (M_T, M_G, N_s, BR2_Map, PSIxR_Global, serial_SigFace_from_DOF, Uhat, UhCommon, gradCommon);
}

extern "C" void LUhat_to_UicbHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_biased_elemwise, scalar* Uhat, scalar* UicbHalf)
{
  //printf("Preparing to call Uhat_to_UicbHalf\n");
#ifdef USE_CPU
  //makeZero(UhCommonHalf, M_T*M_G*N_F*2);
  //no need to zero here, I do it in the subroutine
#endif
#ifdef USE_GPU
  //Start by setting UicbHalf because I'll be using atomicadd (+=) approach
  cudaMemset(UicbHalf, (scalar)0.0, M_T*M_G*2*N_F*2*sizeof(scalar));
#endif
  /*
    printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("t=%d,",t);
  for (int g = 0; g < M_G; g++)
    {
      printf("g=%d,",g);
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
  /*
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock((D+N_N),M_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif
  */
#ifdef USE_CPU
  Uhat_to_UicbHalf arch_args (N_E, Ne_AUG, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, PSIxR_biased_elemwise, Uhat, UicbHalf);
#endif
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;

  //Configuration 1:
  dim3 dimBlock(2 , D+N_N , blkE); //thread organization within each block
  dim3 dimGrid(div+mod, 1); //block organization
  //End Configuration 1.


  /*
  //Configuration 2:
  dim3 dimBlock(M_G , D+N_N , N_E);
  dim3 dimGrid(div+mod,1);
  //End Configuration 2.
  */
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) 
    printf("Error: %s\n", cudaGetErrorString(err));
  Uhat_to_UicbHalf_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, RecoPair, BinarySideAddress, PSIxR_biased_elemwise, Uhat, UicbHalf);
  cudaDeviceSynchronize();
  /*cudaError_t*/ err = cudaGetLastError();
  if (err != cudaSuccess) 
    printf("Error: %s\n", cudaGetErrorString(err));
  //cudaError_t err  = cudaThreadSynchronize();
  //  printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
  /*  
printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
  for (int g = 0; g < M_G; g++)
    {
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
}

extern "C" void LUicbHalf_to_Uicb(int M_T, int M_G, scalar* UicbHalf, scalar* Uicb)
{
  
#ifdef USE_CPU
   UicbHalf_to_Uicb arch_args(M_T, M_G, UicbHalf, Uicb);
#endif

#ifdef USE_GPU
   cudaError_t err = cudaGetLastError();
   if (err != cudaSuccess) 
     printf("Error: %s\n", cudaGetErrorString(err));
   int div = M_T/blkT;
   int mod = 0;
   if (M_T%blkT != 0) mod = 1;
   dim3 dimBlock(2 , M_G , blkT); //thread organization within each block
   //dim3 dimBlock(1,1,1);
   dim3 dimGrid(div+mod, 1); //block organization
   printf("M_G=%d, blkT=%d, div=%d, mod=%d\n", M_G, blkT, div, mod);
   printf("---About to call UicbHalf_to_Uicb_GPU---\n");
   UicbHalf_to_Uicb_GPU<<<dimGrid, dimBlock>>> (M_T, M_G, UicbHalf, Uicb);
   cudaDeviceSynchronize();
#endif
}

extern "C" void LUhat_to_UicbDirect(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_biased_Global, scalar* Uhat, scalar* Uicb)
{
  /*!
    \brief Host C function to launch Uhat_to_UhCommon kernel.
    \param[in] N_F number of field variables
    \param[in] N_s solution nodes per element
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_T number of interfaces
    \param[in] BR2_Map tells each interface who the neighboring elements are
    \param[out] UhCommon common solution value on interfaces
    \param[in] PSIxR_Global is the recovery operator
    \param[in] Uhat is global DG DOF storage
    \section Description
    In GPU mode, launches M_T/blkT blocks of M_G x N_F x blkT
    threads. blkT controls the number of interfaces to set on each block
  */
#ifdef USE_GPU
  printf("LUhat_to_UicbDirect not ready for GPU\n");
#endif

  Uhat_to_UicbDirect arch_args (M_T, M_G, N_s, BR2_Map, PSIxR_biased_Global, Uhat, Uicb);
}

extern "C" void LCorrectUhCommonBC(int M_T, int M_G, int M_B, int* boundaryMap, scalar* UintegF, scalar* UhCommon)
{
  //Don't mess with UintegF, just use it to adjust UhCommon where necessary
#ifdef USE_GPU
  printf("Don't call LCOrrectUhCommonBC on GPU architechture\n");
#endif
#ifdef USE_CPU
  CorrectUhCommonBC arch_args(M_T, M_G, M_B, boundaryMap, UintegF, UhCommon);
#endif
}

extern "C" void LUhat_to_GradCommonBC(int M_B, int M_G, int N_s, int N_N, int* boundaryMap, int* BR2_Map, scalar* serial_Uhat2GradBC, scalar* Uhat, int* RecoPair, int* Alt_FaceFromElem, scalar* UhCommon, scalar* serial_UhCommon2GradBC, scalar* GradCommon)
{
#ifdef USE_GPU
  printf("LUhat_to_GradCommonBC not ready for gpu\n");
#endif
  Uhat_to_GradCommonBC arch_args(M_B, M_G, N_s, N_N, boundaryMap, BR2_Map, serial_Uhat2GradBC, Uhat, RecoPair, Alt_FaceFromElem, UhCommon, serial_UhCommon2GradBC, GradCommon);
}

//Some stuff for artificial dissipation:
//Some stuff for artificial dissipation:
extern "C" void LUhat_to_UhCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf)
{
  //  printf("Preparing to call Uhat_to_UhCommonHalf\n");
#ifdef USE_CPU
  //makeZero(UhCommonHalf, M_T*M_G*N_F*2);
  //no need to zero here, I do it in the subroutine
#endif
#ifdef USE_GPU
  cudaMemset(UhCommonHalf, (scalar)0.0, M_T*M_G*N_F*2*sizeof(scalar));
#endif
  /*
    printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
      printf("t=%d,",t);
  for (int g = 0; g < M_G; g++)
    {
      printf("g=%d,",g);
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
  /*
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock((D+N_N),M_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif
  */
#ifdef USE_CPU
  Uhat_to_UhCommonHalf_AD arch_args (N_E, Ne_AUG, N_N, M_T, M_G, N_s, Cthresh, elemCmax, sensor, RecoPair, BinarySideAddress, PSIxR_elemwise, Uhat, UhCommonHalf);
#endif
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;

  //Configuration 1:
  dim3 dimBlock(M_G , D+N_N , blkE); //thread organization within each block
  dim3 dimGrid(div+mod); //block organization
  //End Configuration 1.


  /*
  //Configuration 2:
  dim3 dimBlock(M_G , D+N_N , N_E);
  dim3 dimGrid(div+mod,1);
  //End Configuration 2.
  */
  
  Uhat_to_UhCommonHalf_AD_GPU<<<dimGrid, dimBlock>>> (N_E, N_N, M_T, M_G, N_s, Cthresh, elemCmax, sensor, RecoPair, BinarySideAddress, PSIxR_elemwise, Uhat, UhCommonHalf);
  cudaError_t err  = cudaThreadSynchronize();
  //  printf("Run kernel: %s\n", cudaGetErrorString(err));
#endif
  /*  
printf("UhCommonHalf Distribution:\n");
  for (int t = 0; t < M_T; t++)
    {
  for (int g = 0; g < M_G; g++)
    {
  for (int fc = 0; fc < N_F; fc++)
    {
  printf("Interface node+f (%d,%d,%d): contributions are %f and %f\n", t, g, fc, UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 0], UhCommonHalf[t*M_G*N_F*2 + g*N_F*2 + fc*2 + 1]);
}
}
}
  */
}

extern "C" void LUhCommonHalf_to_UhCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* UhCommonHalf, scalar* UhCommon)
{
  //  printf("Preparing to call UhCommonHalf_to_UhCommon\n");
#ifdef USE_CPU
  UhCommonHalf_to_UhCommon_AD arch_args(M_T, M_G, Cthresh, elemCmax, sensor, BR2_Map, UhCommonHalf, UhCommon);
#endif
#ifdef USE_GPU
   int div = M_T/blkT;
   int mod = 0;
   if (M_T%blkT != 0) mod = 1;
   dim3 dimBlock(N_F , M_G , blkT); //thread organization within each block
   dim3 dimGrid(div+mod); //block organization
   UhCommonHalf_to_UhCommon_AD_GPU<<<dimGrid, dimBlock>>> (M_T, M_G, Cthresh, elemCmax, sensor, BR2_Map, UhCommonHalf, UhCommon);
   
#endif
  /*
#ifdef USE_GPU
  int div = M_T/blkT;
  int mod = 0;
  if (M_T % blkT != 0) mod = 1;
  dim3 dimBlock(M_G, N_F, blkT);
  dim3 dimGrid(div+mod,1);
#endif
  UhCommonHalf_to_UhCommon arch_args(M_T, M_G, UhCommonHalf, UhCommon);
  */
 }

extern "C" void LUhat_to_Sigma_AD(int N_E, int N_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el)
{
  //For this routine: copying Marc's implementation of Lsolver, which includes
  //the inverse mass matrix of each element multiplied by global residual.
  //I've noticed that these L(whatever) functions always use 3 arguments in
  //dimBlock. I'm not sure what this argument does, though I've concluded that
  //the thread arguments in the kernel have to match
  //with what is passed in to dimBlock

#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(D,N_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  Uhat_to_Sigma_AD arch_args (N_E, N_G, N_s, Cthresh, elemCmax, sensor, serial_MSigVol, Uhat, gradU_el);
}

extern "C" void LUhCommon_to_Sigma_AD(int N_E, int N_G, int N_N, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el)
{
  
#ifdef USE_GPU
  int div = N_E/blkE;
  int mod = 0;
  if (N_E%blkE != 0) mod = 1;
  dim3 dimBlock(D,N_G,blkE);
  dim3 dimGrid(div+mod,1);
#endif

  UhCommon_to_Sigma_AD arch_args (N_E, N_G, N_N, M_G, Cthresh, elemCmax, sensor, BR2_Map, serial_MSigSurf, Alt_FaceFromElem, UhCommon, gradU_el);
}

//Correction to boundary solution distribution when using ICBN algorithm
extern "C" void LBuildBC_UintegF(int M_B, int M_G, int M_s, int* boundaryMap, scalar* psi, scalar* UF, scalar* UintegF)
{
#ifdef USE_GPU
  printf("LBuildBC_UintegF is not GPU-compatible yet\n");
#endif
#ifdef USE_CPU
  BuildBC_UintegF arch_args (M_B, M_G, M_s, boundaryMap, psi, UF, UintegF);
#endif
}

//My combined UhComman and GradCommon operation: Abandon element-wise summation
extern "C" void LUhat_to_CUGUF(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* Sigface_from_Uhat_ROG, scalar* Uhat, scalar* UhCommon, scalar* GradCommon)
{
#ifdef USE_GPU
  printf("DO NOT USE GPU WITH LUhat_to_CUGUF\n");
#endif
#ifdef USE_CPU
  Uhat_to_CUGUF arch_args (M_T, M_G, N_s, BR2_Map, PSIxR_Global, Sigface_from_Uhat_ROG, Uhat, UhCommon, GradCommon);
#endif
}
