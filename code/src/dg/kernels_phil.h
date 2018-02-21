/*!
  \file kernels.h
  \brief Functions to launch some kernels specific to DG
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef KERNELS_PHIL_H
#define KERNELS_PHIL_H
#include "scalar_def.h"
#include "macros.h"

// Here I define the gpu kernels I will be using
extern "C" void LUhat_to_UhCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* Alt_FaceFromElem, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommon);

extern "C" void LUhat_to_GradCommon_v2(int N_E, int N_N, int M_T, int M_G, int N_s, int* RecoPair, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommon);

extern "C" void LUhat_to_UhCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf);

extern "C" void LUhCommonHalf_to_UhCommon(int M_T, int M_G, scalar* UhCommonHalf, scalar* UhCommon);

extern "C" void LUhat_to_GradCommonHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf);

extern "C" void LGradCommonHalf_to_GradCommon(int M_T, int M_G, scalar* GradCommonHalf, scalar* GradCommon);

extern "C" void LUhat_to_GradCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* SigFace_from_Uhat_elemwise, scalar* Uhat, scalar* gradCommonHalf);


extern "C" void LGradCommonHalf_to_GradCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* GradCommonHalf, scalar* GradCommon);

extern "C" void LUhat_to_UhCommon(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* Uhat, scalar* UhCommon);

extern "C" void LUhat_to_Sigma(int N_E, int N_G, int N_s, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el);

extern "C" void LUhCommon_to_Sigma(int N_E, int N_G, int N_N, int M_G, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el);

extern "C" void LUhat_to_CUGUF_OLD(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* serial_SigFace_from_DOF, scalar* Uhat, scalar* UhCommon, scalar* gradCommon);

extern "C" void LUhat_to_UicbHalf(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_biased_elemwise, scalar* Uhat, scalar* UicbHalf);

extern "C" void LUicbHalf_to_Uicb(int M_T, int M_G, scalar* UicbHalf, scalar* Uicb);

extern "C" void LUhat_to_UicbDirect(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_biased_Global, scalar* Uhat, scalar* Uicb);

extern "C" void LCorrectUhCommonBC(int M_T, int M_G, int M_B, int* boundaryMap, scalar* UintegF, scalar* UhCommon);

extern "C"void LUhat_to_GradCommonBC(int M_B, int M_G, int N_s, int N_N, int* boundaryMap, int* BR2_Map, scalar* serial_Uhat2GradBC, scalar* Uhat, int* RecoPair, int* Alt_FaceFromElem, scalar* UhCommon, scalar* serial_UhCommon2GradBC, scalar* GradCommon);

//A bunch of stuff for AD: it's like regular viscous routines with some extra arguments
extern "C" void LUhat_to_UhCommonHalf_AD(int N_E, int Ne_AUG, int N_N, int M_T, int M_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, int* RecoPair, int* BinarySideAddress, scalar* PSIxR_elemwise, scalar* Uhat, scalar* UhCommonHalf);

extern "C" void LUhCommonHalf_to_UhCommon_AD(int M_T, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* UhCommonHalf, scalar* UhCommon);

extern "C" void LUhat_to_Sigma_AD(int N_E, int N_G, int N_s, scalar Cthresh, scalar* elemCmax, int* sensor, scalar* serial_MSigVol, scalar* Uhat, scalar* gradU_el);

extern "C" void LUhCommon_to_Sigma_AD(int N_E, int N_G, int N_N, int M_G, scalar Cthresh, scalar* elemCmax, int* sensor, int* BR2_Map, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* UhCommon, scalar* gradU_el);

extern "C" void LBuildBC_UintegF(int M_B, int M_G, int M_s, int* boundaryMap, scalar* psi, scalar* UF, scalar* UintegF);

extern "C" void LUhat_to_CUGUF(int M_T, int M_G, int N_s, int* BR2_Map, scalar* PSIxR_Global, scalar* Sigface_from_Uhat_ROG, scalar* Uhat, scalar* UhCommon, scalar* GradCommon);

#endif

