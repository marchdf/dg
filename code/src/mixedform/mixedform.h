/*!
  \file mixedform.h
  \Header file for dealing with problems with 2nd order derivatives
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/
#ifndef MIXEDFORM_H
#define MIXEDFORM_H

#include "scalar_def.h"
#include <math.h>
#include "macros.h"
#include <stdio.h>
#include "fullMatrix.h"
#include "polynomialBasis.h"
#include "polynomialsJacobi.h"
#include "Gauss.h"
#include "GmshDefines.h"
#include "simpleMesh.h"
#include "scalar_def.h"
#include "dg_functions.h"
#include "deck.h"

void COPYCOPY_get_element_types(const int order, int &msh_hex, int &msh_qua, int &msh_tri, int &msh_lin);

void SurfForSigma(deck inputs, int N_s, int p, int om, int M_G, int N_N, int M_s, int a, int* invmap, int* order_t, fullMatrix<double> pointsF, fullMatrix<double> weightF, scalar* JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> OutwardNormal, fullMatrix<scalar> psi_w, /*int* indices_support, scalar* Sa_return*/ fullMatrix<scalar> &Sa_return);

//void SurfForSigma(deck inputs, int N_s, int p, int om, int M_G, int N_N, int M_s, int a, int* invmap, fullMatrix<double> pointsF, fullMatrix<double> weightF, scalar* JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> OutwardNormal, fullMatrix<scalar> psi_w, /*int* indices_support, scalar* Sa_return*/ fullMatrix<scalar> &Sa_return);

//void SurfForSigma(deck inputs, int N_s, int p, int om, int M_G, int N_N, int M_s, int a, int* rog_invmap_local, /*int* invmap,*/ fullMatrix<double> pointsF, fullMatrix<double> weightF, scalar* JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> OutwardNormal, fullMatrix<scalar> psi_w, /*int* indices_support, scalar* Sa_return*/ fullMatrix<scalar> &Sa_return);

void VolForSigma(int om, deck inputs, int QuadRule, int N_G, int M_G, int N_s, int p, int a, scalar* invJac, scalar detJ, scalar* detJ_full,  /*scalar* Va_return*/ fullMatrix<scalar> &Va_return);


void BuildSigMatrices(deck inputs, int QuadRule, int N_E, int N_s, int p, int N_G, int M_G, int N_N, int M_s, int* invmap, fullMatrix<double> &pointsF, fullMatrix<double> &weightF, fullMatrix<scalar> &JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> normals, scalar* h_invJac, fullMatrix<scalar> &J, int* Alt_FaceFromElem, int* BR2_map, fullMatrix<scalar> phi, scalar* Minv_Master, fullMatrix<scalar> psi_w, scalar* detJ_full, scalar* serial_SigSurf, scalar* serial_SigVol);

void BuildAuxHatMatrices(deck inputs, int QuadRule, int N_E, int N_s, int p, int N_G, int M_G, int N_N, int M_s, int* invmap, fullMatrix<double> &pointsF, fullMatrix<double> &weightF, fullMatrix<scalar> &JF, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> normals, scalar* h_invJac, fullMatrix<scalar> &J, int* Alt_FaceFromElem, int* BR2_map, fullMatrix<scalar> phi, scalar* Minv_Master, fullMatrix<scalar> psi_w, scalar* detJ_full, scalar* serial_SigSurf, scalar* serial_SigVol);

void SigFaceMatrices(int M_T, int M_G, int N_s, int N_G, int N_N, int M_s, scalar CHI, fullMatrix<scalar> &phi, fullMatrix<scalar> &dphi, fullMatrix<double> &weight, fullMatrix<scalar> &J, fullMatrix<scalar> &psi, fullMatrix<scalar> &psi_w, fullMatrix<scalar> normals, fullMatrix<scalar> &JF, int* Alt_FaceFromElem, int* BR2_Map, int* invmap, int* Better_InvMap, scalar* invJac, scalar* Minv_global, scalar* PSIxR_Global, scalar* detJ_full, scalar* serial_SigFace_from_DOF);

#endif
