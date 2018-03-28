/*!
  \file recovery_tools.h
  \Header file for Recovery operations, for both full and biased recovery
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup gradient
*/
#ifndef RECOVERY_TOOLS_H
#define RECOVERY_TOOLS_H

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
#include "timers.h"

void COPY_get_element_types(const int order, int &msh_hex, int &msh_qua, int &msh_tri, int &msh_lin);

scalar Leg_L2(int index, int p, scalar x);

void RecoveryIndices_2DSimplex(int N_s, int p, int N_basis, int* subIndex_relay);

void PsiAtPoints(int N_s, int p, int N_basis, int N_N, int Npoints, scalar* x, scalar* PsiOut);

void PeriodicityFix(int N_s, int GQresFace, scalar* xNodes_A, scalar* xNodes_B, scalar* xFace, int* FlagPeri_out, int* sign_shift_relay, scalar* Ldomain_relay);

scalar distance(scalar* X1, scalar* X2);

void get_CrossProduct(scalar* U, scalar* V, scalar* perp);

void TransformPhysical2Recovery(int N_s, int GQresA, int GQresFace, scalar* xNodes_A, scalar* xNodes_B, scalar* xGQ_A, scalar* xGQ_B, scalar* xFace, scalar* FaceNormal, scalar* rA_out, scalar* rB_out, scalar* rFace_out);

void PSIxR_shell(deck inputs, TIMERS &timers, int N_N, int N_s, int p, /*simpleMesh &m,*/ int M_G, int Style, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, int BoundaryTag, fullMatrix<scalar> &psi, int* Better_InvMap, int t, int M_s, int omA, int omB,scalar* PSIxR, scalar* PSIxR_A, scalar* PSIxR_B);

void GetPSIxR_oneface(int N_N, TIMERS &timers, int N_s, int p, int N_basis, int M_G, int GradMethod, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, int BoundaryTag, fullMatrix<scalar> &psi, int* Better_InvMap, int t, int M_s, int omA, int omB, scalar* PSIxR, scalar* PSIxR_A, scalar* PSIxR_B);




scalar Dist2Edge_2D(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2);

scalar Dist2Plane_3D(scalar x1, scalar y1, scalar z1, 
		     scalar xA, scalar yA, scalar zA, 
		     scalar xB, scalar yB, scalar zB, 
		     scalar xC, scalar yC, scalar zC);

void GrabPlaneVertices_3D(int GQresFace, scalar* xGeo_A, scalar* xFace, scalar* XV1, scalar* XV2, scalar* XV3, scalar* XV4);

void Reco3DVector_From_Vertices(scalar* X1, scalar* X2, scalar* X3, scalar* X4, scalar* vr, scalar* vs, scalar* vt);

void NodalEqCord(deck inputs, int N_s, int p, int M_G, int N_EP, int N_N, const polynomialBasis* RefBasis, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, scalar* xirgA, scalar* xirgB, scalar* RrgA, scalar* RrgB);

void ICB_Cords_from_Physical(int N_s, int N_N, int GQresA, int GQresFace, int N_EP, scalar* xFace, scalar* xGeo_A, scalar* xGeo_B, scalar* xGQ_A, scalar* xGQ_B, scalar* XrgA, scalar* XrgB, /*begin output*/ scalar* rA_omA, scalar* rB_omB, scalar* rFace_omA, scalar* rFace_omB, scalar* r_rgB_omA, scalar* r_rgA_omB);

void GetPSIxR_Biased_oneface_Nov2017(deck inputs, int N_N, TIMERS &timers, int N_s, int p, /*simpleMesh &m,*/ int N_basis, int M_G, const polynomialBasis *basis, int N_superG, fullMatrix<scalar> &phiRef, fullMatrix<scalar> &dphi, fullMatrix<double> &points, fullMatrix<double> &weight, scalar* xGeo_A, scalar* xGeo_B, scalar* xFace, scalar* FaceNormal, scalar* PSIxR_biA, scalar* PSIxR_biB, scalar* PSIxR_biA_A, scalar* PSIxR_biA_B, scalar* PSIxR_biB_A, scalar* PSIxR_biB_B);

#endif
