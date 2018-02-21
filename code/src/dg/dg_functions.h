/*!
  \file dg_functions.h
  \brief Functions for the DG method
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef DG_FUNCTIONS_H
#define DG_FUNCTIONS_H
#include "fullMatrix.h"
#include "simpleMesh.h"
#include "polynomialBasis.h"
#include "Gauss.h"
#include "scalar_def.h"
#include <map>
#include <vector>

void dg_jac_elements_fast(const int N_G, const int N_E, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &J);
void dg_detJ_OneElement(const int N_G, scalar* xGeo, fullMatrix<scalar> &dphi, scalar* detJ);
void dg_jacobians_elements(const int N_G, const int N_E, const int Ne_AUG, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &XYZNodes_extended, fullMatrix<scalar> &dphi, fullMatrix<scalar> &Jac, fullMatrix<scalar> &invJac, fullMatrix<scalar> &J, fullMatrix<scalar> &invJ);
void dg_jacobians_face(const int M_T, fullMatrix<scalar> &XYZNodesF, fullMatrix<scalar> &dpsi, const fullMatrix<scalar> &normals, fullMatrix<scalar> &JacF, fullMatrix<scalar> &JF, fullMatrix<scalar> &invJF);
int LocalSeekGhost(const std::vector<int> &ghostElement_Key, const std::vector<int> &ghostElement_Val, const std::vector<int> &ghostElement_HNe, int GmshId_on, int GmshId_off);
void dg_inverse_mass_matrix(const int order, const int elem_type, const std::string getElemType, const int N_s, const int N_E, const int Ne_AUG, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &XYZNodes_extended, int QuadRule, scalar* detJ_full, scalar* Minv);
void dg_mappings(const int myid, const int M_s, const int M_T, const int N_s, const int N_E, const int N_N, const std::vector<simpleInterface> &interfaces, const std::map<int,int> &ElementMap, const std::vector<int> &ghostElement_Key, const std::vector<int> &ghostElement_Val, const std::vector<int> &ghostElement_HNe, const std::vector<std::vector<int> > &closures, int* map, int* invmap, int* Alt_FaceFromElem_serial);
void dg_mappings_PERIFIX(const int myid, const int M_s, const int M_T, const int N_s, const int N_E, const int N_N, const std::vector<simpleInterface> &interfaces, const std::map<int,int> &ElementMap, const std::map<int,int> &ghostElementMap, const std::vector<std::vector<int> > &closures, int* map, int* invmap, int* Alt_FaceFromElem_serial);

#endif
