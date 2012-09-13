#ifndef DG_FUNCTIONS_H
#define DG_FUNCTIONS_H
#include "fullMatrix.h"
#include "simpleMesh.h"
#include "polynomialBasis.h"
#include "quadratures/Gauss.h"
#include <scalar_def.h>
#include <map>
#include <vector>

void dg_jac_elements_fast(const int N_G, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &J);
void dg_jacobians_elements(const int N_G, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &dphi, fullMatrix<scalar> &Jac, fullMatrix<scalar> &invJac, fullMatrix<scalar> &J, fullMatrix<scalar> &invJ);
void dg_jacobians_face(const int M_T, const int D, fullMatrix<scalar> &XYZNodesF, fullMatrix<scalar> &dpsi, fullMatrix<scalar> &JacF, fullMatrix<scalar> &JF, fullMatrix<scalar> &invJF);
void dg_inverse_mass_matrix(const int order, const int elem_type, const std::string getElemType, const int N_s, const int N_E, const int D, fullMatrix<scalar> &XYZNodes, scalar* Minv);
void dg_mappings(const int M_s, const int M_T, const int N_F, const int N_s, const int N_E, const std::vector<simpleInterface> &interfaces, std::map<int,int> &ElementMap, const std::vector<std::vector<int> > &closures, int* map, int* invmap);

#endif
