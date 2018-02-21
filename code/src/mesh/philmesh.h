/*!
  \file philMesh.h  
  \brief Class deals with the mesh.
  \copyright Copyright (C) 2012-2020, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Philip E. Johnson <phedjohn@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  Subroutine(s) to supplement mesh usage. Created to deal with HiOCFD5 vortex problem
*/
#ifndef _PHIL_MESH_H_
#define _PHIL_MESH_H_
#include <vector>
#include "fullMatrix.h"
#include "scalar_def.h"
#include <string>
#include <map>
#include "macros.h"
#include <math.h>
#include "simpleMesh.h"
#include "polynomialBasis.h"

fullMatrix<double> PeriFixNodes(int PeriCo2D, const fullMatrix<double> &nodesCONST);

fullMatrix<double> PeriFixNodes_HiOMesh(int PeriCo2D, fullMatrix<double> nodes, /*const std::vector<simpleElement> &elements*/ std::vector<std::vector<int> > elemNodes, int elem_type, int order, int N_N);

#endif
