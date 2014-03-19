/*!
  \file Gauss.cc
  \brief Generate gaussian quadratures
  \copyright Gmsh - Copyright (C) 1997-2010
  \authors C. Geuzaine, J.-F. Remacle
  
  See the LICENSE.txt file for license information. Please report all
  bugs and problems to <gmsh@geuz.org>.
*/
#include "Gauss.h"
static void pts2fullMatrix(int npts, IntPt *pts, fullMatrix<double> &pMatrix, fullMatrix<double> &wMatrix) {
  pMatrix.resize(npts,3);
  wMatrix.resize(npts,1);
  for (int i=0;i<npts;i++) {
    pMatrix(i,0) = pts[i].pt[0];
    pMatrix(i,1) = pts[i].pt[1];
    pMatrix(i,2) = pts[i].pt[2];
    wMatrix(i,0) = pts[i].weight;
  }
}

void gaussIntegration::getTriangle(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQTPts(order),getGQTPts(order),pts,weights);
}
void gaussIntegration::getLine(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQLPts(order),getGQLPts(order),pts,weights);
}
void gaussIntegration::getQuad(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQQPts(order),getGQQPts(order),pts,weights);
}
void gaussIntegration::getTetrahedron(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQTetPts(order),getGQTetPts(order),pts,weights);
}
void gaussIntegration::getHexahedron(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQHPts(order),getGQHPts(order),pts,weights);
}
void gaussIntegration::getPrism(int order, fullMatrix<double> &pts, fullMatrix<double> &weights) {
  pts2fullMatrix(getNGQPriPts(order),getGQPriPts(order),pts,weights);
}
