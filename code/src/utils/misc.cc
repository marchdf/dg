/*!
  \file misc.cc
  \brief Miscellaneous functions used by host functions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "misc.h"

void blasScopy(int N, float* x, int INCX, float* y, int INCY){
  /*!\brief Call float BLAS copy*/
  F77NAME(scopy)(&N, x, &INCX, y, &INCY);
}
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY){
  /*!\brief Call float BLAS axpy*/
  F77NAME(saxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC){
  /*!\brief Call float BLAS gemm*/
  F77NAME(sgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
int blasIsamax(int N, float* x, int INCX){
  /*!\brief Call float BLAS amax*/
  F77NAME(isamax)(&N, x, &INCX);
}
void blasDcopy(int N, double* x, int INCX, double* y, int INCY){
  /*!\brief Call double BLAS copy*/
  F77NAME(dcopy)(&N, x, &INCX, y, &INCY);
}
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY){
  /*!\brief Call double BLAS axpy*/
  F77NAME(daxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC){
  /*!\brief Call double BLAS gemm*/
  F77NAME(dgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
int blasIdamax(int N, double* x, int INCX){
  /*!\brief Call int BLAS amax*/
  F77NAME(idamax)(&N, x, &INCX);
}

void makeZero(scalar* A, int size){
  /*!\brief Set a vector entries to 0*/
  for(int k=0; k < size; k++) A[k] = 0.0;
}

//template<typename T>
//void hostDeepCopyArray(int* src, int* dest, int size){for(int k=0; k < size; k++) dest[k] = src[k];}

int factorial(int n){
  /*! \brief Recursive factorial function*/
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void readTable(const char *fileName, fullMatrix<scalar> &XWUP, scalar &gamma, scalar &alpha, scalar &Q){
  /*! \brief Read the data from a table in a txt file (used by init_cond.cc) */

  // Open the file
  std::ifstream table;
  table.open(fileName,std::ifstream::in);
  if(table.is_open()==0){ printf("No file named %s. Exiting.\n",fileName); exit(1);}

  int nrows;

  // First line contains gamma, alpha, number of rows
  table >> gamma >> alpha >> Q >> nrows;

  // Read the rest of the table, populate XWUP
  XWUP.resize(nrows,4);
  for (int k = 0; k < nrows; k++) table >> XWUP(k,0) >> XWUP(k,1) >> XWUP(k,2) >> XWUP(k,3);
  
  table.close();

}

scalar interpolate(scalar x, std::vector<std::pair<scalar, scalar> > table, scalar BCL, scalar BCR) {
  /*!
    \brief Linear interpolation of the data in table.
    \section Description
    From http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation

    Assumes that "table" is sorted by .first

    Linear interpolation of the data in table. If x is out of bounds,
    replace with boundary condition of the left (BCL) or right (BCR).

    It is used by init_cond.cc.    
  */

  // Make sure the table is ordered
  std::sort(table.begin(), table.end());

  // Check if x is out of bound
  if (x > table.back().first) return BCR;
  if (x < table[0].first) return BCL;

  std::vector<std::pair<scalar, scalar> >::iterator it, it2;
  it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, (scalar)-INFINITY));   // INFINITY is defined in math.h in the glibc implementation
  if (it == table.begin()) return it->second;   // Corner case
  it2 = it;
  --it2;
  return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}


void get_node_rgb_map(const char *fileName,std::map<int,std::vector<int> > &node_rgb_map){
  /*!
    \brief Read the data from the jpeg_mesh.py script linking the node id to an rgb value of the image
    \param[in] filename name of file containing the data
    \param[out] node_rgb_map map between the node ID and the RGB value at that node
  */

  // Open the file
  std::ifstream ifile;
  ifile.open(fileName,std::ifstream::in);
  if(ifile.is_open()==0){ printf("No file named %s. Exiting.\n",fileName); exit(1);}

  // First line contains the number of rows
  int nrows;
  ifile >> nrows;

  // loop on all the lines
  int node_id;
  int r,g,b=0;
  bool found = false;
  for (int k = 0; k < nrows; k++){

    // read the data
    ifile >> node_id >> r >> g >> b;
    std::vector<int> rgb(3);
    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;
    node_rgb_map[node_id] = rgb;
  }

}
