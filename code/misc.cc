#include <misc.h>

void blasScopy(int N, float* x, int INCX, float* y, int INCY){
  F77NAME(scopy)(&N, x, &INCX, y, &INCY);
}
void blasSaxpy(int M, float alpha, float* x, int INCX, float* y, int INCY){
  F77NAME(saxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasSgemm(char or1, char or2, int M , int N, int K, float alpha, float* A, int LDA, float* B, int LDB, float beta, float* C, int LDC){
  F77NAME(sgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
int blasIsamax(int N, float* x, int INCX){
  F77NAME(isamax)(&N, x, &INCX);
}
void blasDcopy(int N, double* x, int INCX, double* y, int INCY){
  F77NAME(dcopy)(&N, x, &INCX, y, &INCY);
}
void blasDaxpy(int M, double alpha, double* x, int INCX, double* y, int INCY){
  F77NAME(daxpy)(&M, &alpha, x ,&INCX, y, &INCY);
}
void blasDgemm(char or1, char or2, int M , int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double beta, double* C, int LDC){
  F77NAME(dgemm)(&or1, &or2, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
int blasIdamax(int N, double* x, int INCX){
  F77NAME(idamax)(&N, x, &INCX);
}

void makeZero(scalar* A, int size){
  for(int k=0; k < size; k++) A[k] = 0.0;
}

int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void readTable(const char *fileName, fullMatrix<scalar> &XWUP, scalar &gamma, scalar &alpha){
  /* Read the data from a table in a txt file

     It is used by init_cond.cc
     
   */

  // Open the file
  std::ifstream table;
  table.open(fileName,std::ifstream::in);
  if(table.is_open()==0){ printf("No file named %s. Exiting.\n",fileName);}

  int nrows;

  // First line contains gamma, alpha, number of rows
  table >> gamma >> alpha >> nrows;

  // Read the rest of the table, populate XWUP
  XWUP.resize(nrows,4);
  for (int k = 0; k < nrows; k++) table >> XWUP(k,0) >> XWUP(k,1) >> XWUP(k,2) >> XWUP(k,3);
  
  table.close();

}

scalar interpolate(scalar x, std::vector<std::pair<scalar, scalar> > table, scalar BCL, scalar BCR) {
  /*
    From http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation

    Assumes that "table" is sorted by .first

    Linear interpolation of the data in table. If x is out of bounds,
    replace with boundary condition of the left (BCL) or right (BCR).

    It is used by init_cond.cc
    
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
