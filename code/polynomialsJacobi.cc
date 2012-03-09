#include <polynomialsJacobi.h>


void JacobiP(const fullMatrix<scalar> x, const int alpha, const int beta, const int N, fullMatrix<scalar> &P){

  // Evaluate Jacobi polynomials of type (alpha, beta) > -1
  // at points x, of order N
  // Normalized so that they are orthornormal to each other
  // Inspired by the book: Nodal DG... p446
  int Nx = x.size1();
  fullMatrix<scalar> PL(N+1, Nx);

  // Initial values for p0 and p1
  scalar gamma0 = pow(2,alpha+beta+1)/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
  for(int j=0;j<Nx;j++) PL(0,j) = 1.0/sqrt(gamma0);
  if (N==0){
    for(int j=0;j<Nx;j++) P(j,0) = PL(0,j);
    return;
  }
  scalar gamma1 = (scalar)(alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  for(int j=0;j<Nx;j++) PL(1,j) = ((alpha+beta+2)*x(j,0)*0.5 + (alpha-beta)*0.5)/sqrt(gamma1);
  if (N==1){
    for(int j=0;j<Nx;j++) P(j,0) = PL(1,j);
    return;
  }
  
  // Repeat value in recurrence.
  scalar aold = 2.0/(2+alpha+beta)*sqrt((scalar)(alpha+1)*(beta+1)/(alpha+beta+3));
  // Forward recurrence using the symmetry of the recurrence.
  for(int i=1;i<N;i++){      
    scalar h1 = 2*i+alpha+beta;
    scalar anew = 2.0/(h1+2)*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    scalar bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    for(int j=0;j<Nx;j++) PL(i+1,j) = 1/anew*(-aold*PL(i-1,j) + (x(j,0)-bnew)*PL(i,j));
    aold = anew;
  }

  for(int j=0;j<Nx;j++) P(j,0) = PL(N,j);
  return;
  
}

void JacobiGQ(const int alpha, const int beta, const int N, fullMatrix<scalar> &x, fullMatrix<scalar> &w){

  //Purpose: Compute the N’th order Gauss quadrature points, x,
  //         and weights, w, associated with the Jacobi
  //         polynomial, of type (alpha,beta) > -1 ( <> -0.5).

  if (N==0){
    x(0,0) = (scalar)(alpha-beta)/(alpha+beta+2);
    w(0,0) = 2;
    return;
  }

  // Form symmetric matrix from recurrence.
  fullMatrix<scalar> J(N+1,N+1); J.scale(0.0);
  fullMatrix<scalar> h1(1,N+1);
  for(int j=0;j<N+1;j++) h1(0,j) = 2*j+alpha+beta;

  for(int j=0;j<N+1;j++) J(j,j) = -0.5*(alpha*alpha-beta*beta)/(h1(0,j)+2)/h1(0,j);
  for(int j=0;j<N;j++)   J(j,j+1) = 2.0/(h1(0,j)+2)*sqrt((j+1)*((j+1)+alpha+beta)*((j+1)+alpha)*((j+1)+beta)/(h1(0,j)+1)/(h1(0,j)+3));

  if (alpha+beta<10*2.2204e-16) J(0,0)=0.0;

  J.add(J.transpose());

  // Compute quadrature by eigenvalue solve
  fullVector<double> eigenValReal(N+1);
  fullVector<double> eigenValImag(N+1);
  fullMatrix<scalar> leftEigenVect(N+1,N+1);
  fullMatrix<scalar> rightEigenVect(N+1,N+1);
  bool sortRealPart=true;
  J.eig(eigenValReal, eigenValImag, leftEigenVect, rightEigenVect, sortRealPart);

  for(int j=0;j<N+1;j++){
    x(j,0) = (scalar)eigenValReal(j);
    w(j,0) = (scalar)rightEigenVect(0,j)*rightEigenVect(0,j)*pow(2,alpha+beta+1)/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
  }
  return;
}


void JacobiGL(const int alpha, const int beta, const int N, fullMatrix<scalar> &x){

  // Purpose: Compute the N’th order Gauss Lobatto quadrature
  //          points, x, associated with the Jacobi polynomial,
  //          of type (alpha,beta) > -1 ( <> -0.5).

  if (N==1){
    x(0,0) =-1.0;
    x(1,0) = 1.0;
    return;
  }
  fullMatrix<scalar> xint(N-2+1,1);
  fullMatrix<scalar> w(N-2+1,1);
  JacobiGQ(alpha+1,beta+1,N-2,xint,w);

  x(0,0) =-1;
  x(N,0) = 1;
  for(int j=1;j<N;j++) x(j,0)=xint(j-1,0);
  return;
}
