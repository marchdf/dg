#ifndef LIMITING_H
#define LIMITING_H

#include <stdio.h>
#include <scalar_def.h>
#include <macros.h>
#include <fullMatrix.h>
#include <misc.h>
#include <cpu_kernels.h>

class Limiting
{
 private:
  int     _order;  // RK order (only implemented RK4)
  scalar* _Lag2Mono;
  scalar* _Mono2Lag;
  scalar* _V1D;
  scalar* _weight; // integration weights
  scalar* _A;      // monomial solution
  scalar* _Alim;   // holds the limited monomial solution
  bool    _limit;  // limit?
  int     _method; // HR=1
  int     _N_s;
  int     _N_E;
  int     _N_F;
  int     _N_G;
  int     _boundaryMap;
  
 public:
  // constructor
 Limiting(bool limit) : _limit(limit){}
  
 Limiting(bool limit, int method, int N_s, int N_E, int N_F, int N_G, int boundaryMap,
	  fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, fullMatrix<scalar> &V1D, fullMatrix<scalar> &weight)
   : _limit(limit), _method(method), _N_s(N_s), _N_E(N_E), _N_F(N_F), _N_G(N_G), _boundaryMap(boundaryMap){
    switch (_method){
    case 1:
      _limit    = true;
#ifdef USE_CPU
      _Lag2Mono = new scalar[_N_s*_N_s];     copyMatrixToPointer(Lag2Mono,_Lag2Mono);
      _Mono2Lag = new scalar[_N_s*_N_s];     copyMatrixToPointer(Mono2Lag,_Mono2Lag);
      _V1D      = new scalar[_N_G*_N_s];     copyMatrixToPointer(V1D,_V1D);
      _weight   = new scalar[_N_G];          copyMatrixToPointer(weight,_weight);
      _A        = new scalar[_N_s*_N_E*N_F];
      _Alim     = new scalar[_N_s*_N_E*N_F];
#elif USE_GPU
      
#endif
      break;
    default:
      printf("Bad method input. Defaulting to no limiting.\n");
      _limit = false;
    }
  }
  
  // destructor
  ~Limiting(){
    if(_limit){
#ifdef USE_CPU
      if(_Lag2Mono) delete[] _Lag2Mono;
      if(_Mono2Lag) delete[] _Mono2Lag;
      if(_V1D)      delete[] _V1D;
      if(_weight)   delete[] _weight;
      if(_A)        delete[] _A;
      if(_Alim)     delete[] _Alim;
#elif USE_GPU
      if(_Lag2Mono) delete[] _Lag2Mono;
      if(_Mono2Lag) delete[] _Mono2Lag;
      if(_V1D)      delete[] _V1D;
      if(_weight)   delete[] _weight;
      if(_A)        delete[] _A;
      if(_Alim)     delete[] _Alim;
#endif      
    }	      
  }
  
  bool getLimitingStatus() const {return _limit;}

  void HRlimiting(scalar* U){
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);
  }

  void copyMatrixToPointer(fullMatrix<scalar> &A, scalar* h_A){
  
    // Column major sorting
    for(int j = 0; j < A.size2(); j++){
      for(int i = 0; i < A.size1(); i++){
	h_A[j*A.size1()+i] = A(i,j);
      }
    }
  }

  void copyMatrixToPointer(const fullMatrix<scalar> &A, scalar* h_A){
  
    // Column major sorting
    for(int j = 0; j < A.size2(); j++){
      for(int i = 0; i < A.size1(); i++){
	h_A[j*A.size1()+i] = A(i,j);
      }
    }
  }

};
#endif
