#ifndef LIMITING_H
#define LIMITING_H

#include <misc.h>

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
  int     _method; // no limit= 0; HR=1
  int     _N_s;
  int     _N_E;
  int     _N_F;
  int     _N_G;
  int     _boundaryMap;
  
 public:
  // constructor
 Limiting(int method) : _method(method){}
  
 Limiting(int method, int N_s, int N_E, int N_F, int N_G, int boundaryMap,
	  fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, fullMatrix<scalar> &V1D, scalar* weight)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_F(N_F), _N_G(N_G), _boundaryMap(boundaryMap){
    switch (_method){
    case 1:{
#ifdef USE_CPU
      _Lag2Mono = new scalar[_N_s*_N_s];     copyMatrixToPointer(Lag2Mono,_Lag2Mono);
      _Mono2Lag = new scalar[_N_s*_N_s];     copyMatrixToPointer(Mono2Lag,_Mono2Lag);
      _V1D      = new scalar[_N_G*_N_s];     copyMatrixToPointer(V1D,_V1D);
      _weight   = new scalar[_N_G];          memcpy(_weight,weight,N_G*sizeof(scalar));
      _A        = new scalar[_N_s*_N_E*N_F];
      _Alim     = new scalar[_N_s*_N_E*N_F];
#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2Mono = new scalar[_N_s*_N_s];     copyMatrixToPointer(Lag2Mono,tmpLag2Mono);
      scalar* tmpMono2Lag = new scalar[_N_s*_N_s];     copyMatrixToPointer(Mono2Lag,tmpMono2Lag);
      scalar* tmpV1D      = new scalar[_N_G*_N_s];     copyMatrixToPointer(V1D,tmpV1D);

      // Allocate on GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Lag2Mono,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Mono2Lag,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_V1D,_N_G*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_weight,_N_G*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_A,_N_s*_N_E*_N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Alim,_N_s*_N_E*_N_F*sizeof(scalar)));

      // Copy data to GPU
      CUDA_SAFE_CALL(cudaMemcpy(_Lag2Mono, tmpLag2Mono, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_Mono2Lag, tmpMono2Lag, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_V1D, tmpV1D, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_weight,weight, N_G*sizeof(scalar), cudaMemcpyHostToDevice));
	    
      delete[] tmpLag2Mono;	
      delete[] tmpMono2Lag;	
      delete[] tmpV1D;	
#endif
    }
      break;
    default:
      printf("No limiting.\n");
      _method = 0;
    }
  }
  
  // destructor
  ~Limiting(){
    if(_method!=0){
#ifdef USE_CPU
      if(_Lag2Mono) delete[] _Lag2Mono;
      if(_Mono2Lag) delete[] _Mono2Lag;
      if(_V1D)      delete[] _V1D;
      if(_weight)   delete[] _weight;
      if(_A)        delete[] _A;
      if(_Alim)     delete[] _Alim;
#elif USE_GPU
      if(_Lag2Mono) CUDA_SAFE_CALL(cudaFree(_Lag2Mono));
      if(_Mono2Lag) CUDA_SAFE_CALL(cudaFree(_Mono2Lag));
      if(_V1D)      CUDA_SAFE_CALL(cudaFree(_V1D));
      if(_weight)   CUDA_SAFE_CALL(cudaFree(_weight));
      if(_A)        CUDA_SAFE_CALL(cudaFree(_A));
      if(_Alim)     CUDA_SAFE_CALL(cudaFree(_Alim));
#endif      
    }	      
  }
  
  int getLimitingMethod() const {return _method;}

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
