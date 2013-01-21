#ifndef LIMITING_H
#define LIMITING_H

#include <misc.h>
#include <constants.h>
#include <physics.h>

class Limiting
{
 private:
  scalar* _Lag2Mono;
  scalar* _Mono2Lag;
  scalar* _powersXYZG;
  scalar* _V1D;
  scalar* _weight; // integration weights
  scalar* _A;      // monomial solution
  scalar* _Alim;   // holds the limited monomial solution
  scalar* _pressure;
  scalar* _pressureMono;
  scalar* _pressureLim;
  scalar* _u;
  scalar* _uMono;
  scalar* _uLim;
  int     _method; // no limit= 0; HR=1
  int     _N_s;
  int     _N_E;
  int     _N_F;
  int     _N_G;
  int     _boundaryMap;
  int*    _neighbors;
  
 public:
  // constructor
 Limiting(int method) : _method(method){}

  // 1D limiting constructor
 Limiting(int method, int N_s, int N_E, int N_F, int N_G, int boundaryMap, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, fullMatrix<scalar> &V1D, scalar* weight)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_F(N_F), _N_G(N_G), _boundaryMap(boundaryMap){
    switch (_method){
    case 1:
    case 2:
    case 3:{
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
    
    // For case specific stuff:
    switch (_method){
    case 2:{
#ifdef USE_CPU
      _pressure     = new scalar[_N_s*_N_E];
      _pressureMono = new scalar[_N_s*_N_E];
      _pressureLim  = new scalar[_N_s*_N_E];
#elif USE_GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressure,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureMono,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar)));
#endif
      }
      break;
    case 3:{
#ifdef USE_CPU
      _pressure     = new scalar[_N_s*_N_E];
      _pressureMono = new scalar[_N_s*_N_E];
      _pressureLim  = new scalar[_N_s*_N_E];
      _u            = new scalar[_N_s*_N_E];
      _uMono        = new scalar[_N_s*_N_E];
      _uLim         = new scalar[_N_s*_N_E];
#elif USE_GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressure,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureMono,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_u,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_uMono,_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_uLim,_N_s*_N_E*sizeof(scalar)));
#endif
      }
      break;
    }
  } // end 1D constructor

  // 2D limiting constructor
 Limiting(int method, int N_s, int N_E, int N_F, int N_G, int* neighbors, fullMatrix<scalar> Lag2Mono, fullMatrix<scalar> Mono2Lag, scalar* powersXYZG, scalar* weight)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_F(N_F), _N_G(N_G){
    switch (_method){
    case 1:
    case 2:
    case 3:{
#ifdef USE_CPU
      // Allocate
      _Lag2Mono   = new scalar[_N_s*_N_s*_N_E];
      _Mono2Lag   = new scalar[_N_s*_N_s*_N_E];
      _powersXYZG = new scalar[_N_s*_N_G*(_N_N+1)*_N_E];
      _neigbors   = new int[_N_N*_N_E];
      _weight     = new scalar[_N_G];          
      _A          = new scalar[_N_s*_N_E*_N_F];
      _Alim       = new scalar[_N_s*_N_E*_N_F];

      // Copy the data to these new pointers
      for(int e = 0; e < _N_E; e++){
	for(int i = 0; i < _N_s; i++){
	  for(int ii = 0; ii < _N_s; ii++){
	    _Lag2Mono[(e*_N_s+i)*_N_s+ii] = Lag2Mono(e,i*_N_s+ii);
	    _Mono2Lag[(e*_N_s+i)*_N_s+ii] = Mono2Lag(e,i*_N_s+ii);}}}
      memcpy(_powersXYZG, powersXYZG, _N_s*_N_G*(_N_N+1)*_N_E*sizeof(scalar));
      memcpy(_neighbors,  neighbors,  _N_N*_N_E*sizeof(scalar));
      memcpy(_weight,     weight,     _N_G*sizeof(scalar));
      
#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2Mono   = new scalar[_N_s*_N_s*_N_E];
      scalar* tmpMono2Lag   = new scalar[_N_s*_N_s*_N_E];
      for(int e = 0; e < _N_E; e++){
	for(int i = 0; i < _N_s; i++){
	  for(int ii = 0; ii < _N_s; ii++){
	    _Lag2Mono[(e*_N_s+i)*_N_s+ii] = Lag2Mono(e,i*_N_s+ii);
	    _Mono2Lag[(e*_N_s+i)*_N_s+ii] = Mono2Lag(e,i*_N_s+ii);}}}

      // Allocate on GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Lag2Mono,_N_s*_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Mono2Lag,_N_s*_N_s*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_powersXYZG,_N_s*_N_G*(_N_N+1)*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_neighbors,_N_N*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_weight,_N_G*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_A,_N_s*_N_E*_N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Alim,_N_s*_N_E*_N_F*sizeof(scalar)));

      // Copy data to GPU
      CUDA_SAFE_CALL(cudaMemcpy(_Lag2Mono,   tmpLag2Mono, _N_s*_N_s*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_Mono2Lag,   tmpMono2Lag, _N_s*_N_s*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_powersXYZG, powersXYZG,  _N_s*_N_G*(_N_N+1)*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_neighbors,  neighbors,   _N_N*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_weight,     weight,      _N_G*sizeof(scalar), cudaMemcpyHostToDevice));
	    
      delete[] tmpLag2Mono;
      delete[] tmpMono2Lag;
#endif
      }
      break;
    default:
      printf("No limiting.\n");
      _method = 0;
    }
    
/*     // For case specific stuff: */
/*     switch (_method){ */
/*     case 2:{ */
/* #ifdef USE_CPU */
/*       _pressure     = new scalar[_N_s*_N_E]; */
/*       _pressureMono = new scalar[_N_s*_N_E]; */
/*       _pressureLim  = new scalar[_N_s*_N_E]; */
/* #elif USE_GPU */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressure,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureMono,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar))); */
/* #endif */
/*       } */
/*       break; */
/*     case 3:{ */
/* #ifdef USE_CPU */
/*       _pressure     = new scalar[_N_s*_N_E]; */
/*       _pressureMono = new scalar[_N_s*_N_E]; */
/*       _pressureLim  = new scalar[_N_s*_N_E]; */
/*       _u            = new scalar[_N_s*_N_E]; */
/*       _uMono        = new scalar[_N_s*_N_E]; */
/*       _uLim         = new scalar[_N_s*_N_E]; */
/* #elif USE_GPU */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressure,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureMono,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_u,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_uMono,_N_s*_N_E*sizeof(scalar))); */
/*       CUDA_SAFE_CALL(cudaMalloc((void**) &_uLim,_N_s*_N_E*sizeof(scalar))); */
/* #endif */
/*       } */
/*       break; */
/*     } */
  } // end 2D constructor
  
  // destructor
  ~Limiting(){
    if(_method!=0){
      if(_Lag2Mono)     del(_Lag2Mono);
      if(_Mono2Lag)     del(_Mono2Lag);
      if(_weight)       del(_weight);
      if(_A)            del(_A);
      if(_Alim)         del(_Alim);
#ifdef ONED
      if(_V1D)          del(_V1D);
#elif TWOD
      if(_powersXYZG)   del(_powersXYZG);
      if(_neighbors)    del(_neigbors);
#endif
    }
    if((_method==2)||(_method==3)){
      if(_pressure)     del(_pressure);
      if(_pressureMono) del(_pressureMono);
      if(_pressureLim)  del(_pressureLim);
    }
    if(_method==3){
      if(_u)            del(_u);
      if(_uMono)        del(_uMono);
      if(_uLim)         del(_uLim);
    }
  }	      
  
  int getLimitingMethod() const {return _method;}

  void HRlimiting(scalar* U){
    // Change of variables
    //Lcpu_Cons2Prim(_N_s, _N_E, _N_F, U, _multifluid, _passive, _model, _gamma0);
    //Lcpu_Cons2Half(_N_s, _N_E, _N_F, U, _multifluid, _passive, _model, _gamma0);

    // Go from lagrange to monomial representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

    // Go back to initial variables
    //Lcpu_Prim2Cons(_N_s, _N_E, _N_F, U, _multifluid, _passive, _model, _gamma0);
    //Lcpu_Half2Cons(_N_s, _N_E, _N_F, U, _multifluid, _passive, _model, _gamma0);
  }

  void MYlimiting(scalar* U){

    // Get the pressure field
    Lpressure(_N_s, _N_E, _N_F, U, _pressure);

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, _pressureMono, _pressureLim);

    // Go from lagrange to monomial representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);

    // My modification
    Llimmodif(_N_s, _N_E, _N_F, _A, _pressureLim, _Alim);
    
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

  }

  void M2limiting(scalar* U){

    // Get the pressure and velocity fields
    Lpressure_u(_N_s, _N_E, _N_F, U, _pressure, _u);

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, _pressureMono, _pressureLim);

    // Limit velocity
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _u, _N_s, 0.0, _uMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, _uMono, _uLim);

    // Limit the other variables
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);

    // My modification
    Llimmodif2(_N_s, _N_E, _N_F, _A, _pressureLim, _uLim, _Alim);
    
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
