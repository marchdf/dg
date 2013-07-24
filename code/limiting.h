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
  scalar* _Lag2MonoX;
  scalar* _MonoX2MonoY;
  scalar* _MonoY2Lag;
  scalar* _XYZCen;
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
  bool    _cartesian;
  int     _D;
  int     _N_s;
  int     _N_E;
  int     _N_G;
  int     _N_N;
  int     _N_ghosts;
  int     _L;  // size of TaylorDxIdx
  int     _order;
  int     _N_s1D;
  int     _N_G1D;
  int     _L2Msize1;
  int     _L2Msize2;
  int     _boundaryMap;
  int*    _neighbors;
  int*    _ghostElementSend;
  int*    _ghostElementRecv;
  int*    _TaylorDxIdx;
  int*    _TaylorDyIdx;
  scalar _refArea;
  
 public:
  // constructor
 Limiting(int method,bool cartesian) : _method(method), _cartesian(cartesian){}

  // 1D limiting constructor
 Limiting(int method, int N_s, int N_E, int N_G, int N_N, int* neighbors, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, fullMatrix<scalar> &V1D, scalar* weight)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_G(N_G), _N_N(N_N){

    _Lag2Mono=NULL;
    _Mono2Lag=NULL;
    _Lag2MonoX=NULL;
    _MonoX2MonoY=NULL;
    _MonoY2Lag=NULL;
    _XYZCen=NULL;
    _powersXYZG=NULL;
    _V1D=NULL;
    _weight=NULL; // integration weights
    _A=NULL;      // monomial solution
    _Alim=NULL;   // holds the limited monomial solution
    _pressure=NULL;
    _pressureMono=NULL;
    _pressureLim=NULL;
    _u=NULL;
    _uMono=NULL;
    _uLim=NULL;
    _neighbors=NULL;
    _ghostElementSend=NULL;
    _ghostElementRecv=NULL;
    _TaylorDxIdx=NULL;
    _TaylorDyIdx=NULL;

    switch (_method){
    case 1:
    case 2:
    case 3:{
#ifdef USE_CPU
      _Lag2Mono = new scalar[_N_s*_N_s];     Lag2Mono.copyMatrixToPointer(_Lag2Mono);
      _Mono2Lag = new scalar[_N_s*_N_s];     Mono2Lag.copyMatrixToPointer(_Mono2Lag);
      _V1D      = new scalar[_N_G*_N_s];     V1D.copyMatrixToPointer(_V1D);
      _weight   = new scalar[_N_G];          memcpy(_weight,weight,N_G*sizeof(scalar));
      _A        = new scalar[_N_s*_N_E*N_F]; 
      _Alim     = new scalar[_N_s*_N_E*N_F]; 
      _neighbors  = new int[_N_N*_N_E];
      memcpy(_neighbors,   neighbors,   _N_N*_N_E*sizeof(int));

#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2Mono = new scalar[_N_s*_N_s];     Lag2Mono.copyMatrixToPointer(tmpLag2Mono);
      scalar* tmpMono2Lag = new scalar[_N_s*_N_s];     Mono2Lag.copyMatrixToPointer(tmpMono2Lag);
      scalar* tmpV1D      = new scalar[_N_G*_N_s];     V1D.copyMatrixToPointer(tmpV1D);

      // Allocate on GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Lag2Mono,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Mono2Lag,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_V1D,_N_G*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_weight,_N_G*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_A,_N_s*_N_E*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Alim,_N_s*_N_E*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_neighbors  ,_N_N*_N_E*sizeof(int)));

      // Copy data to GPU
      CUDA_SAFE_CALL(cudaMemcpy(_Lag2Mono, tmpLag2Mono, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_Mono2Lag, tmpMono2Lag, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_V1D, tmpV1D, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_weight,weight, N_G*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_neighbors,  neighbors,_N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice));
      
      delete[] tmpLag2Mono;	
      delete[] tmpMono2Lag;	
      delete[] tmpV1D;	
#endif
      }
      break;
    default:
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

  // 2D limiting constructor for structured mesh
 Limiting(int method, int N_s, int N_E, int N_G, int order, bool cartesian, int N_N, int M_ghosts, int* neighbors, fullMatrix<scalar> &Lag2MonoX, fullMatrix<scalar> &MonoX2MonoY, fullMatrix<scalar> &MonoY2Lag, fullMatrix<scalar> &V1D, int* ghostElementSend, int* ghostElementRecv, scalar* weight) : _method(method), _N_s(N_s), _N_E(N_E), _N_G(N_G), _order(order), _cartesian(cartesian), _N_N(N_N), _N_ghosts(M_ghosts){

    _N_s1D = (order+1);
    _N_G1D = (order+1), 
    _Lag2Mono=NULL;
    _Mono2Lag=NULL;
    _Lag2MonoX=NULL;
    _MonoX2MonoY=NULL;
    _MonoY2Lag=NULL;
    _XYZCen=NULL;
    _powersXYZG=NULL;
    _V1D=NULL;
    _weight=NULL; // integration weights
    _A=NULL;      // monomial solution
    _Alim=NULL;   // holds the limited monomial solution
    _pressure=NULL;
    _pressureMono=NULL;
    _pressureLim=NULL;
    _u=NULL;
    _uMono=NULL;
    _uLim=NULL;
    _neighbors=NULL;
    _ghostElementSend=NULL;
    _ghostElementRecv=NULL;
    _TaylorDxIdx=NULL;
    _TaylorDyIdx=NULL;

    switch (_method){
    case 1:
    case 2:
    case 3:{
#ifdef USE_CPU
      _Lag2MonoX   = new scalar[_N_s*_N_s];     Lag2MonoX.copyMatrixToPointer(_Lag2MonoX);
      _MonoX2MonoY = new scalar[_N_s*_N_s];     MonoX2MonoY.copyMatrixToPointer(_MonoX2MonoY);
      _MonoY2Lag   = new scalar[_N_s*_N_s];     MonoY2Lag.copyMatrixToPointer(_MonoY2Lag);
      _V1D      = new scalar[_N_G1D*_N_s1D];    V1D.copyMatrixToPointer(_V1D);
      _weight   = new scalar[_N_G1D];           memcpy(_weight,weight,_N_G1D*sizeof(scalar));
      _A        = new scalar[_N_s*(_N_E+_N_ghosts)*N_F];    makeZero(_A,     _N_s*(_N_E+_N_ghosts)*N_F);
      _Alim     = new scalar[_N_s*(_N_E+_N_ghosts)*N_F];    makeZero(_Alim,  _N_s*(_N_E+_N_ghosts)*N_F);
      _neighbors  = new int[_N_N*_N_E]; memcpy(_neighbors,   neighbors,   _N_N*_N_E*sizeof(int));
      _ghostElementSend  = new int[_N_ghosts*3]; memcpy(_ghostElementSend,   ghostElementSend,   _N_ghosts*3*sizeof(int));
      _ghostElementRecv  = new int[_N_ghosts*3]; memcpy(_ghostElementRecv,   ghostElementRecv,   _N_ghosts*3*sizeof(int));
     
#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2MonoX   = new scalar[_N_s*_N_s];     Lag2MonoX.copyMatrixToPointer(tmpLag2MonoX);
      scalar* tmpMonoX2MonoY = new scalar[_N_s*_N_s];     MonoX2MonoY.copyMatrixToPointer(tmpMonoX2MonoY);
      scalar* tmpMonoY2Lag   = new scalar[_N_s*_N_s];     MonoY2Lag.copyMatrixToPointer(tmpMonoY2Lag);
      scalar* tmpV1D         = new scalar[_N_G1D*_N_s1D]; V1D.copyMatrixToPointer(tmpV1D);

      // Allocate on GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Lag2MonoX  ,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_MonoX2MonoY,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_MonoY2Lag  ,_N_s*_N_s*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_V1D        ,_N_G1D*_N_s1D*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_weight     ,_N_G1D*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_A          ,_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Alim       ,_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_neighbors  ,_N_N*_N_E*sizeof(int)));
	    
      // Copy data to GPU
      CUDA_SAFE_CALL(cudaMemcpy(_Lag2MonoX,   tmpLag2MonoX, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_MonoX2MonoY, tmpMonoX2MonoY, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_MonoY2Lag,   tmpMonoY2Lag, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_V1D, tmpV1D, _N_G1D*_N_s1D*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_weight,weight, _N_G1D*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_neighbors,  neighbors,_N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice));

      delete[] tmpLag2MonoX;
      delete[] tmpMonoX2MonoY;
      delete[] tmpMonoY2Lag;
      delete[] tmpV1D;
#endif
      }
      break;
    default:
      _method = 0;
    }
    
    /* // For case specific stuff: */
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
  }// end 2D constructor for structured mesh
  
  // 2D limiting constructor for unstructured mesh
 Limiting(int method, int N_s, int N_E, int N_G, int N_N, int L, int order, int L2Msize1, int L2Msize2, int* neighbors, fullMatrix<scalar> Lag2Mono, fullMatrix<scalar> Mono2Lag, fullMatrix<scalar> XYZCen, scalar* powersXYZG, scalar* weight, scalar refArea, int* TaylorDxIdx, int* TaylorDyIdx)
   : _method(method), _D(D),_N_s(N_s), _N_E(N_E), _N_G(N_G), _N_N(N_N), _L(L), _order(order), _L2Msize1(L2Msize1), _L2Msize2(L2Msize2), _refArea(refArea){

    _Lag2Mono=NULL;
    _Mono2Lag=NULL;
    _Lag2MonoX=NULL;
    _MonoX2MonoY=NULL;
    _MonoY2Lag=NULL;
    _XYZCen=NULL;
    _powersXYZG=NULL;
    _V1D=NULL;
    _weight=NULL; // integration weights
    _A=NULL;      // monomial solution
    _Alim=NULL;   // holds the limited monomial solution
    _pressure=NULL;
    _pressureMono=NULL;
    _pressureLim=NULL;
    _u=NULL;
    _uMono=NULL;
    _uLim=NULL;
    _neighbors=NULL;
    _TaylorDxIdx=NULL;
    _TaylorDyIdx=NULL;
    
    switch (_method){
    case 1:
    case 2:
    case 3:{
#ifdef USE_CPU
      // Allocate
      _Lag2Mono   = new scalar[_L2Msize1*_L2Msize2*_N_E];
      _Mono2Lag   = new scalar[_L2Msize2*_L2Msize1*_N_E];
      _XYZCen     = new scalar[_N_E*_D];
      _powersXYZG = new scalar[_L2Msize1*_N_G*(_N_N+1)*_N_E];
      _neighbors  = new int[_N_N*_N_E];
      _weight     = new scalar[_N_G];
      _TaylorDxIdx= new int[_L];
      _TaylorDyIdx= new int[_L];
      _A          = new scalar[_L2Msize1*_N_E*N_F]; makeZero(_A,    _L2Msize1*_N_E*N_F);
      _Alim       = new scalar[_L2Msize1*_N_E*N_F]; makeZero(_Alim, _L2Msize1*_N_E*N_F);

      // Copy the data to these new pointers
      for(int e = 0; e < _N_E; e++){
	for(int alpha = 0; alpha < D; alpha++){ _XYZCen[e*D+alpha] = XYZCen(e,alpha);}
      	for(int i = 0; i < _L2Msize1; i++){
      	  for(int j = 0; j < _L2Msize2; j++){
      	    _Lag2Mono[(e*_L2Msize1+i)*_L2Msize2+j] = Lag2Mono(e,i*_L2Msize2+j);
      	    _Mono2Lag[(e*_L2Msize2+j)*_L2Msize1+i] = Mono2Lag(e,j*_L2Msize1+i);}}}
      memcpy(_powersXYZG,  powersXYZG,  _L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar));
      memcpy(_neighbors,   neighbors,   _N_N*_N_E*sizeof(int));
      memcpy(_weight,      weight,      _N_G*sizeof(scalar));
      memcpy(_TaylorDxIdx, TaylorDxIdx, _L*sizeof(int));
      memcpy(_TaylorDyIdx, TaylorDyIdx, _L*sizeof(int));
      
#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2Mono   = new scalar[_L2Msize1*_L2Msize2*_N_E];
      scalar* tmpMono2Lag   = new scalar[_L2Msize2*_L2Msize1*_N_E];
      scalar* tmpXYZCen     = new scalar[_N_E*_D];
      for(int e = 0; e < _N_E; e++){
	for(int alpha = 0; alpha < D; alpha++){ tmpXYZCen[e*D+alpha] = XYZCen(e,alpha);}
	for(int i = 0; i < _L2Msize1; i++){
	  for(int j = 0; j < _L2Msize2; j++){
      	    tmpLag2Mono[(e*_L2Msize1+i)*_L2Msize2+j] = Lag2Mono(e,i*_L2Msize2+j);
      	    tmpMono2Lag[(e*_L2Msize2+j)*_L2Msize1+i] = Mono2Lag(e,j*_L2Msize1+i);}}}

      // Allocate on GPU
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Lag2Mono,_L2Msize1*_L2Msize2*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Mono2Lag,_L2Msize2*_L2Msize1*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_XYZCen,_D*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_powersXYZG,_L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_neighbors,_N_N*_N_E*sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_weight,_N_G*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_TaylorDxIdx,_L*sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_TaylorDyIdx,_L*sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_A,_L2Msize1*_N_E*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMalloc((void**) &_Alim,_L2Msize1*_N_E*N_F*sizeof(scalar)));

      // Copy data to GPU
      CUDA_SAFE_CALL(cudaMemcpy(_Lag2Mono,   tmpLag2Mono, _L2Msize1*_L2Msize2*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_Mono2Lag,   tmpMono2Lag, _L2Msize2*_L2Msize1*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_XYZCen,     tmpXYZCen,   _D*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_powersXYZG, powersXYZG,  _L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_neighbors,  neighbors,   _N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_weight,     weight,      _N_G*sizeof(scalar), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_TaylorDxIdx,TaylorDxIdx, _L*sizeof(int), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(_TaylorDyIdx,TaylorDyIdx, _L*sizeof(int), cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemset(_A,    (scalar)0.0, _L2Msize1*_N_E*N_F*sizeof(scalar)));
      CUDA_SAFE_CALL(cudaMemset(_Alim, (scalar)0.0, _L2Msize1*_N_E*N_F*sizeof(scalar)));
      
      delete[] tmpLag2Mono;
      delete[] tmpMono2Lag;
      delete[] tmpXYZCen;
#endif
      }
      break;
    default:
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
    if(_Lag2Mono)     del(_Lag2Mono);
    if(_Mono2Lag)     del(_Mono2Lag);
    if(_Lag2MonoX)    del(_Lag2MonoX);
    if(_MonoX2MonoY)  del(_MonoX2MonoY);
    if(_MonoY2Lag)    del(_MonoY2Lag);
    if(_weight)       del(_weight);
    if(_A)            del(_A);
    if(_Alim)         del(_Alim);
    if(_V1D)          del(_V1D);
    if(_XYZCen)       del(_XYZCen);
    if(_powersXYZG)   del(_powersXYZG);
    if(_neighbors)    del(_neighbors);
    if(_ghostElementRecv) del(_ghostElementRecv);
    if(_ghostElementSend) del(_ghostElementSend);
    if(_TaylorDxIdx)  del(_TaylorDxIdx);
    if(_TaylorDyIdx)  del(_TaylorDyIdx);
    if(_pressure)     del(_pressure);
    if(_pressureMono) del(_pressureMono);
    if(_pressureLim)  del(_pressureLim);
    if(_u)            del(_u);
    if(_uMono)        del(_uMono);
    if(_uLim)         del(_uLim);
  }	      
  
  int getLimitingMethod() const {return _method;}
  
  void HRlimiting(scalar* U){
#ifdef ONED
    // Go from lagrange to monomial representation
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _A, _Alim);
    
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

#elif TWOD
    if(_cartesian){
      // Go from lagrange to monomial representation wrt x
      blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Lag2MonoX, _N_s, U, _N_s, 0.0, _A, _N_s);

      // Communicate the elements on different partitions
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      Lcpu_CommunicateGhosts(_N_s, _N_E, _N_ghosts, _ghostElementSend, _ghostElementRecv, _A);
#endif
      
      // Limit the solution according to Liu (for each x slice)
      Lcpu_hrl1D(_N_s1D, _N_E, _N_G1D, _N_N, _N_s1D, _neighbors, 0, _weight, _V1D, _A, _Alim);
      
      // Go to the monomial representation wrt y
      blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _MonoX2MonoY, _N_s, _Alim, _N_s, 0.0, _A, _N_s);

      // Communicate the elements on different partitions
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      Lcpu_CommunicateGhosts(_N_s, _N_E, _N_ghosts, _ghostElementSend, _ghostElementRecv, _A);
#endif

      // Limit the solution according to Liu (for each y slice)
      Lcpu_hrl1D(_N_s1D, _N_E, _N_G1D, _N_N, _N_s1D, _neighbors, 2, _weight, _V1D, _A, _Alim);
      
      // Go back to lagrange
      blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _MonoY2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);
      
    }
    
    if(!_cartesian){
      /* for(int e=0; e<_N_E;e++){ */
      /*   for(int fc=0;fc<N_F;fc++){ */
      /* 	U[(e*N_F+fc)*_N_s+0] = 2; */
      /* 	U[(e*N_F+fc)*_N_s+1] = 3; */
      /* 	U[(e*N_F+fc)*_N_s+2] = 3; */
      /* 	U[(e*N_F+fc)*_N_s+3] = 3; */
      /* 	U[(e*N_F+fc)*_N_s+4] = 3; */
      /* 	U[(e*N_F+fc)*_N_s+5] = 3; */
      /*   } */
      /* } */
      /* printf("Before limiting:\n"); */
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("U(e=%i,i=%i) = %f\n",8,i,U[(8*N_F+0)*_N_s+i]); */
      /* } */

      // Go from lagrange to monomial representation
      LChangeBasis(_L2Msize1, _L2Msize2, _N_E, _Lag2Mono, U, _A);
      
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("A(e=%i,i=%i) = %.12f\n",8,i,_A[(8*N_F+0)*_N_s+i]); */
      /* } */
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("A(e=%i,i=%i) = %.12f\n",5,i,_A[(5*N_F+0)*_N_s+i]); */
      /* } */
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("A(e=%i,i=%i) = %.12f\n",13,i,_A[(13*N_F+0)*_N_s+i]); */
      /* } */
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("A(e=%i,i=%i) = %.12f\n",1,i,_A[(1*N_F+0)*_N_s+i]); */
      /* } */
  
      // Limit the solution according to Liu
      Lcpu_hrl2D(_L2Msize1, _N_E, _N_G, _N_N, _order, _XYZCen, _powersXYZG, _neighbors, _TaylorDxIdx, _TaylorDyIdx, _weight, _refArea, _A, _Alim);

      /* printf("After limiting:\n"); */
      /* for(int i = 0; i < _N_s; i++){ */
      /* 	printf("Alim(e=%i,i=%i) = %f\n",8,i,_Alim[(8*N_F+0)*_N_s+i]); */
      /* } */
      // Go back to lagrange representation
      LChangeBasis(_L2Msize2, _L2Msize1, _N_E, _Mono2Lag, _Alim, U);
      /* for(int i = 0; i < _N_s; i++){ */
      /*   printf("U(e=%i,i=%i) = %f\n",8,i,U[(8*N_F+0)*_N_s+i]); */
      /* } */
    }
#endif

  }

  void MYlimiting(scalar* U){

    // Get the pressure field
    Lpressure(_N_s, _N_E, U, _pressure);

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _pressureMono, _pressureLim);

    // Go from lagrange to monomial representation
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _A, _Alim);

    // My modification
    Llimmodif(_N_s, _N_E, _A, _pressureLim, _Alim);
    
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

  }

  void M2limiting(scalar* U){

    // Get the pressure and velocity fields
    Lpressure_u(_N_s, _N_E, U, _pressure, _u);

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _pressure, _N_s, 0.0, _pressureMono, _N_s);
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _pressureMono, _pressureLim);

    // Limit velocity
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, _u, _N_s, 0.0, _uMono, _N_s);
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _uMono, _uLim);

    // Limit the other variables
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    Lcpu_hrl1D(_N_s, _N_E, _N_G, _N_N, 1, _neighbors, 0, _weight, _V1D, _A, _Alim);

    // My modification
    Llimmodif2(_N_s, _N_E, _A, _pressureLim, _uLim, _Alim);
    
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

  }
};
#endif
