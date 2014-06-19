/*!
  \file limiting.h
  \brief Class deals with solution limiting.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup limiting Limiting
  \ingroup limiting
*/
#ifndef LIMITING_H
#define LIMITING_H

#include <misc.h>
#include <constants.h>
#include <physics.h>
#include <limiting_kernels.h>
#include <communicator.h>
#include <sensor.h>
#include "simpleMesh.h"
#include "timers.h"
#include "mem_counter.h"
#ifdef USE_GPU
#include <cublas.h>
#endif

// Used to define dynamically variables (mass fractions)
#define _Y(x) _Y ##x
#define _YMono(x) _YMono ##x
#define _YLim(x) _YLim ##x

class Limiting
{
 private:
  scalar* _Lag2Mono;
  scalar* _Mono2Lag;
  scalar* _Lag2MonoX;
  scalar* _Lag2MonoY;
  scalar* _MonoX2MonoY;
  scalar* _MonoX2Lag;
  scalar* _MonoY2Lag;
  scalar* _XYZCen;
  scalar* _powersXYZG;

  scalar* _A;      // monomial solution
  scalar* _Alim;   // limited monomial solution
  scalar* _Utmp;   // for hri/m2i, need a tmp to store limited solution
  
  scalar* _rho, * _rhoMono, * _rhoLim;
  scalar* _rhou, * _rhouMono, * _rhouLim;
  scalar* _rhov, * _rhovMono, * _rhovLim;
  scalar* _pressure, * _pressureMono, * _pressureLim;
  scalar* _K, * _KLim;
  scalar* _E, * _EMono, * _ELim;
  scalar* _gamma, * _gammaMono, * _gammaLim;
  scalar* _beta, * _betaMono, * _betaLim; //=gamma*pinf/(gamma-1)
  scalar* _rhoeLim;
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) scalar* _Y(x), * _YMono(x), * _YLim(x);
#include "loop.h"
  
  int     _method; // no limiting= 0; HR=1; M2L=2
  bool    _cartesian;
  int     _D;
  int     _N_s;
  int     _N_E;
  int     _N_N;
  int     _N_ghosts;
  int     _L;  // size of TaylorDxIdx
  int     _order;
  int     _N_s1D;
  int     _L2Msize1;
  int     _L2Msize2;
  int     _boundaryMap;
  int*    _neighbors;
  int*    _TaylorDxIdx;
  int*    _TaylorDyIdx;
  scalar _refArea;

  TIMERS &_timers;
  
  void common_ctor(){
    /*\brief Common contructor*/
    _Lag2Mono=NULL;
    _Mono2Lag=NULL;
    _Lag2MonoX=NULL;
    _Lag2MonoY=NULL;
    _MonoX2MonoY=NULL;
    _MonoX2Lag=NULL;
    _MonoY2Lag=NULL;
    _XYZCen=NULL;
    _powersXYZG=NULL;
    _A=NULL;      // monomial solution
    _Alim=NULL;   // holds the limited monomial solution
    _Utmp=NULL;
    
    _rho=NULL; _rhoMono=NULL; _rhoLim=NULL;
    _rhou=NULL; _rhouMono=NULL; _rhouLim=NULL;
    _rhov=NULL; _rhovMono=NULL; _rhovLim=NULL;
    _pressure=NULL; _pressureMono=NULL; _pressureLim=NULL;
    _K=NULL; _KLim=NULL;
    _E=NULL; _EMono=NULL, _ELim=NULL;
    _gamma=NULL; _gammaMono=NULL; _gammaLim=NULL;
    _beta=NULL; _betaMono=NULL; _betaLim=NULL;
    _rhoeLim=NULL;
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) _Y(x) = NULL; _YMono(x) = NULL; _YLim(x) = NULL;
#include "loop.h"

    _neighbors=NULL;
    _TaylorDxIdx=NULL;
    _TaylorDyIdx=NULL;
  }
  
 public:
  /*\brief Constructor with method and bool for a cartesian mesh NOT USING ANYMORE*/
  //Limiting(int method,bool cartesian) : _method(method), _cartesian(cartesian){}

  /*!\brief Constructor for 1D limiting*/
 Limiting(int method, int N_s, int N_E, int N_N, simpleMesh &m, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, TIMERS &timers, MEM_COUNTER &mem_counter)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_N(N_N), _timers(timers){
    common_ctor();

    switch (_method){
    case 1:
    case 2:
    case 3:
    case 4:{
#ifdef USE_CPU
      _Lag2Mono = new scalar[_N_s*_N_s];     Lag2Mono.copyMatrixToPointer(_Lag2Mono);  mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _Mono2Lag = new scalar[_N_s*_N_s];     Mono2Lag.copyMatrixToPointer(_Mono2Lag);  mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _neighbors  = new int[_N_N*_N_E];                                                mem_counter.addToCPUCounter(_N_N*_N_E*sizeof(int));   
      memcpy(_neighbors,   m.getNeighbors(),   _N_N*_N_E*sizeof(int));

#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2Mono = new scalar[_N_s*_N_s];     Lag2Mono.copyMatrixToPointer(tmpLag2Mono);
      scalar* tmpMono2Lag = new scalar[_N_s*_N_s];     Mono2Lag.copyMatrixToPointer(tmpMono2Lag);

      // Allocate on GPU
      cudaMalloc((void**) &_Lag2Mono,_N_s*_N_s*sizeof(scalar));                        mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_Mono2Lag,_N_s*_N_s*sizeof(scalar));			       mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_neighbors  ,_N_N*_N_E*sizeof(int));			       mem_counter.addToGPUCounter(_N_N*_N_E*sizeof(int));

      // Copy data to GPU
      cudaMemcpy(_Lag2Mono, tmpLag2Mono, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice);
      cudaMemcpy(_Mono2Lag, tmpMono2Lag, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice);
      cudaMemcpy(_neighbors,  m.getNeighbors(),_N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice);
      
      delete[] tmpLag2Mono;	
      delete[] tmpMono2Lag;	
#endif
      }
      break;
    default:
      _method = 0;
    }

    // For case specific stuff:
    switch (_method){
    case 1:{
#ifdef USE_CPU
      _A        = new scalar[_N_s*_N_E*N_F]; 
      _Alim     = new scalar[_N_s*_N_E*N_F];
      mem_counter.addToCPUCounter(2*_N_s*_N_E*N_F*sizeof(scalar));
#elif USE_GPU
      cudaMalloc((void**) &_A,_N_s*_N_E*N_F*sizeof(scalar));
      cudaMalloc((void**) &_Alim,_N_s*_N_E*N_F*sizeof(scalar));
      mem_counter.addToGPUCounter(2*_N_s*_N_E*N_F*sizeof(scalar));
#endif
    }
      break;
    case 2:{
#ifdef USE_CPU
      _rho      = new scalar[_N_s*_N_E]; _rhoMono      = new scalar[_N_s*_N_E]; _rhoLim      = new scalar[_N_s*_N_E];
      _rhou     = new scalar[_N_s*_N_E]; _rhouMono     = new scalar[_N_s*_N_E]; _rhouLim     = new scalar[_N_s*_N_E];
      _pressure = new scalar[_N_s*_N_E]; _pressureMono = new scalar[_N_s*_N_E]; _pressureLim = new scalar[_N_s*_N_E];
      _K        = new scalar[_N_s*N_E];  _KLim         = new scalar[_N_s*N_E];
      _E        = new scalar[_N_s*N_E];  _EMono        = new scalar[_N_s*N_E];  _ELim        = new scalar[_N_s*N_E]; makeZero(_ELim, _N_s*N_E);
      _gamma    = new scalar[_N_s*_N_E]; _gammaMono    = new scalar[_N_s*_N_E]; _gammaLim    = new scalar[_N_s*_N_E];
      _rhoeLim  = new scalar[_N_s*N_E];
      mem_counter.addToCPUCounter(18*_N_s*_N_E*sizeof(scalar));
#ifdef STIFFENED
      _beta    = new scalar[_N_s*_N_E]; _betaMono    = new scalar[_N_s*_N_E]; _betaLim    = new scalar[_N_s*_N_E];
      mem_counter.addToCPUCounter(3*_N_s*_N_E*sizeof(scalar));
#endif
      // Mass fractions initializations
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) _Y(x) = new scalar[_N_s*_N_E]; _YMono(x) = new scalar[_N_s*_N_E]; _YLim(x) = new scalar[_N_s*_N_E]; mem_counter.addToCPUCounter(3*_N_s*_N_E*sizeof(scalar));
#include "loop.h"
      
#elif USE_GPU
      cudaMalloc((void**) &_rho,_N_s*_N_E*sizeof(scalar));                 mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhoMono,_N_s*_N_E*sizeof(scalar));             mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhoLim,_N_s*_N_E*sizeof(scalar));              mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhou,_N_s*_N_E*sizeof(scalar));                mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhouMono,_N_s*_N_E*sizeof(scalar));            mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhouLim,_N_s*_N_E*sizeof(scalar));             mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_pressure,_N_s*_N_E*sizeof(scalar));            mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_pressureMono,_N_s*_N_E*sizeof(scalar));        mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar));         mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_K,_N_s*_N_E*sizeof(scalar));                   mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_KLim,_N_s*_N_E*sizeof(scalar));                mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));  
      cudaMalloc((void**) &_E,_N_s*_N_E*sizeof(scalar));                   mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_EMono,_N_s*_N_E*sizeof(scalar));               mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_ELim,_N_s*_N_E*sizeof(scalar));       cudaMemset(_ELim, (scalar)0.0, _N_s*N_E*sizeof(scalar)); mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_gamma,_N_s*_N_E*sizeof(scalar));               mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_gammaMono,_N_s*_N_E*sizeof(scalar));           mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_gammaLim,_N_s*_N_E*sizeof(scalar));            mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhoeLim,_N_s*_N_E*sizeof(scalar));             mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
#ifdef STIFFENED
      cudaMalloc((void**) &_beta,_N_s*_N_E*sizeof(scalar));                mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_betaMono,_N_s*_N_E*sizeof(scalar));            mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_betaLim,_N_s*_N_E*sizeof(scalar));             mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
#endif
      // Mass fractions initializations
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) cudaMalloc((void**) &_Y(x),_N_s*_N_E*sizeof(scalar)); \
      cudaMalloc((void**) &_YMono(x),_N_s*_N_E*sizeof(scalar));        \
      cudaMalloc((void**) &_YLim(x),_N_s*_N_E*sizeof(scalar));	       \
      mem_counter.addToGPUCounter(3*_N_s*_N_E*sizeof(scalar));
#include "loop.h"

#endif
    }
      break;
    case 3:
    case 4:{
#ifdef USE_CPU
      _Utmp     = new scalar[_N_s*_N_E*N_F];                        mem_counter.addToCPUCounter(_N_s*_N_E*N_F*sizeof(scalar));
#elif USE_GPU
      cudaMalloc((void**) &_Utmp,_N_s*_N_E*N_F*sizeof(scalar));     mem_counter.addToGPUCounter(_N_s*_N_E*N_F*sizeof(scalar));
#endif
    }
      break;
    }
  } // end 1D constructor

  /*!\brief Constructor for 2D limiting for structured mesh*/
 Limiting(int method, int N_s, int N_E, int order, bool cartesian, int N_N, int N_ghosts, simpleMesh &m, fullMatrix<scalar> &Lag2MonoX, fullMatrix<scalar> &MonoX2MonoY, fullMatrix<scalar> &MonoY2Lag, TIMERS &timers, MEM_COUNTER &mem_counter) : _method(method), _N_s(N_s), _N_E(N_E), _order(order), _cartesian(cartesian), _N_N(N_N), _N_ghosts(N_ghosts), _timers(timers){

    _N_s1D = (order+1);
    common_ctor();
    
    switch (_method){
    case 1:
    case 2:
    case 3:
    case 4:{
#ifdef USE_CPU
      _Lag2MonoX   = new scalar[_N_s*_N_s];     Lag2MonoX.copyMatrixToPointer(_Lag2MonoX);                         mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _Lag2MonoY   = new scalar[_N_s*_N_s];                                                                        mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _MonoX2MonoY = new scalar[_N_s*_N_s];     MonoX2MonoY.copyMatrixToPointer(_MonoX2MonoY);                     mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar)); 
      _MonoX2Lag   = new scalar[_N_s*_N_s];                                                                        mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _MonoY2Lag   = new scalar[_N_s*_N_s];     MonoY2Lag.copyMatrixToPointer(_MonoY2Lag);                         mem_counter.addToCPUCounter(_N_s*_N_s*sizeof(scalar));
      _neighbors  = new int[_N_N*_N_E];         memcpy(_neighbors,   m.getNeighbors(),   _N_N*_N_E*sizeof(int));   mem_counter.addToCPUCounter(_N_N*_N_E*sizeof(int));   
     
#elif USE_GPU
      // tmp host pointers to copy data to gpu
      scalar* tmpLag2MonoX   = new scalar[_N_s*_N_s];     Lag2MonoX.copyMatrixToPointer(tmpLag2MonoX);
      scalar* tmpMonoX2MonoY = new scalar[_N_s*_N_s];     MonoX2MonoY.copyMatrixToPointer(tmpMonoX2MonoY);
      scalar* tmpMonoY2Lag   = new scalar[_N_s*_N_s];     MonoY2Lag.copyMatrixToPointer(tmpMonoY2Lag);

      // Allocate on GPU
      cudaMalloc((void**) &_Lag2MonoX  ,_N_s*_N_s*sizeof(scalar));                                                 mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_Lag2MonoY  ,_N_s*_N_s*sizeof(scalar));						   mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_MonoX2MonoY,_N_s*_N_s*sizeof(scalar));						   mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_MonoX2Lag  ,_N_s*_N_s*sizeof(scalar));						   mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_MonoY2Lag  ,_N_s*_N_s*sizeof(scalar));						   mem_counter.addToGPUCounter(_N_s*_N_s*sizeof(scalar));
      cudaMalloc((void**) &_neighbors  ,_N_N*_N_E*sizeof(int));							   mem_counter.addToGPUCounter(_N_N*_N_E*sizeof(int));   
	    
      // Copy data to GPU
      cudaMemcpy(_Lag2MonoX,   tmpLag2MonoX, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice);
      cudaMemcpy(_MonoX2MonoY, tmpMonoX2MonoY, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice);
      cudaMemcpy(_MonoY2Lag,   tmpMonoY2Lag, N_s*N_s*sizeof(scalar), cudaMemcpyHostToDevice);
      cudaMemcpy(_neighbors,  m.getNeighbors(),_N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice);

      delete[] tmpLag2MonoX;
      delete[] tmpMonoX2MonoY;
      delete[] tmpMonoY2Lag;
#endif

      }
      break;
    default:
      _method = 0;
    }
    
    // For case specific stuff:
    switch (_method){
    case 1:{
#ifdef USE_CPU
      _A        = new scalar[_N_s*(_N_E+_N_ghosts)*N_F];    makeZero(_A,     _N_s*(_N_E+_N_ghosts)*N_F);
      _Alim     = new scalar[_N_s*(_N_E+_N_ghosts)*N_F];    makeZero(_Alim,  _N_s*(_N_E+_N_ghosts)*N_F);
      mem_counter.addToCPUCounter(2*_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar));
#elif USE_GPU
      cudaMalloc((void**) &_A          ,_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar));
      cudaMalloc((void**) &_Alim       ,_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar));
      mem_counter.addToGPUCounter(2*_N_s*(_N_E+_N_ghosts)*N_F*sizeof(scalar));
#endif
      }
      break;
    case 2:{

      // Some extra modal-nodal transforms
      blasGemm('N','N', _N_s, _N_s, _N_s, 1, _MonoY2Lag, _N_s, _MonoX2MonoY, _N_s, 0.0, _MonoX2Lag, _N_s);
      blasGemm('N','N', _N_s, _N_s, _N_s, 1, _MonoX2MonoY, _N_s, _Lag2MonoX, _N_s, 0.0, _Lag2MonoY, _N_s);

#ifdef USE_CPU
      _rho      = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhoMono      = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhoLim      = new scalar[_N_s*_N_E];
      _rhou     = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhouMono     = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhouLim     = new scalar[_N_s*_N_E];
      _rhov     = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhovMono     = new scalar[_N_s*(_N_E+_N_ghosts)]; _rhovLim     = new scalar[_N_s*_N_E];
      _pressure = new scalar[_N_s*(_N_E+_N_ghosts)]; _pressureMono = new scalar[_N_s*(_N_E+_N_ghosts)]; _pressureLim = new scalar[_N_s*_N_E];
      _K        = new scalar[_N_s*_N_E];  _KLim         = new scalar[_N_s*_N_E];
      _E        = new scalar[_N_s*_N_E];  _EMono        = new scalar[_N_s*_N_E];  _ELim        = new scalar[_N_s*_N_E]; makeZero(_ELim, _N_s*_N_E);
      _gamma    = new scalar[_N_s*(_N_E+_N_ghosts)]; _gammaMono    = new scalar[_N_s*(_N_E+_N_ghosts)]; _gammaLim    = new scalar[_N_s*_N_E];
      _rhoeLim  = new scalar[_N_s*_N_E];
      mem_counter.addToCPUCounter(10*_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      mem_counter.addToCPUCounter(12*_N_s*_N_E*sizeof(scalar));
#ifdef STIFFENED
      _beta    = new scalar[_N_s*(_N_E+_N_ghosts)]; _betaMono    = new scalar[_N_s*(_N_E+_N_ghosts)]; _betaLim    = new scalar[_N_s*_N_E];
      mem_counter.addToCPUCounter(2*_N_s*(_N_E+_N_ghosts)*sizeof(scalar)); mem_counter.addToCPUCounter(_N_s*_N_E*sizeof(scalar));
#endif
      // Mass fractions initializations
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) _Y(x) = new scalar[_N_s*(_N_E+_N_ghosts)]; _YMono(x) = new scalar[_N_s*(_N_E+_N_ghosts)]; _YLim(x) = new scalar[_N_s*_N_E];  mem_counter.addToCPUCounter(2*_N_s*(_N_E+_N_ghosts)*sizeof(scalar)); mem_counter.addToCPUCounter(_N_s*_N_E*sizeof(scalar));
#include "loop.h"

#elif USE_GPU
      cudaMalloc((void**) &_rho,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                      mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhoMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                  mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhoLim,_N_s*_N_E*sizeof(scalar));                               mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhou,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                     mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhouMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                 mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhouLim,_N_s*_N_E*sizeof(scalar));                              mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhov,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                     mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhovMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                 mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_rhovLim,_N_s*_N_E*sizeof(scalar));                              mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_pressure,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                 mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_pressureMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));             mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_pressureLim,_N_s*_N_E*sizeof(scalar));                          mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_K,_N_s*_N_E*sizeof(scalar));                                    mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_KLim,_N_s*_N_E*sizeof(scalar));                                 mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_E,_N_s*_N_E*sizeof(scalar));                                    mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_EMono,_N_s*_N_E*sizeof(scalar));                                mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_ELim,_N_s*_N_E*sizeof(scalar));       cudaMemset(_ELim, (scalar)0.0, _N_s*N_E*sizeof(scalar)); mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_gamma,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                    mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_gammaMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_gammaLim,_N_s*_N_E*sizeof(scalar));                             mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
      cudaMalloc((void**) &_rhoeLim,_N_s*_N_E*sizeof(scalar));                              mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
#ifdef STIFFENED
      cudaMalloc((void**) &_beta,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                     mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_betaMono,_N_s*(_N_E+_N_ghosts)*sizeof(scalar));                 mem_counter.addToGPUCounter(_N_s*(_N_E+_N_ghosts)*sizeof(scalar));
      cudaMalloc((void**) &_betaLim,_N_s*_N_E*sizeof(scalar));                              mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
#endif
      // Mass fractions initializations
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) cudaMalloc((void**) &_Y(x),_N_s*(_N_E+_N_ghosts)*sizeof(scalar)); \
      cudaMalloc((void**) &_YMono(x),_N_s*(_N_E+_N_ghosts)*sizeof(scalar));        \
      cudaMalloc((void**) &_YLim(x),_N_s*_N_E*sizeof(scalar));		\
      mem_counter.addToGPUCounter(2*_N_s*(_N_E+_N_ghosts)*sizeof(scalar)); mem_counter.addToGPUCounter(_N_s*_N_E*sizeof(scalar));
#include "loop.h"

#endif
    }
      break;
    case 3:
    case 4:{
      // Some extra modal-nodal transforms
      blasGemm('N','N', _N_s, _N_s, _N_s, 1, _MonoY2Lag, _N_s, _MonoX2MonoY, _N_s, 0.0, _MonoX2Lag, _N_s);
      blasGemm('N','N', _N_s, _N_s, _N_s, 1, _MonoX2MonoY, _N_s, _Lag2MonoX, _N_s, 0.0, _Lag2MonoY, _N_s);
#ifdef USE_CPU
      _Utmp     = new scalar[_N_s*_N_E*N_F];                        mem_counter.addToCPUCounter(_N_s*_N_E*N_F*sizeof(scalar));
#elif USE_GPU
      cudaMalloc((void**) &_Utmp,_N_s*_N_E*N_F*sizeof(scalar));     mem_counter.addToGPUCounter(_N_s*_N_E*N_F*sizeof(scalar));
#endif
    }
      break;
    }
  }// end 2D constructor for structured mesh
  
/*   /\*!\brief Constructor for 2D limiting for unstructured mesh*\/ */
/*  Limiting(int method, int N_s, int N_E, int N_G, int N_N, int L, int order, int L2Msize1, int L2Msize2, simpleMesh &m, fullMatrix<scalar> Lag2Mono, fullMatrix<scalar> Mono2Lag, fullMatrix<scalar> XYZCen, scalar* powersXYZG, scalar* weight, scalar refArea, int* TaylorDxIdx, int* TaylorDyIdx) */
/*    : _method(method), _D(D),_N_s(N_s), _N_E(N_E), _N_G(N_G), _N_N(N_N), _L(L), _order(order), _L2Msize1(L2Msize1), _L2Msize2(L2Msize2), _refArea(refArea){ */

/*     common_ctor(); */
    
/*     switch (_method){ */
/*     case 1: */
/*     case 2: */
/*     case 3:{ */
/* #ifdef USE_CPU */
/*       // Allocate */
/*       _Lag2Mono   = new scalar[_L2Msize1*_L2Msize2*_N_E]; */
/*       _Mono2Lag   = new scalar[_L2Msize2*_L2Msize1*_N_E]; */
/*       _XYZCen     = new scalar[_N_E*_D]; */
/*       _powersXYZG = new scalar[_L2Msize1*_N_G*(_N_N+1)*_N_E]; */
/*       _neighbors  = new int[_N_N*_N_E]; */
/*       _weight     = new scalar[_N_G]; */
/*       _TaylorDxIdx= new int[_L]; */
/*       _TaylorDyIdx= new int[_L]; */
/*       _A          = new scalar[_L2Msize1*_N_E*N_F]; makeZero(_A,    _L2Msize1*_N_E*N_F); */
/*       _Alim       = new scalar[_L2Msize1*_N_E*N_F]; makeZero(_Alim, _L2Msize1*_N_E*N_F); */

/*       // Copy the data to these new pointers */
/*       for(int e = 0; e < _N_E; e++){ */
/* 	for(int alpha = 0; alpha < D; alpha++){ _XYZCen[e*D+alpha] = XYZCen(e,alpha);} */
/*       	for(int i = 0; i < _L2Msize1; i++){ */
/*       	  for(int j = 0; j < _L2Msize2; j++){ */
/*       	    _Lag2Mono[(e*_L2Msize1+i)*_L2Msize2+j] = Lag2Mono(e,i*_L2Msize2+j); */
/*       	    _Mono2Lag[(e*_L2Msize2+j)*_L2Msize1+i] = Mono2Lag(e,j*_L2Msize1+i);}}} */
/*       memcpy(_powersXYZG,  powersXYZG,  _L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar)); */
/*       memcpy(_neighbors,   m.getNeighbors(),   _N_N*_N_E*sizeof(int)); */
/*       memcpy(_weight,      weight,      _N_G*sizeof(scalar)); */
/*       memcpy(_TaylorDxIdx, TaylorDxIdx, _L*sizeof(int)); */
/*       memcpy(_TaylorDyIdx, TaylorDyIdx, _L*sizeof(int)); */
      
/* #elif USE_GPU */
/*       // tmp host pointers to copy data to gpu */
/*       scalar* tmpLag2Mono   = new scalar[_L2Msize1*_L2Msize2*_N_E]; */
/*       scalar* tmpMono2Lag   = new scalar[_L2Msize2*_L2Msize1*_N_E]; */
/*       scalar* tmpXYZCen     = new scalar[_N_E*_D]; */
/*       for(int e = 0; e < _N_E; e++){ */
/* 	for(int alpha = 0; alpha < D; alpha++){ tmpXYZCen[e*D+alpha] = XYZCen(e,alpha);} */
/* 	for(int i = 0; i < _L2Msize1; i++){ */
/* 	  for(int j = 0; j < _L2Msize2; j++){ */
/*       	    tmpLag2Mono[(e*_L2Msize1+i)*_L2Msize2+j] = Lag2Mono(e,i*_L2Msize2+j); */
/*       	    tmpMono2Lag[(e*_L2Msize2+j)*_L2Msize1+i] = Mono2Lag(e,j*_L2Msize1+i);}}} */

/*       // Allocate on GPU */
/*       cudaMalloc((void**) &_Lag2Mono,_L2Msize1*_L2Msize2*_N_E*sizeof(scalar)); */
/*       cudaMalloc((void**) &_Mono2Lag,_L2Msize2*_L2Msize1*_N_E*sizeof(scalar)); */
/*       cudaMalloc((void**) &_XYZCen,_D*_N_E*sizeof(scalar)); */
/*       cudaMalloc((void**) &_powersXYZG,_L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar)); */
/*       cudaMalloc((void**) &_neighbors,_N_N*_N_E*sizeof(int)); */
/*       cudaMalloc((void**) &_weight,_N_G*sizeof(scalar)); */
/*       cudaMalloc((void**) &_TaylorDxIdx,_L*sizeof(int)); */
/*       cudaMalloc((void**) &_TaylorDyIdx,_L*sizeof(int)); */
/*       cudaMalloc((void**) &_A,_L2Msize1*_N_E*N_F*sizeof(scalar)); */
/*       cudaMalloc((void**) &_Alim,_L2Msize1*_N_E*N_F*sizeof(scalar)); */

/*       // Copy data to GPU */
/*       cudaMemcpy(_Lag2Mono,   tmpLag2Mono, _L2Msize1*_L2Msize2*_N_E*sizeof(scalar), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_Mono2Lag,   tmpMono2Lag, _L2Msize2*_L2Msize1*_N_E*sizeof(scalar), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_XYZCen,     tmpXYZCen,   _D*_N_E*sizeof(scalar), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_powersXYZG, powersXYZG,  _L2Msize1*_N_G*(_N_N+1)*_N_E*sizeof(scalar), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_neighbors,  m.getNeighbors(),   _N_N*_N_E*sizeof(int), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_weight,     weight,      _N_G*sizeof(scalar), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_TaylorDxIdx,TaylorDxIdx, _L*sizeof(int), cudaMemcpyHostToDevice); */
/*       cudaMemcpy(_TaylorDyIdx,TaylorDyIdx, _L*sizeof(int), cudaMemcpyHostToDevice); */
/*       cudaMemset(_A,    (scalar)0.0, _L2Msize1*_N_E*N_F*sizeof(scalar)); */
/*       cudaMemset(_Alim, (scalar)0.0, _L2Msize1*_N_E*N_F*sizeof(scalar)); */
      
/*       delete[] tmpLag2Mono; */
/*       delete[] tmpMono2Lag; */
/*       delete[] tmpXYZCen; */
/* #endif */
/*       } */
/*       break; */
/*     default: */
/*       _method = 0; */
/*     } */
    
/*     // For case specific stuff: */
/*   } // end 2D constructor */
  
  /*! \brief Destructor */
  ~Limiting(){
    if(_Lag2Mono)     del(_Lag2Mono);
    if(_Mono2Lag)     del(_Mono2Lag);
    if(_Lag2MonoX)    del(_Lag2MonoX);
    if(_Lag2MonoY)    del(_Lag2MonoY);
    if(_MonoX2MonoY)  del(_MonoX2MonoY);
    if(_MonoX2Lag)    del(_MonoX2Lag);
    if(_MonoY2Lag)    del(_MonoY2Lag);
    if(_A)            del(_A);
    if(_Alim)         del(_Alim);
    if(_Utmp)         del(_Utmp);
    if(_rho){         del(_rho); del(_rhoMono); del(_rhoLim);}
    if(_rhou){        del(_rhou); del(_rhouMono); del(_rhouLim);}
    if(_rhov){        del(_rhov); del(_rhovMono); del(_rhovLim);}
    if(_pressure){    del(_pressure); del(_pressureMono); del(_pressureLim);}
    if(_K){           del(_K); del(_KLim);}
    if(_E){           del(_E); del(_EMono); del(_ELim);}
    if(_gamma){       del(_gamma); del(_gammaMono); del(_gammaLim);}
    if(_beta){        del(_beta); del(_betaMono); del(_betaLim);}
    if(_rhoeLim){     del(_rhoeLim);}
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) if(_Y(x)) del(_Y(x)); del(_YMono(x)); del(_YLim(x));
#include "loop.h"
    if(_XYZCen)       del(_XYZCen);
    if(_powersXYZG)   del(_powersXYZG);
    if(_neighbors)    del(_neighbors);
    if(_TaylorDxIdx)  del(_TaylorDxIdx);
    if(_TaylorDyIdx)  del(_TaylorDyIdx);
  }	      
  
  int getLimitingMethod() const {/*!Return limiting method*/return _method;}
  void HRlimiting(COMMUNICATOR &communicator, scalar* U);
  void M2limiting(COMMUNICATOR &communicator, scalar* U);
  void HRIlimiting(COMMUNICATOR &communicator, SENSOR &sensor, scalar* U);
  void M2Ilimiting(COMMUNICATOR &communicator, SENSOR &sensor, scalar* U);
};
#endif
