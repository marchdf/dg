#ifndef LIMITING_H
#define LIMITING_H

#include <misc.h>
#include <constants.h>

class Limiting
{
 private:
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
  }
  
  // destructor
  ~Limiting(){
    if(_method!=0){
      if(_Lag2Mono) del(_Lag2Mono);
      if(_Mono2Lag) del(_Mono2Lag);
      if(_V1D)      del(_V1D);
      if(_weight)   del(_weight);
      if(_A)        del(_A);
      if(_Alim)     del(_Alim);
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
    scalar* pressure     = new scalar[_N_s*_N_E];
    scalar* pressureMono = new scalar[_N_s*_N_E];
    scalar* pressureLim  = new scalar[_N_s*_N_E];
    for(int e = 0; e < _N_E; e++){
      for(int i = 0; i < _N_s; i++){
	scalar rho  = U[(e*_N_F+0)*_N_s+i];
	scalar rhou = U[(e*_N_F+1)*_N_s+i];
	scalar E    = U[(e*_N_F+2)*_N_s+i];
#ifdef MULTIFLUID
#ifdef GAMCONS
	scalar gamma=1.0+rho/U[(e*_N_F+3)*_N_s+i];
#elif GAMNCON
	scalar gamma=1.0+1.0/U[(e*_N_F+3)*_N_s+i];
#endif
#elif PASSIVE
	scalar gamma = constants::GLOBAL_GAMMA;
#endif
	pressure[e*_N_s+i] = (gamma-1)*(E - 0.5*rhou*rhou/rho);
      }
    }

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, pressure, _N_s, 0.0, pressureMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, pressureMono, pressureLim);

    // Go from lagrange to monomial representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    // Limit the solution according to Liu
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);

    // My modification:
    if (_N_s == 2){
      for(int e = 0; e < _N_E; e++){
#ifdef MULTIFLUID
	  scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	  scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	  scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	  scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	  scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	  scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	  scalar d0 = _Alim[(e*_N_F+3)*_N_s+0];
	  scalar d1 = _Alim[(e*_N_F+3)*_N_s+1];
#ifdef GAMCONS
	  scalar D1 = 0.5*((d0+d1)/(a0+a1)-(d0-d1)/(a0-a1));
	  scalar D0 = 0.5*((d0-d1)/(a0-a1)+(d0+d1)/(a0+a1));
	  d0 = D0; d1 = D1;
#endif
	  scalar g1lim = d1/d0*(g0-0.25*((b0+b1)*(b0+b1)/(a0+a1) + (b0-b1)*(b0-b1)/(a0-a1))) + 0.25*((b0+b1)*(b0+b1)/(a0+a1) - (b0-b1)*(b0-b1)/(a0-a1)) + pressureLim[e*_N_s+1]*d0;
	  scalar g0lim = g0;
	  _Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	  _Alim[(e*_N_F+2)*_N_s+0] = g0lim;
#elif PASSIVE
	  scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	  scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	  scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	  scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	  scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	  scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	  scalar gamma = constants::GLOBAL_GAMMA;
	  scalar g1lim = 0.25*((b0+b1)*(b0+b1)/(a0+a1) - (b0-b1)*(b0-b1)/(a0-a1)) + pressureLim[e*_N_s+1]*gamma;
	  scalar g0lim = g0;
	  _Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	  _Alim[(e*_N_F+2)*_N_s+0] = g0lim;
#endif
      }
    }

    if (_N_s == 3){
      for(int e = 0; e < _N_E; e++){
#ifdef MULTIFLUID
	scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	scalar a2 = _Alim[(e*_N_F+0)*_N_s+2];
	scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	scalar b2 = _Alim[(e*_N_F+1)*_N_s+2];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	scalar g2 = _A[(e*_N_F+2)*_N_s+2];
	scalar d0 = _Alim[(e*_N_F+3)*_N_s+0];
	scalar d1 = _Alim[(e*_N_F+3)*_N_s+1];
	scalar d2 = _Alim[(e*_N_F+3)*_N_s+2];
#ifdef GAMCONS
	scalar D0 = d0/a0;
	scalar D1 = 0.5*((d0+d1+0.5*d2)/(a0+a1+0.5*a2) - (d0-d1+0.5*d2)/(a0-a1+0.5*a2));
	scalar D2 = (d0+d1+0.5*d2)/(a0+a1+0.5*a2) + (d0-d1+0.5*d2)/(a0-a1+0.5*a2) - 2*D0;
	d0 = D0; d1 = D1; d2 = D2;
#endif
	scalar g2lim = (1/(1+1.0/6.0*d2/d0))
	  *(d2/d0*(g0+1.0/6.0*g2 -0.5*b0*b0/a0)
	    + 0.5*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) + (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2) - 2*b0*b0/a0)
	    + pressureLim[e*_N_s+2]*d0);
	scalar g1lim = d1/d0*(g0+1.0/6.0*(g2-g2lim)-0.5*b0*b0/a0)
	  + 0.25*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) - (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2))
	  + pressureLim[e*_N_s+1]*d0;
	scalar g0lim = g0 + 1.0/6.0*(g2-g2lim);
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+2] = g2lim;
#elif PASSIVE
	scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	scalar a2 = _Alim[(e*_N_F+0)*_N_s+2];
	scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	scalar b2 = _Alim[(e*_N_F+1)*_N_s+2];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	scalar g2 = _A[(e*_N_F+2)*_N_s+2];
	scalar gamma = constants::GLOBAL_GAMMA;
	scalar g2lim = 0.5*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) + (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2) - 2*b0*b0/a0) + pressureLim[e*_N_s+2]*gamma;
	scalar g1lim = 0.25*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) - (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2))
	  + pressureLim[e*_N_s+1]*gamma;
	scalar g0lim = g0 + 1.0/6.0*(g2-g2lim);
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+2] = g2lim;
#endif
      }
    }
    
    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

    delete[] pressure;
    delete[] pressureMono;
    delete[] pressureLim;
	    
  }

  void M2limiting(scalar* U){

    // Get the pressure field
    scalar* pressure     = new scalar[_N_s*_N_E];
    scalar* pressureMono = new scalar[_N_s*_N_E];
    scalar* pressureLim  = new scalar[_N_s*_N_E];
    scalar* u            = new scalar[_N_s*_N_E];
    scalar* uMono        = new scalar[_N_s*_N_E];
    scalar* uLim         = new scalar[_N_s*_N_E];
    for(int e = 0; e < _N_E; e++){
      for(int i = 0; i < _N_s; i++){
	scalar rho  = U[(e*_N_F+0)*_N_s+i];
	scalar rhou = U[(e*_N_F+1)*_N_s+i];
	scalar E    = U[(e*_N_F+2)*_N_s+i];
	u[e*_N_s+i] = rhou/rho;
#ifdef MULTIFLUID
#ifdef GAMCONS
	scalar gamma=1.0+rho/U[(e*_N_F+3)*_N_s+i];
#elif GAMNCON
	scalar gamma=1.0+1.0/U[(e*_N_F+3)*_N_s+i];
#endif
#elif PASSIVE
	scalar gamma = constants::GLOBAL_GAMMA;
#endif
	pressure[e*_N_s+i] = (gamma-1.0)*(E - 0.5*rhou*rhou/rho);
      }
    }

    // Limit pressure
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, pressure, _N_s, 0.0, pressureMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, pressureMono, pressureLim);

    // Limit velocity
    blasGemm('N','N', _N_s, _N_E, _N_s, 1, _Lag2Mono, _N_s, u, _N_s, 0.0, uMono, _N_s);
    Lcpu_hrl(_N_s, _N_E, 1, _N_G, _boundaryMap, _weight, _V1D, uMono, uLim);
    
    // Limit the other variables
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Lag2Mono, _N_s, U, _N_s, 0.0, _A, _N_s);
    Lcpu_hrl(_N_s, _N_E, _N_F, _N_G, _boundaryMap, _weight, _V1D, _A, _Alim);

    // M2 modification:
    if (_N_s == 2){
      for(int e = 0; e < _N_E; e++){
#ifdef MULTIFLUID
	scalar a0lim = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1lim = _Alim[(e*_N_F+0)*_N_s+1];
	scalar d0lim = _Alim[(e*_N_F+3)*_N_s+0];
	scalar d1lim = _Alim[(e*_N_F+3)*_N_s+1];
	scalar b0 = _A[(e*_N_F+1)*_N_s+0];
	scalar b1 = _A[(e*_N_F+1)*_N_s+1];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];

	// Momentum (to avoid u oscillations)
	scalar b0lim = b0;
	scalar b1lim = b0lim/a0lim*a1lim + uLim[e*_N_s+1]*a0lim;
	_Alim[(e*_N_F+1)*_N_s+1] = b1lim;
	_Alim[(e*_N_F+1)*_N_s+0] = b0lim;
	  
	// Energy (like MYL)
#ifdef GAMCONS
	scalar D1lim = 0.5*((d0lim+d1lim)/(a0lim+a1lim)-(d0lim-d1lim)/(a0lim-a1lim));
	scalar D0lim = 0.5*((d0lim-d1lim)/(a0lim-a1lim)+(d0lim+d1lim)/(a0lim+a1lim));
	d0lim = D0lim; d1lim = D1lim;
#endif
	scalar g1lim = d1lim/d0lim*(g0-0.25*((b0lim+b1lim)*(b0lim+b1lim)/(a0lim+a1lim) + (b0lim-b1lim)*(b0lim-b1lim)/(a0lim-a1lim))) + 0.25*((b0lim+b1lim)*(b0lim+b1lim)/(a0lim+a1lim) - (b0lim-b1lim)*(b0lim-b1lim)/(a0lim-a1lim)) + pressureLim[e*_N_s+1]*d0lim;
	scalar g0lim = g0;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
#elif PASSIVE
	scalar a0lim = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1lim = _Alim[(e*_N_F+0)*_N_s+1];
	scalar b0 = _A[(e*_N_F+1)*_N_s+0];
	scalar b1 = _A[(e*_N_F+1)*_N_s+1];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];

	// Momentum (to avoid u oscillations)
	scalar b0lim = b0;
	scalar b1lim = b0lim/a0lim*a1lim + uLim[e*_N_s+1]*a0lim;
	_Alim[(e*_N_F+1)*_N_s+1] = b1lim;
	_Alim[(e*_N_F+1)*_N_s+0] = b0lim;

	// Energy (like MYL)
	scalar gamma = constants::GLOBAL_GAMMA;
	scalar g1lim = 0.25*((b0lim+b1lim)*(b0lim+b1lim)/(a0lim+a1lim) - (b0lim-b1lim)*(b0lim-b1lim)/(a0lim-a1lim)) + pressureLim[e*_N_s+1]*gamma;
	scalar g0lim = g0;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
#endif
      }
    }

    if (_N_s == 3){
      for(int e = 0; e < _N_E; e++){
#ifdef MULTIFLUID
	scalar a0lim = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1lim = _Alim[(e*_N_F+0)*_N_s+1];
	scalar a2lim = _Alim[(e*_N_F+0)*_N_s+2];
	scalar d0lim = _Alim[(e*_N_F+3)*_N_s+0];
	scalar d1lim = _Alim[(e*_N_F+3)*_N_s+1];
	scalar d2lim = _Alim[(e*_N_F+3)*_N_s+2];
	scalar b0 = _A[(e*_N_F+1)*_N_s+0];
	scalar b1 = _A[(e*_N_F+1)*_N_s+1];
	scalar b2 = _A[(e*_N_F+1)*_N_s+2];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	scalar g2 = _A[(e*_N_F+2)*_N_s+2];

	// Momentum (to avoid u oscillations)
	scalar b2lim = (1/(1+1.0/6.0*a2lim/a0lim))*((b0+1.0/6.0*b2)/a0lim*a2lim + uLim[e*_N_s+2]*a0lim);
	scalar b0lim = b0 + 1.0/6.0*(b2-b2lim);
	scalar b1lim = b0lim/a0lim*a1lim + uLim[e*_N_s+1]*a0lim;
	_Alim[(e*_N_F+1)*_N_s+0] = b0lim;
	_Alim[(e*_N_F+1)*_N_s+1] = b1lim;
	_Alim[(e*_N_F+1)*_N_s+2] = b2lim;

	// Energy (like MYL)
#ifdef GAMCONS
	scalar D0lim = d0lim/a0lim;
	scalar D1lim = 0.5*((d0lim+d1lim+0.5*d2lim)/(a0lim+a1lim+0.5*a2lim) - (d0lim-d1lim+0.5*d2lim)/(a0lim-a1lim+0.5*a2lim));
	scalar D2lim = (d0lim+d1lim+0.5*d2lim)/(a0lim+a1lim+0.5*a2lim) + (d0lim-d1lim+0.5*d2lim)/(a0lim-a1lim+0.5*a2lim) - 2*D0lim;
	d0lim = D0lim; d1lim = D1lim; d2lim = D2lim;
#endif
	scalar g2lim = (1/(1+1.0/6.0*d2lim/d0lim))
	  *(d2lim/d0lim*(g0+1.0/6.0*g2 -0.5*b0lim*b0lim/a0lim)
	    + 0.5*((b0lim+b1lim+0.5*b2lim)*(b0lim+b1lim+0.5*b2lim)/(a0lim+a1lim+0.5*a2lim) + (b0lim-b1lim+0.5*b2lim)*(b0lim-b1lim+0.5*b2lim)/(a0lim-a1lim+0.5*a2lim) - 2*b0lim*b0lim/a0lim)
	    + pressureLim[e*_N_s+2]*d0lim);
	scalar g1lim = d1lim/d0lim*(g0+1.0/6.0*(g2-g2lim)-0.5*b0lim*b0lim/a0lim)
	  + 0.25*((b0lim+b1lim+0.5*b2lim)*(b0lim+b1lim+0.5*b2lim)/(a0lim+a1lim+0.5*a2lim) - (b0lim-b1lim+0.5*b2lim)*(b0lim-b1lim+0.5*b2lim)/(a0lim-a1lim+0.5*a2lim))
	  + pressureLim[e*_N_s+1]*d0lim;
	scalar g0lim = g0 + 1.0/6.0*(g2-g2lim);
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+2] = g2lim;
#elif PASSIVE
	scalar a0lim = _Alim[(e*_N_F+0)*_N_s+0];
	scalar a1lim = _Alim[(e*_N_F+0)*_N_s+1];
	scalar a2lim = _Alim[(e*_N_F+0)*_N_s+2];
	scalar b0 = _A[(e*_N_F+1)*_N_s+0];
	scalar b1 = _A[(e*_N_F+1)*_N_s+1];
	scalar b2 = _A[(e*_N_F+1)*_N_s+2];
	scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	scalar g2 = _A[(e*_N_F+2)*_N_s+2];

	// Momentum (to avoid u oscillations)
	scalar b2lim = (1/(1+1.0/6.0*a2lim/a0lim))*((b0+1.0/6.0*b2)/a0lim*a2lim + uLim[e*_N_s+2]*a0lim);
	scalar b0lim = b0 + 1.0/6.0*(b2-b2lim);
	scalar b1lim = b0lim/a0lim*a1lim + uLim[e*_N_s+1]*a0lim;
	_Alim[(e*_N_F+1)*_N_s+0] = b0lim;
	_Alim[(e*_N_F+1)*_N_s+1] = b1lim;
	_Alim[(e*_N_F+1)*_N_s+2] = b2lim;

	// Energy (like MYL)
	scalar gamma = constants::GLOBAL_GAMMA;
	scalar g2lim = 0.5*((b0lim+b1lim+0.5*b2lim)*(b0lim+b1lim+0.5*b2lim)/(a0lim+a1lim+0.5*a2lim) + (b0lim-b1lim+0.5*b2lim)*(b0lim-b1lim+0.5*b2lim)/(a0lim-a1lim+0.5*a2lim) - 2*b0lim*b0lim/a0lim)+ pressureLim[e*_N_s+2]*gamma;
	scalar g1lim =  0.25*((b0lim+b1lim+0.5*b2lim)*(b0lim+b1lim+0.5*b2lim)/(a0lim+a1lim+0.5*a2lim) - (b0lim-b1lim+0.5*b2lim)*(b0lim-b1lim+0.5*b2lim)/(a0lim-a1lim+0.5*a2lim)) + pressureLim[e*_N_s+1]*gamma;
	scalar g0lim = g0 + 1.0/6.0*(g2-g2lim);
	_Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	_Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	_Alim[(e*_N_F+2)*_N_s+2] = g2lim;
#endif
      }
    }


    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

    delete[] pressure;
    delete[] pressureMono;
    delete[] pressureLim;
    delete[] u;
    delete[] uMono;
    delete[] uLim;
	    
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
