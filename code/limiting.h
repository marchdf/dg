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
  bool    _multifluid;
  bool    _passive;
  int     _model;
  scalar  _gamma0;
  
 public:
  // constructor
 Limiting(int method) : _method(method){}
  
 Limiting(int method, int N_s, int N_E, int N_F, int N_G, int boundaryMap, bool multifluid, bool passive, int model, scalar gamma0, fullMatrix<scalar> &Lag2Mono, fullMatrix<scalar> &Mono2Lag, fullMatrix<scalar> &V1D, scalar* weight)
   : _method(method), _N_s(N_s), _N_E(N_E), _N_F(N_F), _N_G(N_G), _boundaryMap(boundaryMap), _multifluid(multifluid), _passive(passive), _model(model), _gamma0(gamma0){
    switch (_method){
    case 1:
    case 2:
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
	if(_multifluid){
	  scalar gamma=0;
	  if      (_model==0) gamma=1.0+rho/U[(e*_N_F+3)*_N_s+i];
	  else if (_model==1) gamma=1.0+1.0/U[(e*_N_F+3)*_N_s+i];
	  pressure[e*_N_s+i] = (gamma-1.0)*(E - 0.5*rhou*rhou/rho);
	}
	else if(_passive){
	  pressure[e*_N_s+i] = (_gamma0-1)*(E - 0.5*rhou*rhou/rho);
	}
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
	if (_multifluid){
	  scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	  scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	  scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	  scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	  scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	  scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	  scalar d0 = _Alim[(e*_N_F+3)*_N_s+0];
	  scalar d1 = _Alim[(e*_N_F+3)*_N_s+1];
	  if (_model==0){
	    scalar D1 = 0.5*((d0+d1)/(a0+a1)-(d0-d1)/(a0-a1));
	    scalar D0 = 0.5*((d0-d1)/(a0-a1)+(d0+d1)/(a0+a1));
	    d0 = D0; d1 = D1;
	  }
	  scalar g1lim = d1/d0*(g0-0.25*((b0+b1)*(b0+b1)/(a0+a1) + (b0-b1)*(b0-b1)/(a0-a1))) + 0.25*((b0+b1)*(b0+b1)/(a0+a1) - (b0-b1)*(b0-b1)/(a0-a1)) + pressureLim[e*_N_s+1]*d0;
	  scalar g0lim = g0;
	  _Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	  _Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	}
	else if (_passive){
	  scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	  scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	  scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	  scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	  scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	  scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	  scalar g1lim = 0.25*((b0+b1)*(b0+b1)/(a0+a1) - (b0-b1)*(b0-b1)/(a0-a1)) + pressureLim[e*_N_s+1]*_gamma0;
	  scalar g0lim = g0;
	  _Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	  _Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	}
      }
    }

    if (_N_s == 3){
      for(int e = 0; e < _N_E; e++){
	if (_multifluid){
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
	  if (_model==0){
	    scalar D0 = d0/a0;
	    scalar D1 = 0.5*((d0+d1+0.5*d2)/(a0+a1+0.5*a2) - (d0-d1+0.5*d2)/(a0-a1+0.5*a2));
	    scalar D2 = (d0+d1+0.5*d2)/(a0+a1+0.5*a2) + (d0-d1+0.5*d2)/(a0-a1+0.5*a2) - 2*D0;
	    d0 = D0; d1 = D1; d2 = D2;
	  }
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
	}
	else if (_passive){
	  scalar a0 = _Alim[(e*_N_F+0)*_N_s+0];
	  scalar a1 = _Alim[(e*_N_F+0)*_N_s+1];
	  scalar a2 = _Alim[(e*_N_F+0)*_N_s+2];
	  scalar b0 = _Alim[(e*_N_F+1)*_N_s+0];
	  scalar b1 = _Alim[(e*_N_F+1)*_N_s+1];
	  scalar b2 = _Alim[(e*_N_F+1)*_N_s+2];
	  scalar g0 = _A[(e*_N_F+2)*_N_s+0];
	  scalar g1 = _A[(e*_N_F+2)*_N_s+1];
	  scalar g2 = _A[(e*_N_F+2)*_N_s+2];
	  scalar g2lim = 0.5*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) + (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2) - 2*b0*b0/a0) + pressureLim[e*_N_s+2]*_gamma0;
	  scalar g1lim = 0.25*((b0+b1+0.5*b2)*(b0+b1+0.5*b2)/(a0+a1+0.5*a2) - (b0-b1+0.5*b2)*(b0-b1+0.5*b2)/(a0-a1+0.5*a2))
	    + pressureLim[e*_N_s+1]*_gamma0;
	  scalar g0lim = g0 + 1.0/6.0*(g2-g2lim);
	  _Alim[(e*_N_F+2)*_N_s+0] = g0lim;
	  _Alim[(e*_N_F+2)*_N_s+1] = g1lim;
	  _Alim[(e*_N_F+2)*_N_s+2] = g2lim;
	}
      }
    }


    // Go back to lagrange representation
    blasGemm('N','N', _N_s, _N_E*_N_F, _N_s, 1, _Mono2Lag, _N_s, _Alim, _N_s, 0.0, U, _N_s);

    delete[] pressure;
    delete[] pressureMono;
    delete[] pressureLim;
	    
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
