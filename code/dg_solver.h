//
// DG solver class
//
#ifndef DG_SOLVER_H
#define DG_SOLVER_H

#include <stdio.h>
#include <scalar_def.h>
#include <macros.h>

class DG_SOLVER
{
 private:
  int _D;
  int _N_F;
  int _N_E;
  int _N_s;
  int _N_G;
  int _M_T;
  int _M_s;
  int _M_G;
  int _boundaryMap;
  int _flux;
  int _model;
  scalar _gamma0;
  int _blas;
  bool _multifluid;
  bool _passive;
  scalar* h_phi     ; 
  scalar* h_phi_w   ; 
  scalar* h_dphi    ; 
  scalar* h_dphi_w  ;
  scalar* h_J;
  scalar* h_invJac  ; 
  scalar* h_UF      ; 
  scalar* h_Uinteg  ; 
  scalar* h_dUinteg ; 
  scalar* h_UintegF ; 
  scalar* h_s       ; 
  scalar* h_sJ      ; 
  scalar* h_S       ; 
  scalar* h_f       ; 
  scalar* h_fJ      ; 
  scalar* h_F       ; 
  scalar* h_q       ; 
  scalar* h_Q       ; 

  
#ifdef USE_GPU
  scalar* d_phi, *d_phi_w, *d_dphi, *d_dphi_w, *d_invJac;
  scalar* d_UF;
  scalar* d_Uinteg, *d_dUinteg, *d_UintegF;
#endif

  
 public:
  // constructor
 DG_SOLVER(int D, int N_F, int N_E, int N_s, int N_G, int M_T, int M_s, int M_G,
	   scalar* phi, scalar* dphi, scalar* phi_w, scalar* dphi_w, scalar* J, scalar* invJac,
	   int boundaryMap, int flux, int model, scalar gamma0, int blas, bool multifluid, bool passive) :
  _D(D), _N_F(N_F), _N_E(N_E), _N_s(N_s), _N_G(N_G), _M_T(M_T), _M_s(M_s), _M_G(M_G),
    _boundaryMap(boundaryMap), _flux(flux), _model(model), _gamma0(gamma0), _blas(blas), _multifluid(multifluid), _passive(passive) {

#ifdef USE_CPU

    h_phi     = new scalar[N_G*N_s];          makeZero(h_phi,N_G*N_s);
    h_phi_w   = new scalar[N_G*N_s];          makeZero(h_phi_w,N_G*N_s);          
    h_dphi    = new scalar[D*N_G*N_s];	      makeZero(h_dphi,D*N_G*N_s);	 
    h_dphi_w  = new scalar[D*N_G*N_s];	      makeZero(h_dphi_w,D*N_G*N_s);
    h_J       = new scalar[N_E];              makeZero(h_J,N_E);                                 // not same as J!!
    h_invJac  = new scalar[N_G*D*N_E*D];      makeZero(h_invJac,N_G*D*N_E*D);                    // not same as invJac!!
    h_UF      = new scalar[2*N_F*M_s*M_T];    makeZero(h_UF,2*N_F*M_s*M_T); 
    h_Uinteg  = new scalar[N_F*N_G*N_E];      makeZero(h_Uinteg,N_F*N_G*N_E);	 
    h_dUinteg = new scalar[D*N_G*N_E*N_F];    makeZero(h_dUinteg,D*N_G*N_E*N_F); 
    h_UintegF = new scalar[2*N_F*M_G*M_T];    makeZero(h_UintegF,2*N_F*M_G*M_T); 
    h_s       = new scalar[N_G*N_E*N_F];      makeZero(h_s,N_G*N_E*N_F);	 
    h_sJ      = new scalar[N_G*N_E*N_F];      makeZero(h_sJ,N_G*N_E*N_F);	 
    h_S       = new scalar[N_s*N_E*N_F];      makeZero(h_S,N_s*N_E*N_F);	 
    h_f       = new scalar[D*N_F*N_G*N_E];    makeZero(h_f,D*N_F*N_G*N_E); 
    h_fJ      = new scalar[D*N_G*N_E*N_F];    makeZero(h_fJ,D*N_G*N_E*N_F); 
    h_F       = new scalar[N_s*N_E*N_F];      makeZero(h_F,N_s*N_E*N_F);	 
    h_q       = new scalar[M_G*M_T*N_F*2];    makeZero(h_q,M_G*M_T*N_F*2); 
    h_Q       = new scalar[N_s*N_E*N_F];      makeZero(h_Q,N_s*N_E*N_F);   
    
    memcpy(h_phi,    phi,    N_G*N_s*sizeof(scalar));
    memcpy(h_phi_w,  phi_w,  N_G*N_s*sizeof(scalar));
    memcpy(h_dphi,   dphi,   D*N_G*N_s*sizeof(scalar));
    memcpy(h_dphi_w, dphi_w, D*N_G*N_s*sizeof(scalar));
    memcpy(h_J, J, N_E*sizeof(scalar));
    memcpy(h_invJac, invJac, N_G*D*N_E*D*sizeof(scalar));
	   
#elif USE_GPU
    // Allocate space on the GPU
    scalar* d_phi, *d_phi_w, *d_dphi, *d_dphi_w;
    scalar* d_J, *d_invJac;
    scalar* d_UF;
    scalar* d_Uinteg, *d_dUinteg, *d_UintegF;
    scalar* d_s, *d_f, *d_q; 
    scalar* d_sJ, *d_fJ; 
    scalar* d_S, *d_F, *d_Q;
  
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_phi,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_phi_w,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dphi,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dphi_w,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_J,N_E*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_invJac,N_G*D*N_E*D*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_UF,M_s*M_T*N_F*2*sizeof(scalar)));

    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Uinteg,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_dUinteg,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_UintegF,M_G*M_T*N_F*2*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_s,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_sJ,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_S,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_f,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_fJ,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_F,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_q,M_G*M_T*N_F*2*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_Q,N_s*N_E*N_F*sizeof(scalar)));


    // Send the stuff to the device
    CUDA_SAFE_CALL(cudaMemcpy(d_phi, h_phi, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_phi_w, h_phi_w, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_dphi, h_dphi, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_dphi_w, h_dphi_w, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_J, h_J, N_E*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_invJac, h_invJac, N_G*D*N_E*D*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemset(d_Q, (scalar)0.0, N_E*N_F*N_s*sizeof(scalar)));


#endif


  };

  // destructor
  ~DG_SOLVER(){
#ifdef USE_CPU
    delete[] h_phi;
    delete[] h_phi_w;
    delete[] h_dphi;
    delete[] h_dphi_w;
    delete[] h_J;
    delete[] h_invJac;
    delete[] h_UF;
    delete[] h_Uinteg;
    delete[] h_dUinteg;
    delete[] h_UintegF;
    delete[] h_s;
    delete[] h_sJ;
    delete[] h_S;
    delete[] h_f;
    delete[] h_fJ;
    delete[] h_F;
    delete[] h_q;
    delete[] h_Q;
#elif USE_GPU
    CUDA_SAFE_CALL(cudaFree(d_phi));
    CUDA_SAFE_CALL(cudaFree(d_phi_w));
    CUDA_SAFE_CALL(cudaFree(d_dphi));
    CUDA_SAFE_CALL(cudaFree(d_dphi_w));
    CUDA_SAFE_CALL(cudaFree(d_J));
    CUDA_SAFE_CALL(cudaFree(d_invJac));
    CUDA_SAFE_CALL(cudaFree(d_UF));
    CUDA_SAFE_CALL(cudaFree(d_Uinteg));
    CUDA_SAFE_CALL(cudaFree(d_dUinteg));
    CUDA_SAFE_CALL(cudaFree(d_UintegF));
    CUDA_SAFE_CALL(cudaFree(d_s));
    CUDA_SAFE_CALL(cudaFree(d_sJ));
    CUDA_SAFE_CALL(cudaFree(d_S));
    CUDA_SAFE_CALL(cudaFree(d_f));
    CUDA_SAFE_CALL(cudaFree(d_fJ));
    CUDA_SAFE_CALL(cudaFree(d_F));
    CUDA_SAFE_CALL(cudaFree(d_q));
    CUDA_SAFE_CALL(cudaFree(d_Q));
    status = cublasShutdown();
#endif
  };

  // Main solver function
  void dg_solver(scalar* U, scalar* f_rk){

    // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
    if(_multifluid)   Lcpu_mapToFace_multifluid(_M_s, _M_T, _N_F, _N_s, _boundaryMap, U, arch(UF));
    else if(_passive) Lcpu_mapToFace_passive(_M_s, _M_T, _N_F, _N_s, _boundaryMap, U, arch(UF));
    
    // collocationU: requires phi, dphi, Ustar, Uinteg, dUinteg and some sizes
    if (_blas==1) {
      blasGemm('N','N', _N_G   , _N_E*_N_F, _N_s, 1, arch(phi),  _N_G   , U, _N_s, 0.0, arch(Uinteg), _N_G);
      blasGemm('N','N', _N_G*_D, _N_E*_N_F, _N_s, 1, arch(dphi), _N_G*_D, U, _N_s, 0.0, arch(dUinteg), _N_G*_D);}
    else Lcpu_collocationU(_D, _N_G, _N_s, _N_E, _N_F, arch(Uinteg), arch(dUinteg), arch(phi), arch(dphi), U);
    
    // collocationUF: requires psi, UF, UintegF and some sizes
    blasCopy(2*_N_F*_M_T, arch(UF), 1, arch(UintegF), 1);
      
    // evaluate_sf: requires Uinteg, (dUintegR), H0, G0, s,f
    if(_multifluid) Lcpu_evaluate_sf_multifluid(_D, _N_G, _N_E, _N_F, _model, arch(s), arch(f), arch(Uinteg), arch(dUinteg), arch(invJac));
    if(_passive)    Lcpu_evaluate_sf_passive(_D, _N_G, _N_E, _N_F, _gamma0, arch(s), arch(f), arch(Uinteg), arch(dUinteg), arch(invJac));
      
    // evaluate_q: requires UintegF, normals, q, H0, G0
    if(_multifluid) Lcpu_evaluate_q_multifluid(_M_G, _M_T, _N_F, _flux, _model, arch(q), arch(UintegF));
    if(_passive)    Lcpu_evaluate_q_passive(_M_G, _M_T, _N_F, _flux, _gamma0, arch(q), arch(UintegF));
      
    // redistribute_sf: requires J, invJac, s, f, phi_w, dphi_w, sJ, fJ, S, F
    Lcpu_redistribute_sf(_D, _N_G, _N_E, _N_F, arch(sJ), arch(fJ), arch(s), arch(f), arch(J), arch(invJac));
      
    // matrix-matrix multiply for sf
    if (_blas==1)  {
      blasGemm('T','N', _N_s, _N_E*_N_F, _N_G   , 1, arch(phi_w) , _N_G   , arch(sJ), _N_G  , 0.0, arch(S), _N_s);
      blasGemm('T','N', _N_s, _N_E*_N_F, _N_G*_D, 1, arch(dphi_w), _N_G*_D, arch(fJ), _N_G*_D, 0.0, arch(F), _N_s);}
    else Lcpu_gemm_sf(_D, _N_G, _N_s, _N_E, _N_F, arch(S), arch(F), arch(sJ), arch(fJ), arch(phi_w), arch(dphi_w));
      
    // map_q: requires map, Qtcj, Q (might want to do this in the previous step)
    Lcpu_mapToElement(_N_s, _N_E, _N_F, arch(Q), arch(q));

    // Make f_rk = S+F+Q
    Lcpu_addSFQ(_N_s, _N_E, _N_F, f_rk, arch(S), arch(F), arch(Q));
    
  }; // end solver function

};

#endif // DG_SOLVER_H
