//
// DG solver class
//
#ifndef DG_SOLVER_H
#define DG_SOLVER_H

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
  scalar* _phi     ; 
  scalar* _phi_w   ; 
  scalar* _dphi    ; 
  scalar* _dphi_w  ;
  scalar* _J       ;
  scalar* _invJac  ; 
  scalar* _UF      ; 
  scalar* _Uinteg  ; 
  scalar* _dUinteg ; 
  scalar* _UintegF ; 
  scalar* _s       ; 
  scalar* _sJ      ; 
  scalar* _S       ; 
  scalar* _f       ; 
  scalar* _fJ      ; 
  scalar* _F       ; 
  scalar* _q       ; 
  scalar* _Q       ; 

  // To calculate the conservation of certain fields
  std::string consfile;
  FILE *consf;
  scalar* _UgC;
  scalar* _phiC;
  scalar* _JC;
  scalar* _I;
  scalar* _weight; // integration weights
  
 public:
  // constructor
 DG_SOLVER(int D, int N_F, int N_E, int N_s, int N_G, int M_T, int M_s, int M_G,
	   scalar* phi, scalar* dphi, scalar* phi_w, scalar* dphi_w, scalar* J, scalar* invJac, scalar* weight,
	   int boundaryMap, int flux, int model, scalar gamma0, int blas, bool multifluid, bool passive) :
  _D(D), _N_F(N_F), _N_E(N_E), _N_s(N_s), _N_G(N_G), _M_T(M_T), _M_s(M_s), _M_G(M_G),
    _boundaryMap(boundaryMap), _flux(flux), _model(model), _gamma0(gamma0), _blas(blas), _multifluid(multifluid), _passive(passive) {

#ifdef USE_CPU
    
    _phi     = new scalar[N_G*N_s];          makeZero(_phi,N_G*N_s);
    _phi_w   = new scalar[N_G*N_s];          makeZero(_phi_w,N_G*N_s);          
    _dphi    = new scalar[D*N_G*N_s];	     makeZero(_dphi,D*N_G*N_s);	 
    _dphi_w  = new scalar[D*N_G*N_s];	     makeZero(_dphi_w,D*N_G*N_s);
    _J       = new scalar[N_E];              makeZero(_J,N_E);                                 // not same as J!!
    _invJac  = new scalar[N_G*D*N_E*D];      makeZero(_invJac,N_G*D*N_E*D);                    // not same as invJac!!
    _UF      = new scalar[2*N_F*M_s*M_T];    makeZero(_UF,2*N_F*M_s*M_T); 
    _Uinteg  = new scalar[N_F*N_G*N_E];      makeZero(_Uinteg,N_F*N_G*N_E);	 
    _dUinteg = new scalar[D*N_G*N_E*N_F];    makeZero(_dUinteg,D*N_G*N_E*N_F); 
    _UintegF = new scalar[2*N_F*M_G*M_T];    makeZero(_UintegF,2*N_F*M_G*M_T); 
    _s       = new scalar[N_G*N_E*N_F];      makeZero(_s,N_G*N_E*N_F);	 
    _sJ      = new scalar[N_G*N_E*N_F];      makeZero(_sJ,N_G*N_E*N_F);	 
    _S       = new scalar[N_s*N_E*N_F];      makeZero(_S,N_s*N_E*N_F);	 
    _f       = new scalar[D*N_F*N_G*N_E];    makeZero(_f,D*N_F*N_G*N_E); 
    _fJ      = new scalar[D*N_G*N_E*N_F];    makeZero(_fJ,D*N_G*N_E*N_F); 
    _F       = new scalar[N_s*N_E*N_F];      makeZero(_F,N_s*N_E*N_F);	 
    _q       = new scalar[M_G*M_T*N_F*2];    makeZero(_q,M_G*M_T*N_F*2); 
    _Q       = new scalar[N_s*N_E*N_F];      makeZero(_Q,N_s*N_E*N_F);   
    
    memcpy(_phi,    phi,    N_G*N_s*sizeof(scalar));
    memcpy(_phi_w,  phi_w,  N_G*N_s*sizeof(scalar));
    memcpy(_dphi,   dphi,   D*N_G*N_s*sizeof(scalar));
    memcpy(_dphi_w, dphi_w, D*N_G*N_s*sizeof(scalar));
    memcpy(_J, J, N_E*sizeof(scalar));
    memcpy(_invJac, invJac, N_G*D*N_E*D*sizeof(scalar));
	   
#elif USE_GPU
    // Allocate space on the GPU
    CUDA_SAFE_CALL(cudaMalloc((void**) &_phi,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_phi_w,N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_dphi,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_dphi_w,D*N_G*N_s*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_J,N_E*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_invJac,N_G*D*N_E*D*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &_UF,M_s*M_T*N_F*2*sizeof(scalar)));

    CUDA_SAFE_CALL(cudaMalloc((void**) &_Uinteg,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_dUinteg,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_UintegF,M_G*M_T*N_F*2*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &_s,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_sJ,N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_S,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &_f,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_fJ,D*N_G*N_E*N_F*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_F,N_s*N_E*N_F*sizeof(scalar)));
    
    CUDA_SAFE_CALL(cudaMalloc((void**) &_q,M_G*M_T*N_F*2*sizeof(scalar)));
    CUDA_SAFE_CALL(cudaMalloc((void**) &_Q,N_s*N_E*N_F*sizeof(scalar)));


    // Send the stuff to the device
    CUDA_SAFE_CALL(cudaMemcpy(_phi, phi, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_phi_w, phi_w, N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_dphi, dphi, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_dphi_w, dphi_w, D*N_G*N_s*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_J, J, N_E*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_invJac, invJac, N_G*D*N_E*D*sizeof(scalar), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemset(_Q, (scalar)0.0, N_E*N_F*N_s*sizeof(scalar)));

#endif

    // Initialize some stuff for conservation calculations
    consfile = "conservation.dat";
    consf    = fopen(consfile.c_str(),"w");
    _UgC     = new scalar[N_G*N_E*N_F];  makeZero(_UgC,N_G*N_E*N_F);
    _phiC    = new scalar[N_G*N_s];      memcpy(_phiC,phi,N_G*N_s*sizeof(scalar));
    _JC      = new scalar[N_E];          memcpy(_JC,J,N_E*sizeof(scalar));
    _I       = new scalar[N_F];          makeZero(_I,N_F);
    _weight  = new scalar[N_G];          memcpy(_weight,weight,_N_G*sizeof(scalar));
  };

  // destructor
  ~DG_SOLVER(){
#ifdef USE_CPU
    delete[] _phi;
    delete[] _phi_w;
    delete[] _dphi;
    delete[] _dphi_w;
    delete[] _J;
    delete[] _invJac;
    delete[] _UF;
    delete[] _Uinteg;
    delete[] _dUinteg;
    delete[] _UintegF;
    delete[] _s;
    delete[] _sJ;
    delete[] _S;
    delete[] _f;
    delete[] _fJ;
    delete[] _F;
    delete[] _q;
    delete[] _Q;
#elif USE_GPU
    CUDA_SAFE_CALL(cudaFree(_phi));
    CUDA_SAFE_CALL(cudaFree(_phi_w));
    CUDA_SAFE_CALL(cudaFree(_dphi));
    CUDA_SAFE_CALL(cudaFree(_dphi_w));
    CUDA_SAFE_CALL(cudaFree(_J));
    CUDA_SAFE_CALL(cudaFree(_invJac));
    CUDA_SAFE_CALL(cudaFree(_UF));
    CUDA_SAFE_CALL(cudaFree(_Uinteg));
    CUDA_SAFE_CALL(cudaFree(_dUinteg));
    CUDA_SAFE_CALL(cudaFree(_UintegF));
    CUDA_SAFE_CALL(cudaFree(_s));
    CUDA_SAFE_CALL(cudaFree(_sJ));
    CUDA_SAFE_CALL(cudaFree(_S));
    CUDA_SAFE_CALL(cudaFree(_f));
    CUDA_SAFE_CALL(cudaFree(_fJ));
    CUDA_SAFE_CALL(cudaFree(_F));
    CUDA_SAFE_CALL(cudaFree(_q));
    CUDA_SAFE_CALL(cudaFree(_Q));
#endif
    delete[] _UgC;
    delete[] _I;
    delete[] _weight;
    fclose(consf);
  };

  // Main solver function
  void dg_solver(scalar* U, scalar* f_rk){

    // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
    if(_multifluid)   Lcpu_mapToFace_multifluid(_M_s, _M_T, _N_F, _N_s, _boundaryMap, U, _UF);
    else if(_passive) Lcpu_mapToFace_passive(_M_s, _M_T, _N_F, _N_s, _boundaryMap, U, _UF);
    
    // collocationU: requires phi, dphi, Ustar, Uinteg, dUinteg and some sizes
    if (_blas==1) {
      blasGemm('N','N', _N_G   , _N_E*_N_F, _N_s, 1, _phi,  _N_G   , U, _N_s, 0.0, _Uinteg, _N_G);
      blasGemm('N','N', _N_G*_D, _N_E*_N_F, _N_s, 1, _dphi, _N_G*_D, U, _N_s, 0.0, _dUinteg, _N_G*_D);}
    else Lcpu_collocationU(_D, _N_G, _N_s, _N_E, _N_F, _Uinteg, _dUinteg, _phi, _dphi, U);
    
    // collocationUF: requires psi, UF, UintegF and some sizes
    blasCopy(2*_N_F*_M_T, _UF, 1, _UintegF, 1);
      
    // evaluate_sf: requires Uinteg, (dUintegR), H0, G0, s,f
    if(_multifluid) Lcpu_evaluate_sf_multifluid(_D, _N_G, _N_E, _N_F, _model, _s, _f, _Uinteg, _dUinteg, _invJac);
    if(_passive)    Lcpu_evaluate_sf_passive(_D, _N_G, _N_E, _N_F, _gamma0, _s, _f, _Uinteg, _dUinteg, _invJac);
      
    // evaluate_q: requires UintegF, normals, q, H0, G0
    if(_multifluid) Lcpu_evaluate_q_multifluid(_M_G, _M_T, _N_F, _flux, _model, _q, _UintegF);
    if(_passive)    Lcpu_evaluate_q_passive(_M_G, _M_T, _N_F, _flux, _gamma0, _q, _UintegF);
      
    // redistribute_sf: requires J, invJac, s, f, phi_w, dphi_w, sJ, fJ, S, F
    Lcpu_redistribute_sf(_D, _N_G, _N_E, _N_F, _sJ, _fJ, _s, _f, _J, _invJac);
      
    // matrix-matrix multiply for sf
    if (_blas==1)  {
      blasGemm('T','N', _N_s, _N_E*_N_F, _N_G   , 1, _phi_w , _N_G   , _sJ, _N_G  , 0.0, _S, _N_s);
      blasGemm('T','N', _N_s, _N_E*_N_F, _N_G*_D, 1, _dphi_w, _N_G*_D, _fJ, _N_G*_D, 0.0, _F, _N_s);}
    else Lcpu_gemm_sf(_D, _N_G, _N_s, _N_E, _N_F, _S, _F, _sJ, _fJ, _phi_w, _dphi_w);
      
    // map_q: requires map, Qtcj, Q (might want to do this in the previous step)
    Lcpu_mapToElement(_N_s, _N_E, _N_F, _Q, _q);

    // Make f_rk = S+F+Q
    Lcpu_addSFQ(_N_s, _N_E, _N_F, f_rk, _S, _F, _Q);
    
  }; // end solver function


  // Function to calculate conservation of certain quantities
  void conservation(scalar* U, double time){

    // Collocate the solution to the integration points
    hostblasGemm('N','N', _N_G, _N_E*_N_F, _N_s, 1, _phiC, _N_G, U, _N_s, 0.0, _UgC, _N_G);

    // Take the cell average of the solution
    makeZero(_I, _N_F);
    for(int fc = 0; fc < _N_F; fc++){
      for(int e = 0; e < _N_E; e++){
    	for(int g = 0; g < _N_G; g++){
    	  _I[fc] += _UgC[(e*_N_F+fc)*_N_G+g]*_JC[e]*_weight[g];
    	}
      }
    }
    // write to file
    fprintf(consf,"%20.16E\t", time); for(int fc = 0; fc < _N_F; fc++) fprintf(consf,"%20.16E\t", _I[fc]); fprintf(consf,"\n");

  };//end conservation function

};

#endif // DG_SOLVER_H
