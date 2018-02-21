/*!
  \file dg_solver.h  
  \brief DG solver class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef DG_SOLVER_H
#define DG_SOLVER_H

#include "physics.h" 
#include "vis_physics.h" //viscous/diffusive flux functions
#include "boundaries.h"
#include <stdio.h>
#include "simpleMesh.h"
#include "timers.h"
#include "mem_counter.h"
#include "kernels.h"
#include "kernels_phil.h" //additional kernel files for mixed-form approach
#include "BR2_tools.h" //This is the header file for BR2 procedures
#include "communicator.h" //12/2/2017: Must have this, unfortunately, for AD approach
#ifdef USE_GPU
#include <cublas.h>
#endif
#ifdef USE_MPI
#include "mpi.h"
#endif


class DG_SOLVER
{
 private:
  
  int _N_E;
  int _N_s;
  int _N_G;
  int _N_N;
  int _M_T;
  int _M_s;
  int _M_G;
  int _M_B;
  int _Ne_AUG; //N_E plus number of ghost elements for the processor
  // COMMUNICATOR _communicator;

 
  int* _map;
  int* _invmap;
  int* _boundaryMap;
  int _rflctiveIdx;
  int _otheroneIdx;
  int _otheronestart;
  //BR2 Edit: No slip boundary considerations:
  int _noslpIdx;
  int _noslpstrt;
  int _noslpcount;
  //zero gradient boundary considerations: (Actually, just farfield; zero gradient may have been a mistake)
  int _nogradIdx;
  int _nogradstrt;
  int _nogradcount;
  //Simple inflow boundary (diesignated as inflow A)
  int _AnflwIdx;
  int _Anflwstrt;
  int _Anflwcount;
  //Kushner Jet inflow boundary
  int _KJetIdx;
  int _KJetstrt;
  int _KJetcount;
  //Homogeneous boundary (should only be used for scalarAD case)
  int _HomoIdx;
  int _Homostrt;
  int _Homocount;
  //Subsonic Ouflow boundary 
  int _SubOutIdx;
  int _SubOutstrt;
  int _SubOutcount;
  //END BR2 Edit
  scalar* _phi     ; 
  scalar* _phi_w   ; 
  scalar* _dphi    ; 
  scalar* _dphi_w  ;
  scalar* _psi     ;
  scalar* _psi_w   ;
  //scalar* _xyz     ;
  //scalar* _xyzf    ;
  scalar* _J       ;
  scalar* _invJac  ;
  scalar* _JF      ;
  scalar* _normals ;
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
  scalar* _qJ      ;
  scalar* _Qtcj    ;
  scalar* _Q       ; 
  //int* _sensor     ; //PEJ 10/23/2017 for smart AD approach
  
  //BR2 Edit: 01/05/2016: Declarations for BR2 considerations:
  //For clarification, talk to Phil
  //scalar* _phigF        ; //trace of element shape functions at each local interface. IMPORTED  
  //scalar* _dphigF       ; //gradient of those traces IMPORTED
  int*    _BR2_Map      ; //a small structure linking faces of reference element to mesh interfaces. IMPORTED
  //scalar* _Jac_trace    ; //the Jacobian matrix of each element at its traces. IMPORTED
  //scalar* _invJac_trace ; //invese of the Jacobian matrix of each element at its traces. IMPORTED
  //scalar* _BR2_resi     ; //the residual for the BR2 system, based on the DG solution's jumps across each interface
  //scalar* _BR2_poly     ; //the DOF for the polynomial expansion of the gradient correction for each element
  //scalar* _dUinteg_phys ; //Corrected gradient of the solution U, at interior quadrature points, wrt physical coordinates
  //scalar* _dUintegF_phys;  //Corrected gradient of the solution U, at interface quadrature points, from each side of the interface, wrt physical coordinates
  //scalar* _M_inv        ; //inverse mass matrix, necessary to calculate BR2_resi. IMPORTED
  //scalar* _weight_trace ; //stores the gaussian quadrature weights along the element intersections, in specific direction based on edge and cell
  //END BR2 Edit

  //PEJ 05/23/2017 - ? : Declarations for general mixed-form approach
  scalar* _PSIxR_Global; //the recovery operator, in serial storage, interface-oriented
  scalar* _PSIxR_elemwise; //the recovery operator again in element-wise storage pattern
  scalar* _SigFace_from_Uhat_elemwise; //the interface gradient calculation in element-wise storage pattern
  scalar* _SigFace_from_Uhat; //the interface gradient calculation in interface-wise storage pattern
  int* _RecoPair; //another element-interface pairing structure
  int* _BinarySideAddress; //either 1 or zero, telling an element if it is elA or elB on a given local face
  scalar* _UhCommon; //the common solution values on each interface
  scalar* _UhCommonHalf; //the common solution values on each interface
  scalar * _sigma; //the gradient approximation over element interiors
  scalar* _serial_MSigVol; //gets element contribution to sigma
  scalar* _serial_MSigSurf; //gets interface contribution to sigma
  int* _Alt_FaceFromElem; //tells each element its interface addresses
  //scalar* _serial_SigFace_from_DOF; //gets interface gradient from DOF of 2 nieghoring elements
  scalar* _gradCommon; //common gradient on interfaces
  scalar* _gradCommonHalf; //partial common gradient on interfaces
  scalar* _PSIxR_biased_Global; //biased recovery operator
  scalar* _PSIxR_biased_elemwise; //biased recovery operator
  //scalar* _Uicb; //the competing icb solution values on each interface
  scalar* _UicbHalf; //the neighboring element contributions to icb solution
  //PEJ 10/01/2017: Bring in detJ by quadrature weights
  scalar* _detJxW_full;
  scalar* _detJ_full;
  //PEJ 10/11/2017; Some boundary gradient stuff:
  scalar* _serial_Uhat2GradBC;
  scalar* _serial_UhCommon2GradBC;
  //PEJ 10/15/2017: Artificial dissipation equipment 
 
  scalar* _ADeps; //copied in from main.cc, D entries per element
  scalar* _ADepsF; //copied in from main.cc, D entries per interface
  scalar* _NaiveGradCommonHalf; //broken gradient of DG solution at interfaces
  scalar* _NaiveGradCommon; //summed at interface
  scalar* _DivMax; //single value, max divergence in domain
  scalar _Beta_S; //copied in from main.cc, it's the single-value AD strenght (set to 0 for no AD)
  scalar _Mew_S; //copied in from main.cc, it's the single-value AD spread factor
  scalar _Cthresh; //copied in from main.cc, it's the cutoff C magnitude for AD treatment
  scalar* _CsMax; //maximum C variable in the domain, D entries
  scalar* _LamMax; //max directional wavespeed in the domain, D entries
  scalar _epsGen; //small epsilon to avoid division by zero, copted in from main.cc
  scalar* _ADepsMax; //D entries, max element size in each direction
  scalar* _elemCmax; //maximum C value (direction-independent) in each element
  
  //All-partition AD parameter maxima:
  scalar* GLO_DivMax;
  scalar* GLO_LamMax;
  scalar* GLO_CsMax;

  TIMERS &_timers;
  COMMUNICATOR &_communicator; //PEJ 12/11/2017

  // To calculate the conservation of certain fields
  std::string consfile;
  FILE *consf;
  scalar* _UgC;
  scalar* _phiC;
  scalar* _JC;
  scalar* _I;
  scalar* _weight; // integration weights

  //PEJ 10/01/2017: To calculate the taylor-green statistics
  scalar* _EKlocal;
#ifdef TGVSTATS
  FILE*TGVstatsf;
#endif
  
 public:
  /*!
    \brief Constructor
    \param[in] N_E number of elements
    \param[in] N_s number of nodes per element
    \param[in] N_G number of gaussian nodes per element
    \param[in] N_N number of neighbors per element
    \param[in] M_T number of interfaces
    \param[in] M_s number of nodes per interface
    \param[in] M_G number of gaussian nodes per interface
    \param[in] M_B number of boundaries
    \param[in] map map from elements to interfaces
    \param[in] invmap map from interfaces to elements
    \param[in] phi element polynomial basis
    \param[in] phi_w element polynomial basis (premultiplied by weights)
    \param[in] dphi element polynomial basis derivative
    \param[in] dphi_w element polynomial basis derivative (premultiplied by weights)
    \param[in] psi interface polynomial basis
    \param[in] psi_w interface polynomial basis (premultiplied by weights)
    \param[in] J jacobian of each element
    \param[in] invJac inverse jacobian of each element
    \param[in] JF jacobian of each interface
    \param[in] weight integration weights
    \param[in] normals normals to all the interfaces
    \param[in] m mesh we are operating on
    \param[in] timers timers to use
    \param[in] mem_counter memory counter to use
  */
  //BR2 Edit: DG_solver with room for extra arguments
 DG_SOLVER(int N_E, int N_s, int N_G,  int N_N, int M_T, int M_s, int M_G, int M_B, int Ne_AUG, COMMUNICATOR &communicator, int myid,
	   int* map, int* invmap, scalar* phi, scalar* dphi, scalar* phi_w, scalar* dphi_w, scalar* psi, scalar* psi_w, //scalar* xyz, scalar* xyzf,
	    scalar* J, scalar* invJac, scalar* JF, scalar* weight, scalar* normals, simpleMesh &m, TIMERS &timers, MEM_COUNTER &mem_counter,
	    /*scalar* Jac_trace,*/ /*scalar* invJac_trace,*/ /*scalar* phigF,*/ /*scalar* dphigF,*/ int* BR2_Map, /*scalar* M_inv,*/ /*scalar* weight_trace,*/  
	   scalar* PSIxR_Global, scalar* PSIxR_elemwise, scalar* SigFace_from_Uhat_elemwise, int* RecoPair, int* BinarySideAddress, scalar* serial_MSigVol, scalar* serial_MSigSurf, int* Alt_FaceFromElem, scalar* SigFace_from_Uhat, scalar* PSIxR_biased_Global, scalar* PSIxR_biased_elemwise, scalar* detJxW_full, scalar* detJ_full, scalar* serial_Uhat2GradBC, scalar* serial_UhCommon2GradBC, scalar* ADeps, scalar* ADepsF, scalar Beta_S, scalar epsGen, scalar Mew_S, scalar Cthresh) :
  
  _N_E(N_E), _N_s(N_s), _N_G(N_G), _N_N(N_N), _M_T(M_T), _M_s(M_s), _M_G(M_G), _M_B(M_B), _Ne_AUG(Ne_AUG), _communicator(communicator), _timers(timers), _Beta_S(Beta_S), _epsGen(epsGen), _Mew_S(Mew_S), _Cthresh(Cthresh){
    //END BR2 Edit
    //The old DG_solver 
    /*
      DG_SOLVER(int N_E, int N_s, int N_G,  int N_N, int M_T, int M_s, int M_G, int M_B, 
      int* map, int* invmap, scalar* phi, scalar* dphi, scalar* phi_w, scalar* dphi_w, scalar* psi, scalar* psi_w, //scalar* xyz, scalar* xyzf,
      scalar* J, scalar* invJac, scalar* JF, scalar* weight, scalar* normals, simpleMesh &m, TIMERS &timers, MEM_COUNTER &mem_counter) :
      _N_E(N_E), _N_s(N_s), _N_G(N_G), _N_N(N_N), _M_T(M_T), _M_s(M_s), _M_G(M_G), _M_B(M_B), _timers(timers) {
    */
    
    // Indexes for boundary conditions
    //printf("Inside DG solver, about to do the boundary thing\n");
    //printf("Supposed number of boundary interfaces = _M_B = %d\n", _M_B);
    int* boundaryIdx = m.getBoundaryIdx();
    //printf("Succesfully retrived the boundaryIdx pointer array\n");
    _rflctiveIdx = boundaryIdx[0];       // number of reflective interfaces
    //printf("Number of reflective interfaces = rflctiveIdx = %d\n", _rflctiveIdx);
    _otheroneIdx = _rflctiveIdx;
    _otheronestart = _rflctiveIdx;
    //_farfieldIdx = boundaryIdx[1]-boundaryIdx[0]; // number of farfield interfaces
    //_farfieldstart = boundaryIdx[1];
    
    
    //PEJ Edit 01/20/2016:  
    //printf("About to populate _noslpstrt, _noslpIdx\n");
    _noslpstrt = _rflctiveIdx; //end of reflective section, start of no slip section
    _noslpIdx = boundaryIdx[1]; //end of no slip section
    _noslpcount = _noslpIdx - _noslpstrt;

    _nogradstrt = _noslpIdx; //end of noslip section, start of zero normal gradient section
    _nogradIdx = boundaryIdx[2]; //end of zero normal gradient section
    _nogradcount = _nogradIdx - _nogradstrt;

    _Anflwstrt = _nogradIdx; //end of nograd section, start of A inflow section
    _AnflwIdx = boundaryIdx[3]; //end of A inflow section
    _Anflwcount = _AnflwIdx - _Anflwstrt; //number of A inflow interfaces

    _KJetstrt = _AnflwIdx; //end of Anflw section, start of KJet inflow section
    _KJetIdx = boundaryIdx[4]; //end of KJet section
    _KJetcount = _KJetIdx - _KJetstrt; //number of KJet inflow interfaces
    
    _Homostrt = _KJetIdx; //end of KJet segment, start of homogeneous segment
    _HomoIdx = boundaryIdx[5]; //end of Homogensous BC section
    _Homocount = _HomoIdx - _Homostrt; //number of homogeneous BC interfaces

    _SubOutstrt = _HomoIdx; //end of Homo segment, start of Subsonic Outflow segment
    _SubOutIdx = boundaryIdx[6]; //end of Subsonic outflow BC section
    _SubOutcount = _SubOutIdx - _SubOutstrt; //number of Subsonic Outflow BC interfaces
    
    boundaryIdx = NULL;
    //Print some boundary information
    printf("\npr=%d: Boundary information, inside DG_Solver\n", myid);
    printf("Supposed number of boundary interfaces = _M_B = %d\n", _M_B);
    printf("Number of reflective interfaces = rflctiveIdx = %d\n", _rflctiveIdx - 0);
    printf("Number of placeholder interfaces = otheroneIdx = %d\n", _otheroneIdx-_otheronestart);
    printf("Number of no-slip interfaces = noslpIdx = %d\n", _noslpcount);
    printf("Number of no-grad interfaces = nogradIdx = %d\n", _nogradcount);
    printf("Number of A inflow interfaces = Anflwcount = %d\n", _Anflwcount);
    printf("Number of KJet inflow interfaces = KJetcount = %d\n", _KJetcount);
    printf("Number of Homogeneous BC interfaces = Homocount = %d\n", _Homocount);
    printf("Number of Subsonic Outflow BC interfaces = SubOutcount = %d\n", _SubOutcount);
    //END PEJ EDIT
    int Autocheck = 0;
    printf("pr=%d: DGsolver check 1\n", myid);
#ifdef USE_CPU
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _map     = new int[M_s*M_T*N_F*2];                                             mem_counter.addToCPUCounter(M_s*M_T*N_F*2*sizeof(int));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _invmap  = new int[M_s*N_N*N_E*N_F*2];                                         mem_counter.addToCPUCounter(M_s*N_N*N_E*N_F*2*sizeof(int)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _boundaryMap  = new int[M_B];                                                  mem_counter.addToCPUCounter(M_B*sizeof(int)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _phi     = new scalar[N_G*N_s];          makeZero(_phi,N_G*N_s);               mem_counter.addToCPUCounter(N_G*N_s*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _phi_w   = new scalar[N_G*N_s];          makeZero(_phi_w,N_G*N_s);             mem_counter.addToCPUCounter(N_G*N_s*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _dphi    = new scalar[D*N_G*N_s];	     makeZero(_dphi,D*N_G*N_s);	           mem_counter.addToCPUCounter(D*N_G*N_s*sizeof(scalar));  
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _dphi_w  = new scalar[D*N_G*N_s];	     makeZero(_dphi_w,D*N_G*N_s);          mem_counter.addToCPUCounter(D*N_G*N_s*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _psi     = new scalar[M_G*M_s];	     makeZero(_psi,M_G*M_s);               mem_counter.addToCPUCounter(M_G*M_s*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _psi_w   = new scalar[M_G*M_s];	     makeZero(_psi_w,M_G*M_s);             mem_counter.addToCPUCounter(M_G*M_s*sizeof(scalar)); 
    //_xyz     = new scalar[D*N_E*N_G];	     makeZero(_xyz,D*N_E*N_G);
    //_xyzf    = new scalar[D*M_T*M_G];	     makeZero(_xyzf,D*M_T*M_G);
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _J       = new scalar[N_E];              makeZero(_J,N_E);                     mem_counter.addToCPUCounter(N_E*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _invJac  = new scalar[N_G*D*N_E*D];      makeZero(_invJac,N_G*D*N_E*D);        mem_counter.addToCPUCounter(N_G*D*N_E*D*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _JF = new scalar[2*M_T];     	     makeZero(_JF,2*M_T);	           mem_counter.addToCPUCounter(2*M_T*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _normals = new scalar[D*M_T];	     makeZero(_normals,D*M_T);	           mem_counter.addToCPUCounter(D*M_T*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _UF      = new scalar[2*N_F*M_s*M_T];    makeZero(_UF,2*N_F*M_s*M_T);          mem_counter.addToCPUCounter(2*N_F*M_s*M_T*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _Uinteg  = new scalar[N_F*N_G*Ne_AUG];      makeZero(_Uinteg,N_F*N_G*Ne_AUG);	   mem_counter.addToCPUCounter(N_F*N_G*Ne_AUG*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _dUinteg = new scalar[D*N_G*N_E*N_F];    makeZero(_dUinteg,D*N_G*N_E*N_F);     mem_counter.addToCPUCounter(D*N_G*N_E*N_F*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _UintegF = new scalar[2*N_F*M_G*M_T];    makeZero(_UintegF,2*N_F*M_G*M_T);     mem_counter.addToCPUCounter(2*N_F*M_G*M_T*sizeof(scalar)); 
    //printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _s       = new scalar[N_G*N_E*N_F];      makeZero(_s,N_G*N_E*N_F);	           mem_counter.addToCPUCounter(N_G*N_E*N_F*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _sJ      = new scalar[N_G*N_E*N_F];      makeZero(_sJ,N_G*N_E*N_F);	           mem_counter.addToCPUCounter(N_G*N_E*N_F*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _S       = new scalar[N_s*N_E*N_F];      makeZero(_S,N_s*N_E*N_F);	           mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar)); 
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _f       = new scalar[D*N_F*N_G*N_E];    makeZero(_f,D*N_F*N_G*N_E);           mem_counter.addToCPUCounter(D*N_F*N_G*N_E*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _fJ      = new scalar[D*N_G*N_E*N_F];    makeZero(_fJ,D*N_G*N_E*N_F);          mem_counter.addToCPUCounter(D*N_G*N_E*N_F*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _F       = new scalar[N_s*N_E*N_F];      makeZero(_F,N_s*N_E*N_F);	           mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _q       = new scalar[M_G*M_T*N_F*2];    makeZero(_q,M_G*M_T*N_F*2);           mem_counter.addToCPUCounter(M_G*M_T*N_F*2*sizeof(scalar)); 
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _qJ      = new scalar[M_G*M_T*N_F*2];    makeZero(_qJ,M_G*M_T*N_F*2);          mem_counter.addToCPUCounter(M_G*M_T*N_F*2*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _Qtcj    = new scalar[M_s*M_T*N_F*2];    makeZero(_Qtcj,M_s*M_T*N_F*2);        mem_counter.addToCPUCounter(M_s*M_T*N_F*2*sizeof(scalar)); 
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _Q       = new scalar[N_s*N_E*N_F];      makeZero(_Q,N_s*N_E*N_F);             mem_counter.addToCPUCounter(N_s*N_E*N_F*sizeof(scalar));       
    //BR2 Edit 01/05/2016: Extra allocation for BR2 applications
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_dUinteg_phys  = new scalar[N_E*N_F*N_G*D];     makeZero(_dUinteg_phys, N_E*N_F*N_G*D);    mem_counter.addToCPUCounter(N_E*N_F*N_G*D*sizeof(scalar));
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_dUintegF_phys = new scalar[M_T*2*N_F*M_G*D];   makeZero(_dUintegF_phys, M_T*2*N_F*M_G*D); mem_counter.addToCPUCounter(M_T*2*N_F*M_G*D*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_BR2_poly      = new scalar[N_E*N_N*N_F*D*N_s]; makeZero(_BR2_poly, N_E*N_N*N_F*D*N_s);    mem_counter.addToCPUCounter(N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_BR2_resi      = new scalar[N_E*N_N*N_F*D*N_s]; makeZero(_BR2_resi, N_E*N_N*N_F*D*N_s);    mem_counter.addToCPUCounter(N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_invJac_trace  = new scalar[M_G*D*D*N_N*N_E*D];   makeZero(_invJac_trace, M_G*D*D*N_N*N_E*D);  mem_counter.addToCPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_Jac_trace     = new scalar[M_G*D*D*N_N*N_E*D];   makeZero(_Jac_trace, M_G*D*D*N_N*N_E*D);  mem_counter.addToCPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_phigF         = new scalar[M_G*D*N_N*N_s];       makeZero(_phigF, M_G*D*N_N*N_s);             mem_counter.addToCPUCounter(M_G*D*N_N*N_s*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_dphigF        = new scalar[D*M_G*D*N_N*N_s];     makeZero(_dphigF, D*M_G*D*N_N*N_s);          mem_counter.addToCPUCounter(D*M_G*D*N_N*N_s*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _BR2_Map       = new int[M_T*4];                                                           mem_counter.addToCPUCounter(M_T*4*sizeof(int));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_M_inv         = new scalar[N_E*N_s*N_s];       makeZero(_M_inv, N_E*N_s*N_s);             mem_counter.addToCPUCounter(N_E*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_weight_trace  = new scalar[D*N_N*M_G];         makeZero(_weight_trace, D*N_N*M_G);        mem_counter.addToCPUCounter(D*N_N*M_G*sizeof(scalar));
    //END BR2 edit
    
    //Create space for new scalars
    //PEJ 05/23/2017-?
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _PSIxR_Global =               new scalar[M_T*M_G*2*N_s];         makeZero(_PSIxR_Global, M_T*M_G*2*N_s);                       mem_counter.addToCPUCounter(M_T*M_G*2*N_s*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _PSIxR_elemwise =             new scalar[Ne_AUG*(D+N_N)*M_G*N_s];   makeZero(_PSIxR_elemwise, Ne_AUG*(D+N_N)*M_G*N_s);               mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _SigFace_from_Uhat_elemwise = new scalar[Ne_AUG*(D+N_N)*D*M_G*N_s]; makeZero(_SigFace_from_Uhat_elemwise, Ne_AUG*(D+N_N)*D*M_G*N_s); mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar));
    _SigFace_from_Uhat =    new scalar[M_T*D*M_G*2*N_s];       makeZero(_SigFace_from_Uhat,M_T*D*M_G*2*N_s);           mem_counter.addToCPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _RecoPair =                   new int[Ne_AUG*(D+N_N)];                                                                            mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _BinarySideAddress =          new int[Ne_AUG*(D+N_N)];                                                                            mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
    //   printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _UhCommon     =               new scalar[M_T*M_G*N_F];           makeZero(_UhCommon    , M_T*M_G*N_F);                         mem_counter.addToCPUCounter(M_T*M_G*N_F*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _UhCommonHalf     =           new scalar[M_T*M_G*N_F*2];         makeZero(_UhCommonHalf    , M_T*M_G*N_F*2);                   mem_counter.addToCPUCounter(M_T*M_G*N_F*2*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _sigma =                      new scalar[N_E*D*N_G*N_F];         makeZero(_sigma, N_E*D*N_G*N_F);                              mem_counter.addToCPUCounter(N_E*D*N_G*N_F*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _serial_MSigVol =             new scalar[N_E*D*N_G*N_s];         makeZero(_serial_MSigVol, N_E*D*N_G*N_s);                     mem_counter.addToCPUCounter(N_E*D*N_G*N_s*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _serial_MSigSurf =            new scalar[N_E*D*N_G*N_N*M_G];     makeZero(_serial_MSigSurf, N_E*D*N_G*N_N*M_G);                mem_counter.addToCPUCounter(N_E*D*N_G*N_N*M_G*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _Alt_FaceFromElem =           new int[N_E*N_N];                                                                                mem_counter.addToCPUCounter(N_E*N_N*sizeof(int));
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_serial_SigFace_from_DOF =    new scalar[M_T*D*M_G*2*N_s];       makeZero(_serial_SigFace_from_DOF,M_T*D*M_G*2*N_s);           mem_counter.addToCPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _gradCommon =                 new scalar[M_T*M_G*N_F*D];         makeZero(_gradCommon,M_T*M_G*N_F*D);                          mem_counter.addToCPUCounter(M_T*M_G*N_F*D*sizeof(scalar));
    // printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _gradCommonHalf =             new scalar[M_T*M_G*N_F*D*2];       makeZero(_gradCommonHalf,M_T*M_G*N_F*D*2);                    mem_counter.addToCPUCounter(M_T*M_G*N_F*D*2*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    //_Uicb =                       new scalar[M_T*M_G*2*N_F];         makeZero(_Uicb, M_T*M_G*2*N_F);                               mem_counter.addToCPUCounter(M_T*M_G*2*N_F*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _UicbHalf =                   new scalar[M_T*M_G*2*N_F*2];       makeZero(_UicbHalf, M_T*M_G*2*N_F*2);                         mem_counter.addToCPUCounter(M_T*M_G*2*N_F*2*sizeof(scalar));
    //  printf("AutoCheck %d\n", Autocheck); Autocheck++;
    _PSIxR_biased_Global   = new scalar[M_T*2*M_G*2*N_s]; makeZero(_PSIxR_biased_Global, M_T*2*M_G*2*N_s); mem_counter.addToCPUCounter(M_T*2*M_G*2*N_s*sizeof(scalar));
    _PSIxR_biased_elemwise =      new scalar[Ne_AUG*(D+N_N)*2*M_G*N_s]; makeZero(_PSIxR_biased_elemwise, Ne_AUG*(D+N_N)*2*M_G*N_s);      mem_counter.addToCPUCounter(Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar));
    _detJxW_full =               new scalar[N_E*N_G];               makeZero(_detJxW_full, N_E*N_G);                             mem_counter.addToCPUCounter(N_E*N_G*sizeof(scalar));
    _detJ_full =               new scalar[N_E*N_G];               makeZero(_detJ_full, N_E*N_G);                             mem_counter.addToCPUCounter(N_E*N_G*sizeof(scalar));
    _serial_Uhat2GradBC =     new scalar[M_B*D*M_G*N_s];            makeZero(_serial_Uhat2GradBC, M_B*D*M_G*N_s); mem_counter.addToCPUCounter(M_B*D*M_G*N_s*sizeof(scalar));
    _serial_UhCommon2GradBC =   new scalar[M_B*D*M_G*(N_N*M_G)];    makeZero(_serial_UhCommon2GradBC, M_B*D*M_G*(N_N*M_G)); mem_counter.addToCPUCounter(M_B*D*M_G*(N_N*M_G)*sizeof(scalar));
   
    _ADeps = new scalar[Ne_AUG*D]; makeZero(_ADeps, Ne_AUG*D); mem_counter.addToCPUCounter(Ne_AUG*D*sizeof(scalar));
    _ADepsF = new scalar[M_T*D]; makeZero(_ADepsF, M_T*D); mem_counter.addToCPUCounter(M_T*D*sizeof(scalar));
    _NaiveGradCommonHalf = new scalar[M_T*M_G*N_F*D*2];       makeZero(_NaiveGradCommonHalf,M_T*M_G*N_F*D*2);       mem_counter.addToCPUCounter(M_T*M_G*N_F*D*2*sizeof(scalar));
    _NaiveGradCommon =            new scalar[M_T*M_G*N_F*D];         makeZero(_NaiveGradCommon,M_T*M_G*N_F*D);     mem_counter.addToCPUCounter(M_T*M_G*N_F*D*sizeof(scalar));
    _CsMax  = new scalar[D]; makeZero(_CsMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    _DivMax = new scalar[D]; makeZero(_DivMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    _LamMax = new scalar[D]; makeZero(_LamMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    _ADepsMax = new scalar[D]; makeZero(_ADepsMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    _elemCmax = new scalar[Ne_AUG]; makeZero(_elemCmax, Ne_AUG); mem_counter.addToCPUCounter(Ne_AUG*sizeof(scalar));
    GLO_CsMax  = new scalar[D]; makeZero(GLO_CsMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    GLO_DivMax = new scalar[D]; makeZero(GLO_DivMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    GLO_LamMax = new scalar[D]; makeZero(GLO_LamMax, D); mem_counter.addToCPUCounter(D*sizeof(scalar));
    //END PEJ EDIT

  
   
    //This section is for copyimg data from the outside
    memcpy(_map        , map        , M_s*M_T*N_F*2*sizeof(int));
    memcpy(_invmap     , invmap     , M_s*N_N*N_E*N_F*2*sizeof(int));
    memcpy(_boundaryMap, m.getBoundaryMap(), M_B*sizeof(int));
    memcpy(_phi        , phi        , N_G*N_s*sizeof(scalar));
    memcpy(_phi_w      , phi_w      , N_G*N_s*sizeof(scalar));
    memcpy(_dphi       , dphi       , D*N_G*N_s*sizeof(scalar));
    memcpy(_dphi_w     , dphi_w     , D*N_G*N_s*sizeof(scalar));
    memcpy(_psi        , psi        , M_G*M_s*sizeof(scalar));
    memcpy(_psi_w      , psi_w      , M_G*M_s*sizeof(scalar));
    //memcpy(_xyz        , xyz        , D*N_E*N_G*sizeof(scalar));
    //memcpy(_xyzf       , xyzf       , D*M_T*M_G*sizeof(scalar));
    memcpy(_J          , J          , N_E*sizeof(scalar));
    memcpy(_invJac     , invJac     , N_G*D*N_E*D*sizeof(scalar));
    memcpy(_JF         , JF         , 2*M_T*sizeof(scalar));
    memcpy(_normals    , normals    , D*M_T*sizeof(scalar));
    //BR2 Edit: Copy some incoming information
    //memcpy(_invJac_trace , invJac_trace , M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //memcpy(_Jac_trace    , Jac_trace    , M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //memcpy(_phigF        , phigF        , M_G*D*N_N*N_s*sizeof(scalar));
    //memcpy(_dphigF       , dphigF       , D*M_G*D*N_N*N_s*sizeof(scalar));
    memcpy(_BR2_Map      , BR2_Map      , M_T*4*sizeof(int));
    //memcpy(_M_inv        , M_inv        , N_E*N_s*N_s*sizeof(scalar));
    //memcpy(_weight_trace , weight_trace , D*N_N*M_G*sizeof(scalar));
    //END BR2 Edit
    //PEJ Edit: PSIxR_Global comes from the outside
    memcpy(_PSIxR_Global              , PSIxR_Global              , M_T*M_G*2*N_s*sizeof(scalar));
    memcpy(_PSIxR_elemwise            , PSIxR_elemwise            , Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar));
 
    memcpy(_SigFace_from_Uhat_elemwise, SigFace_from_Uhat_elemwise, Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar));

    memcpy(_SigFace_from_Uhat   , SigFace_from_Uhat   , M_T*D*M_G*2*N_s*sizeof(scalar));
    memcpy(_RecoPair                  , RecoPair                  , Ne_AUG*(D+N_N)*sizeof(int));
 
    
    memcpy(_BinarySideAddress         , BinarySideAddress         , Ne_AUG*(D+N_N)*sizeof(int));
    memcpy(_serial_MSigSurf           , serial_MSigSurf           , N_E*D*N_G*N_N*M_G*sizeof(scalar));
    memcpy(_serial_MSigVol            , serial_MSigVol            , N_E*D*N_G*N_s*sizeof(scalar));
    memcpy(_Alt_FaceFromElem          , Alt_FaceFromElem          , N_E*N_N*sizeof(int));
 
    //memcpy(_serial_SigFace_from_DOF   , serial_SigFace_from_DOF   , M_T*D*M_G*2*N_s*sizeof(scalar));
    memcpy(_PSIxR_biased_Global       , PSIxR_biased_Global       , M_T*2*M_G*2*N_s*sizeof(scalar));
    memcpy(_PSIxR_biased_elemwise     , PSIxR_biased_elemwise     , Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar));
    memcpy(_detJxW_full               , detJxW_full              , N_E*N_G*sizeof(scalar));
    memcpy(_detJ_full                 , detJ_full              , N_E*N_G*sizeof(scalar));
    memcpy(_serial_Uhat2GradBC        , serial_Uhat2GradBC       , M_B*D*M_G*N_s*sizeof(scalar));
    memcpy(_serial_UhCommon2GradBC    , serial_UhCommon2GradBC   , M_B*D*M_G*(N_N*M_G)*sizeof(scalar));
    
    memcpy(_ADeps , ADeps, Ne_AUG*D*sizeof(scalar));
    memcpy(_ADepsF , ADepsF, M_T*D*sizeof(scalar));
    //END PEJ EDIT
  
#elif USE_GPU
    // Allocate space on the GPU
    cudaMalloc((void**) &_map        , M_s*M_T*N_F*2*sizeof(int));            mem_counter.addToGPUCounter(M_s*M_T*N_F*2*sizeof(int));	  
    cudaMalloc((void**) &_invmap     , M_s*N_N*N_E*N_F*2*sizeof(int));	      mem_counter.addToGPUCounter(M_s*N_N*N_E*N_F*2*sizeof(int)); 	  
    cudaMalloc((void**) &_boundaryMap, M_B*sizeof(int));		      mem_counter.addToGPUCounter(M_B*sizeof(int)); 		  
    cudaMalloc((void**) &_phi        , N_G*N_s*sizeof(scalar));		      mem_counter.addToGPUCounter(N_G*N_s*sizeof(scalar)); 	  
    cudaMalloc((void**) &_phi_w      , N_G*N_s*sizeof(scalar));		      mem_counter.addToGPUCounter(N_G*N_s*sizeof(scalar)); 	  
    cudaMalloc((void**) &_dphi       , D*N_G*N_s*sizeof(scalar));	      mem_counter.addToGPUCounter(D*N_G*N_s*sizeof(scalar));  	  
    cudaMalloc((void**) &_dphi_w     , D*N_G*N_s*sizeof(scalar));	      mem_counter.addToGPUCounter(D*N_G*N_s*sizeof(scalar)); 	  
    cudaMalloc((void**) &_psi        , M_G*M_s*sizeof(scalar));		      mem_counter.addToGPUCounter(M_G*M_s*sizeof(scalar)); 	  
    cudaMalloc((void**) &_psi_w      , M_G*M_s*sizeof(scalar));		      mem_counter.addToGPUCounter(M_G*M_s*sizeof(scalar)); 	  
    //cudaMalloc((void**) &_xyz        , D*N_E*N_G*sizeof(scalar));	                                                                     
    //cudaMalloc((void**) &_xyzf       , D*M_T*M_G*sizeof(scalar));	                                                                     
    cudaMalloc((void**) &_J          , N_E*sizeof(scalar));		      mem_counter.addToGPUCounter(N_E*sizeof(scalar)); 		  
    cudaMalloc((void**) &_invJac     , N_G*D*N_E*D*sizeof(scalar));	      mem_counter.addToGPUCounter(N_G*D*N_E*D*sizeof(scalar)); 	  
    cudaMalloc((void**) &_JF         , 2*M_T*sizeof(scalar));		      mem_counter.addToGPUCounter(2*M_T*sizeof(scalar)); 		  
    cudaMalloc((void**) &_normals    , D*M_T*sizeof(scalar));		      mem_counter.addToGPUCounter(D*M_T*sizeof(scalar)); 		  
    									      	  
    cudaMalloc((void**) &_UF     , M_s*M_T*N_F*2*sizeof(scalar));	      mem_counter.addToGPUCounter(2*N_F*M_s*M_T*sizeof(scalar));      
    cudaMalloc((void**) &_Uinteg , N_G*Ne_AUG*N_F*sizeof(scalar));		      mem_counter.addToGPUCounter(N_F*N_G*Ne_AUG*sizeof(scalar)); 	      	  
    cudaMalloc((void**) &_dUinteg, D*N_G*N_E*N_F*sizeof(scalar));	      mem_counter.addToGPUCounter(D*N_G*N_E*N_F*sizeof(scalar));      	  
    cudaMalloc((void**) &_UintegF, M_G*M_T*N_F*2*sizeof(scalar));	      mem_counter.addToGPUCounter(2*N_F*M_G*M_T*sizeof(scalar));      
    									      
    cudaMalloc((void**) &_s      , N_G*N_E*N_F*sizeof(scalar));		      mem_counter.addToGPUCounter(N_G*N_E*N_F*sizeof(scalar)); 	      
    cudaMalloc((void**) &_sJ     , N_G*N_E*N_F*sizeof(scalar));		      mem_counter.addToGPUCounter(N_G*N_E*N_F*sizeof(scalar)); 	      	  
    cudaMalloc((void**) &_S      , N_s*N_E*N_F*sizeof(scalar));		      mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar)); 	      	  
    									      
    cudaMalloc((void**) &_f      , D*N_G*N_E*N_F*sizeof(scalar));	      mem_counter.addToGPUCounter(D*N_F*N_G*N_E*sizeof(scalar));      	  
    cudaMalloc((void**) &_fJ     , D*N_G*N_E*N_F*sizeof(scalar));	      mem_counter.addToGPUCounter(D*N_G*N_E*N_F*sizeof(scalar));      	  
    cudaMalloc((void**) &_F      , N_s*N_E*N_F*sizeof(scalar));		      mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar)); 	      	  
    									      
    cudaMalloc((void**) &_q      , M_G*M_T*N_F*2*sizeof(scalar));	      mem_counter.addToGPUCounter(M_G*M_T*N_F*2*sizeof(scalar));      
    cudaMalloc((void**) &_qJ     , M_G*M_T*N_F*2*sizeof(scalar));	      mem_counter.addToGPUCounter(M_G*M_T*N_F*2*sizeof(scalar));      
    cudaMalloc((void**) &_Qtcj   , M_s*M_T*N_F*2*sizeof(scalar));	      mem_counter.addToGPUCounter(M_s*M_T*N_F*2*sizeof(scalar));      
    cudaMalloc((void**) &_Q      , N_s*N_E*N_F*sizeof(scalar));  	      mem_counter.addToGPUCounter(N_s*N_E*N_F*sizeof(scalar));      

    //BR2 Edit 01/05/2016: gpu allocation for BR2 resources
    //cudaMalloc((void**) &_dUinteg_phys  , N_E*N_F*N_G*D*sizeof(scalar));      mem_counter.addToGPUCounter(N_E*N_F*N_G*D*sizeof(scalar));
    //cudaMalloc((void**) &_BR2_poly      , N_E*N_N*N_F*D*N_s*sizeof(scalar));  mem_counter.addToGPUCounter(N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //cudaMalloc((void**) &_BR2_resi      , N_E*N_N*N_F*D*N_s*sizeof(scalar));  mem_counter.addToGPUCounter(N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //cudaMalloc((void**) &_invJac_trace  , M_G*D*D*N_N*N_E*D*sizeof(scalar));    mem_counter.addToGPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //cudaMalloc((void**) &_Jac_trace     , M_G*D*D*N_N*N_E*D*sizeof(scalar));    mem_counter.addToGPUCounter(M_G*D*D*N_N*N_E*D*sizeof(scalar));
    //cudaMalloc((void**) &_phigF         , M_G*D*N_N*N_s*sizeof(scalar));        mem_counter.addToGPUCounter(M_G*D*N_N*N_s*sizeof(scalar));
    //cudaMalloc((void**) &_dphigF        , D*M_G*D*N_N*N_s*sizeof(scalar));      mem_counter.addToGPUCounter(D*M_G*D*N_N*N_s*sizeof(scalar));
    cudaMalloc((void**) &_BR2_Map       , M_T*4*sizeof(int));                 mem_counter.addToGPUCounter(M_T*4*sizeof(int));
    //cudaMalloc((void**) &_dUintegF_phys , M_T*2*N_F*M_G*D*sizeof(scalar));    mem_counter.addToGPUCounter(M_T*2*N_F*M_G*D*sizeof(scalar));
    //cudaMalloc((void**) &_M_inv         , N_E*N_s*N_s*sizeof(scalar));        mem_counter.addToGPUCounter(N_E*N_s*N_s*sizeof(scalar));
    //cudaMalloc((void**) &_weight_trace  , D*N_N*M_G*sizeof(scalar));          mem_counter.addToGPUCounter(D*N_N*M_G*sizeof(scalar));
    //Set some stuff to zero
    //cudaMemset(_dUinteg_phys , (scalar)0.0, N_E*N_F*N_G*D*sizeof(scalar));
    //cudaMemset(_BR2_poly     , (scalar)0.0, N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //cudaMemset(_BR2_resi     , (scalar)0.0, N_E*N_N*N_F*D*N_s*sizeof(scalar));
    //cudaMemset(_dUintegF_phys, (scalar)0.0, M_T*2*N_F*M_G*D*sizeof(scalar));
    //END BR2 Edit

    //PEJ 05/23/2017: gpu allocation for mixed-form. Doing this for private and public
    cudaMalloc((void**) &_UhCommon                  , M_T*M_G*N_F*sizeof(scalar));             mem_counter.addToGPUCounter(M_T*M_G*N_F*sizeof(scalar));
    cudaMalloc((void**) &_UhCommonHalf              , M_T*M_G*N_F*2*sizeof(scalar));             mem_counter.addToGPUCounter(M_T*M_G*N_F*2*sizeof(scalar));
    cudaMalloc((void**) &_PSIxR_Global              , M_T*M_G*2*N_s*sizeof(scalar));           mem_counter.addToGPUCounter(M_T*M_G*2*N_s*sizeof(scalar));
    cudaMalloc((void**) &_PSIxR_elemwise            , Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar));     mem_counter.addToGPUCounter(Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar));
    cudaMalloc((void**) &_SigFace_from_Uhat_elemwise, Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar));   mem_counter.addToGPUCounter(Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar));
    cudaMalloc((void**) &_SigFace_from_Uhat   , M_T*D*M_G*2*N_s*sizeof(scalar));         mem_counter.addToGPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
    cudaMalloc((void**) &_RecoPair                  , Ne_AUG*(D+N_N)*sizeof(int));                mem_counter.addToGPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
    cudaMalloc((void**) &_BinarySideAddress         , Ne_AUG*(D+N_N)*sizeof(int));                mem_counter.addToGPUCounter(Ne_AUG*(D+N_N)*sizeof(int));
    cudaMalloc((void**) &_sigma                     , N_E*D*N_G*N_F*sizeof(scalar));           mem_counter.addToGPUCounter(N_E*D*N_G*N_F*sizeof(scalar));
    cudaMalloc((void**) &_serial_MSigSurf           , N_E*D*N_G*N_N*M_G*sizeof(scalar));       mem_counter.addToGPUCounter(N_E * D*N_G * N_N*M_G*sizeof(scalar));
    cudaMalloc((void**) &_serial_MSigVol            , N_E*D*N_G*N_s*sizeof(scalar));           mem_counter.addToGPUCounter(N_E * D*N_G * N_s*sizeof(scalar));
    cudaMalloc((void**) &_Alt_FaceFromElem          , N_E*N_E*sizeof(int));                    mem_counter.addToGPUCounter(N_E*N_N*sizeof(int));
    //cudaMalloc((void**) &_serial_SigFace_from_DOF   , M_T*D*M_G*2*N_s*sizeof(scalar));         mem_counter.addToGPUCounter(M_T*D*M_G*2*N_s*sizeof(scalar));
    cudaMalloc((void**) &_gradCommon                , M_T*M_G*N_F*D*sizeof(scalar));            mem_counter.addToGPUCounter(M_T*M_G*N_F*D*sizeof(scalar));
    cudaMalloc((void**) &_gradCommonHalf            , M_T*M_G*N_F*D*2*sizeof(scalar));           mem_counter.addToGPUCounter(M_T*M_G*N_F*D*2*sizeof(scalar));
    cudaMalloc((void**) &_PSIxR_biased_Global       , M_T*2*M_G*2*N_s*sizeof(scalar));           mem_counter.addToGPUCounter(M_T*2*M_G*2*N_s*sizeof(scalar));
    cudaMalloc((void**) &_PSIxR_biased_elemwise     , Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar));     mem_counter.addToGPUCounter(Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar));
    cudaMalloc((void**) &_detJxW_full               , N_E*N_G*sizeof(scalar));                 mem_counter.addToGPUCounter(N_E*N_G*sizeof(scalar));
    cudaMalloc((void**) &_detJ_full               , N_E*N_G*sizeof(scalar));                 mem_counter.addToGPUCounter(N_E*N_G*sizeof(scalar));
    cudaMalloc((void**) &_serial_Uhat2GradBC        , M_B*D*M_G*N_s*sizeof(scalar));    mem_counter.addToGPUCounter(M_B*D*M_G*N_s*sizeof(scalar));
    cudaMalloc((void**) &_serial_UhCommon2GradBC    , M_B*D*M_G*(N_N*M_G)*sizeof(scalar));    mem_counter.addToGPUCounter(M_B*D*M_G*(N_N*M_G)*sizeof(scalar));

    cudaMalloc((void**) &_ADeps , Ne_AUG*D*sizeof(scalar)); mem_counter.addToGPUCounter(Ne_AUG*D*sizeof(scalar));
    cudaMalloc((void**) &_ADepsF , M_T*D*sizeof(scalar)); mem_counter.addToGPUCounter(M_T*D*sizeof(scalar));
    cudaMalloc((void**) &_DivMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &_CsMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &_LamMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &_ADepsMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &_elemCmax , Ne_AUG*sizeof(scalar)); mem_counter.addToGPUCounter(Ne_AUG*sizeof(scalar));
    cudaMalloc((void**) &GLO_DivMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &GLO_CsMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    cudaMalloc((void**) &GLO_LamMax , D*sizeof(scalar)); mem_counter.addToGPUCounter(D*sizeof(scalar));
    //end PEJ Edit

    // Set some stuff to zero
    cudaMemset(_UF     , (scalar)0.0, M_s*M_T*N_F*2*sizeof(scalar));
    cudaMemset(_Uinteg , (scalar)0.0, N_G*Ne_AUG*N_F*sizeof(scalar));
    cudaMemset(_dUinteg, (scalar)0.0, D*N_G*N_E*N_F*sizeof(scalar));
    cudaMemset(_UintegF, (scalar)0.0, M_G*M_T*N_F*2*sizeof(scalar));
    cudaMemset(_s      , (scalar)0.0, N_G*N_E*N_F*sizeof(scalar));
    cudaMemset(_sJ     , (scalar)0.0, N_G*N_E*N_F*sizeof(scalar));
    cudaMemset(_S      , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
    cudaMemset(_f      , (scalar)0.0, D*N_G*N_E*N_F*sizeof(scalar));
    cudaMemset(_fJ     , (scalar)0.0, D*N_G*N_E*N_F*sizeof(scalar));
    cudaMemset(_F      , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));
    cudaMemset(_q      , (scalar)0.0, M_G*M_T*N_F*2*sizeof(scalar));
    cudaMemset(_qJ     , (scalar)0.0, M_G*M_T*N_F*2*sizeof(scalar));
    cudaMemset(_Qtcj   , (scalar)0.0, M_s*M_T*N_F*2*sizeof(scalar));
    cudaMemset(_Q      , (scalar)0.0, N_s*N_E*N_F*sizeof(scalar));  


    //PEJ 05/23/2017: I guess I'll set UhCommon and sigma to zero
    cudaMemset(_UhCommon       , (scalar)0.0, M_T*M_G*N_F*sizeof(scalar));
    cudaMemset(_UhCommonHalf   , (scalar)0.0, M_T*M_G*N_F*2*sizeof(scalar));
    cudaMemset(_sigma          , (scalar)0.0, N_E*D*N_G*N_F*sizeof(scalar));
    cudaMemset(_gradCommon     , (scalar)0.0, M_T*M_G*N_F*D*sizeof(scalar));
    cudaMemset(_gradCommonHalf , (scalar)0.0, M_T*M_G*N_F*D*2*sizeof(scalar));
    //cudaMemset(_Uicb           , (scalar)0.0, M_T*M_G*2*N_F*sizeof(scalar));
    cudaMemset(_UicbHalf       , (scalar)0.0, M_T*M_G*2*N_F*2*sizeof(scalar));
    //PEJ 10/13/2017; Set AD gradients to zero
    cudaMemset(_NaiveGradCommonHalf, (scalar)0.0, M_T*M_G*N_F*D*2*sizeof(scalar));
    cudaMemset(_NaiveGradCommon     , (scalar)0.0, M_T*M_G*N_F*D*sizeof(scalar));
   //END PEJ Edit

    // Send the stuff to the device
    cudaMemcpy(_map        , map        , M_s*M_T*N_F*2*sizeof(int) , cudaMemcpyHostToDevice);
    cudaMemcpy(_invmap     , invmap     , M_s*N_N*N_E*N_F*2*sizeof(int) , cudaMemcpyHostToDevice);
    cudaMemcpy(_boundaryMap, m.getBoundaryMap(), M_B*sizeof(int)         , cudaMemcpyHostToDevice);
    cudaMemcpy(_phi        , phi        , N_G*N_s*sizeof(scalar)    , cudaMemcpyHostToDevice);
    cudaMemcpy(_phi_w      , phi_w      , N_G*N_s*sizeof(scalar)    , cudaMemcpyHostToDevice);
    cudaMemcpy(_dphi       , dphi       , D*N_G*N_s*sizeof(scalar)  , cudaMemcpyHostToDevice);
    cudaMemcpy(_dphi_w     , dphi_w     , D*N_G*N_s*sizeof(scalar)  , cudaMemcpyHostToDevice);
    cudaMemcpy(_psi        , psi        , M_G*M_s*sizeof(scalar)    , cudaMemcpyHostToDevice);
    cudaMemcpy(_psi_w      , psi_w      , M_G*M_s*sizeof(scalar)    , cudaMemcpyHostToDevice);
    //cudaMemcpy(_xyz        , xyz        , D*N_E*N_G*sizeof(scalar)  , cudaMemcpyHostToDevice);
    //cudaMemcpy(_xyzf       , xyzf       , D*M_T*M_G*sizeof(scalar)  , cudaMemcpyHostToDevice);
    cudaMemcpy(_J          , J          , N_E*sizeof(scalar)        , cudaMemcpyHostToDevice);
    cudaMemcpy(_invJac     , invJac     , N_G*D*N_E*D*sizeof(scalar), cudaMemcpyHostToDevice);
    cudaMemcpy(_JF         , JF         , 2*M_T*sizeof(scalar)      , cudaMemcpyHostToDevice);
    cudaMemcpy(_normals    , normals    , D*M_T*sizeof(scalar)      , cudaMemcpyHostToDevice);
    //BR2 Edit 01/05/2016: send some BR2 information to the device
    //cudaMemcpy(_invJac_trace, invJac_trace , M_G*D*D*N_N*N_E*D*sizeof(scalar) , cudaMemcpyHostToDevice);
    //cudaMemcpy(_Jac_trace   , Jac_trace    , M_G*D*D*N_N*N_E*D*sizeof(scalar) , cudaMemcpyHostToDevice);
    //cudaMemcpy(_phigF       , phigF        , M_G*D*N_N*N_s*sizeof(scalar)     , cudaMemcpyHostToDevice);
    //cudaMemcpy(_dphigF      , dphigF       , D*M_G*D*N_N*N_s*sizeof(scalar)   , cudaMemcpyHostToDevice);
    cudaMemcpy(_BR2_Map     , BR2_Map      , M_T*4*sizeof(int)              , cudaMemcpyHostToDevice);
    //cudaMemcpy(_M_inv       , M_inv        , N_E*N_s*N_s*sizeof(scalar)     , cudaMemcpyHostToDevice);
    //cudaMemcpy(_weight_trace, weight_trace , D*N_N*M_G*sizeof(scalar)       , cudaMemcpyHostToDevice);
    //END BR2 Edit

    //PEJ 05/23/2017: Send public and private to the device.
    //Since Marc doesn't do this for stuff that changes every timestep, I'm only doing it for things that have been imported
    cudaMemcpy(_PSIxR_Global               , PSIxR_Global               , M_T*M_G*2*N_s*sizeof(scalar)         , cudaMemcpyHostToDevice);
    cudaMemcpy(_PSIxR_elemwise             , PSIxR_elemwise             , Ne_AUG*(D+N_N)*M_G*N_s*sizeof(scalar)   , cudaMemcpyHostToDevice);
    cudaMemcpy(_SigFace_from_Uhat_elemwise , SigFace_from_Uhat_elemwise , Ne_AUG*(D+N_N)*D*M_G*N_s*sizeof(scalar) , cudaMemcpyHostToDevice);
    cudaMemcpy(_SigFace_from_Uhat    , SigFace_from_Uhat    , M_T*D*M_G*2*N_s*sizeof(scalar)       , cudaMemcpyHostToDevice);
    cudaMemcpy(_RecoPair                   , RecoPair                   , Ne_AUG*(D+N_N)*sizeof(int)              , cudaMemcpyHostToDevice);
    cudaMemcpy(_BinarySideAddress          , BinarySideAddress          , Ne_AUG*(D+N_N)*sizeof(int)              , cudaMemcpyHostToDevice);
    cudaMemcpy(_serial_MSigVol             , serial_MSigVol             , N_E*D*N_G*N_s*sizeof(scalar)         , cudaMemcpyHostToDevice);
    cudaMemcpy(_serial_MSigSurf            , serial_MSigSurf            , N_E*D*N_G*N_N*M_G*sizeof(scalar)     , cudaMemcpyHostToDevice);
    cudaMemcpy(_Alt_FaceFromElem           , Alt_FaceFromElem           , N_E*N_N*sizeof(int)                  , cudaMemcpyHostToDevice); 
    //cudaMemcpy(_serial_SigFace_from_DOF    , serial_SigFace_from_DOF    , M_T*D*M_G*2*N_s*sizeof(scalar)       , cudaMemcpyHostToDevice);
    cudaMemcpy(_PSIxR_biased_Global        , PSIxR_biased_Global        , M_T*2*M_G*2*N_s*sizeof(scalar)          , cudaMemcpyHostToDevice);
    cudaMemcpy(_PSIxR_biased_elemwise      , PSIxR_biased_elemwise      , Ne_AUG*(D+N_N)*2*M_G*N_s*sizeof(scalar) , cudaMemcpyHostToDevice);
    cudaMemcpy(_detJxW_full                , detJxW_full                , N_E*N_G*sizeof(scalar)               , cudaMemcpyHostToDevice);
    cudaMemcpy(_detJ_full                  , detJ_full                , N_E*N_G*sizeof(scalar)               , cudaMemcpyHostToDevice);
    cudaMemcpy(_serial_Uhat2GradBC                , serial_Uhat2GradBC                , M_B*D*M_G*N_s*sizeof(scalar)               , cudaMemcpyHostToDevice);
    cudaMemcpy(_serial_UhCommon2GradBC                , serial_UhCommon2GradBC                , M_B*D*M_G*(N_N*M_G)*sizeof(scalar)               , cudaMemcpyHostToDevice);
   
    cudaMemcpy(_ADeps, ADeps, Ne_AUG*D*sizeof(scalar) , cudaMemcpyHostToDevice);
    cudaMemcpy(_ADepsF, ADepsF, M_T*D*sizeof(scalar) , cudaMemcpyHostToDevice);
    //END PEJ EDIT
    

    
#endif

    printf("pr=%d: DGsolver check 2: Memory allocation and copying complete\n",myid);
    //To finish dgsolver declaration: Set ADepsMax
    for (int a = 0; a < D; a++)
      {
	_ADepsMax[a] = _ADeps[0*D+a];
      }
    for (int e = 0; e < _Ne_AUG; e++)
      {
	for (int a = 0; a < D; a++)
	  {
	    _ADepsMax[a] = fmax(_ADepsMax[a], _ADeps[e*D+a]);
	  }
      }
    
    printf("pr=%d, In DGsolver: ADepsMax as follows:\n",myid);
    for (int a = 0; a < D; a++)
      {
	printf("ADepsMax[%d]=%f\n", a, _ADepsMax[a]);
      }
    /*
        printf("Inspecting in DG_Solver declaration: Global PSIxR:_M_T=%d, _M_G=%d\n",_M_T,_M_G);
    for (int t = 0; t < _M_T; t++)
      {
	printf("Interface %d:\n",t);
	for (int row = 0; row < _M_G; row++)
	  {
	    printf("Quadrature node %d: ",row);
	    for (int col = 0; col < 2*_N_s; col++)
	      {
		printf("%f, ",_PSIxR_Global[t*_M_G*2*_N_s + row*2*_N_s + col]);
	      }
	    printf("\n");
	  }
      }
    */
    // Initialize some stuff for conservation calculations
    consfile = "conservation.dat";
    consf    = fopen(consfile.c_str(),"a");
    _UgC     = new scalar[N_G*N_E*N_F];  makeZero(_UgC,N_G*N_E*N_F);                   	      mem_counter.addToCPUCounter(N_G*N_E*N_F*sizeof(scalar));
    _phiC    = new scalar[N_G*N_s];      memcpy(_phiC,phi,N_G*N_s*sizeof(scalar));            mem_counter.addToCPUCounter(N_G*N_E*sizeof(scalar));
    _JC      = new scalar[N_E];          memcpy(_JC,J,N_E*sizeof(scalar));                    mem_counter.addToCPUCounter(N_E*sizeof(scalar));
    _I       = new scalar[N_F];          makeZero(_I,N_F);                                    mem_counter.addToCPUCounter(N_F*sizeof(scalar));
    _weight  = new scalar[N_G];          memcpy(_weight,weight,_N_G*sizeof(scalar));          mem_counter.addToCPUCounter(N_G*sizeof(scalar));

    //PEJ 010/01/2017: Also initialize some stuff for the tgv statistics
#ifdef TGVSTATS
    //TGVstatsf = fopen("TGVstatistics.csv","w"); //original naming system

    std::string _fname;
    _fname = "TGVstatistics.csv";
#ifdef USE_MPI
    char myidstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(myidstr, "%d", myid);
    _fname += myidstr;
#endif
    TGVstatsf = fopen(_fname.c_str(),"w");
#endif

    _EKlocal = new scalar[N_E]; makeZero(_EKlocal,N_E); mem_counter.addToCPUCounter(N_E*sizeof(scalar));
  };

  /*! Destructor */
  ~DG_SOLVER(){
    del(_map);
    del(_invmap);
    del(_boundaryMap);
    del(_phi);
    del(_phi_w);
    del(_dphi);
    del(_dphi_w);
    del(_psi);
    del(_psi_w);
    //del(_xyz);
    //del(_xyzf);
    del(_J);
    del(_invJac);
    del(_JF);
    del(_normals);
    del(_UF);
    del(_Uinteg);
    del(_dUinteg);
    del(_UintegF);
    del(_s);
    del(_sJ);
    del(_S);
    del(_f);
    del(_fJ);
    del(_F);
    del(_q);
    del(_qJ);
    del(_Qtcj);
    del(_Q);
    delete[] _UgC;
    delete[] _phiC;
    delete[] _JC;
    delete[] _I;
    delete[] _weight;
    //BR2 Edit 01/05/2016: Destroy the BR2 Additions
    //del(_phigF);
    //del(_dphigF);
    del(_BR2_Map);
    //del(_dUinteg_phys);
    //del(_dUintegF_phys);
    //del(_BR2_poly);
    //del(_BR2_resi);
    //del(_invJac_trace);
    //del(_M_inv);
    //del(_weight_trace);
    //del(_Jac_trace);
    //END BR2 Edit

    //PEJ Edit 05/23/2017
    del(_UhCommon);
    del(_UhCommonHalf);
    del(_PSIxR_Global);
    del(_PSIxR_elemwise);
    del(_SigFace_from_Uhat_elemwise);
    del(_SigFace_from_Uhat);
    del(_RecoPair);
    del(_BinarySideAddress);
    del(_sigma);
    del(_serial_MSigVol);
    del(_serial_MSigSurf);
    del(_Alt_FaceFromElem);
    //del(_serial_SigFace_from_DOF);
    del(_gradCommon);
    del(_gradCommonHalf);
    //del(_Uicb);
    del(_UicbHalf);
    del(_PSIxR_biased_Global);
    del(_PSIxR_biased_elemwise);
    //PEJ Edit 10/01/2017:
    del(_EKlocal);
    del(_detJxW_full);
    del(_detJ_full);
    //PEJ Edit 10/11/2017:
    del(_serial_Uhat2GradBC);
    del(_serial_UhCommon2GradBC);
    //PEJ Edit 10/13/2017:
    
    del(_ADeps);
    del(_ADepsF);
    del(_NaiveGradCommonHalf);
    del(_NaiveGradCommon);
    del(_DivMax);
    del(_CsMax);
    del(_LamMax);
    del(_ADepsMax);
    del(_elemCmax);
    del(GLO_DivMax);
    del(GLO_CsMax);
    del(GLO_LamMax);
    //END PEJ Edit


    fclose(consf);
#ifdef TGVSTATS
    fclose(TGVstatsf);
#endif
  };
  
  

  /*!
    \brief Main solver function
    \param[in] U solution
    \param[out] f_rk f(t,U) with DG for dU/dt = f(t,U)
    \param[in] ruk substep of RK routine, here just for debugging.
  */    
  //BR2 Edit 01/15/2016: The dg_solver function must take in RK substep index, for debug
  //void dg_solver(scalar* U, scalar* f_rk)
  // int count = 0;
  //PEJ edit 10/24/2017: Relay in the array of sensor element tags
  //PEJ 11/28/2017: Also relay in global AD approach maxima
  void dg_solver(scalar* U, scalar* f_rk, scalar ruk, int* SensorTag, scalar* rkLamMax, scalar* rkDivMax, scalar* rkCsMax, int PROC)
  {
    /*
    //PEJ Edit 01/20/2016: Look at the boundary sorting
    if (ruk == 1)
      {
	printf("BC sorting\n");
	printf("reflective interfaces:\n");
	for (int k = 0; k < _rflctiveIdx; k++)
	  {
	    //the zero is for where the reflective interfaces begin in marc's boundary map
	    printf("k=%d, t=%d\n",k,_boundaryMap[0+k]);
	  }
	printf("no-Slip interfaces:\n");
	for (int k = _noslpstrt; k < _noslpIdx; k++)
	  {
	    //the _noslpstrt is for where the no slip interfaces begin in marc's boundary map
	    printf("k=%d, t=%d\n",k,_boundaryMap[k]);
	  }
	printf("no-grad interfaces:\n");
	for (int k = _nogradstrt; k < _nogradIdx; k++)
	  {
	    //the _nogradstrt is for where the no grad interfaces begin in marc's boundary map
	    printf("k=%d, t=%d\n",k,_boundaryMap[k]);
	  }
	printf("A inflow interfaces:\n");
	for (int k = _Anflwstrt; k < _AnflwIdx; k++)
	  {
	    //the _Anflw is for where the A inflow interfaces begin in marc's boundary map
	    printf("k=%d, t=%d\n",k,_boundaryMap[k]);
	  }
      }
    //END PEJ Edit
    */
    int autocheck = 0;
    //   printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //0
    _timers.start_timer(5);
    

    //Before anything else, deal with the items that RADSINGLEFLUID needs for communication

    //  printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //2
    // collocationU: requires phi, dphi, Ustar, Uinteg, dUinteg and some sizes
    //Matrix versionL take phi@integration_points * U_solution_weights = U@integration_points
    //Example of column-major usage: leading row argument of U is N_s, so we tell
    //blassGemm that U has N_E*N_F columns; each element gets a new set of N_F columns,
    //and each column has N_s entries in it
#ifdef RADSINGLEFLUID
    _timers.start_timer(8); //Get solution at element quadrature points, flesh and ghost elements:
    blasGemm('N','N', _N_G   , _Ne_AUG*N_F, _N_s, 1, _phi,  _N_G   , U, _N_s, 0.0, _Uinteg, _N_G);
    _timers.stop_timer(8);
#endif
    //The quadrature node solution (_Uinteg) in ghost elements is necessary only
    //if using the RADSINGLEFLUID approach, where I need to know Cmax in ghost elements.
    //Otherwise, just do the collocation in flesh elements:
#ifndef RADSINGLEFLUID
    _timers.start_timer(8); //Get solution at element quadrature points, flesh and ghost elements:
    //Can't see any solution difference between using Ne_AUG and N_e here.
    //blasGemm('N','N', _N_G   , _Ne_AUG*N_F, _N_s, 1, _phi,  _N_G   , U, _N_s, 0.0, _Uinteg, _N_G);
    blasGemm('N','N', _N_G   , _N_E*N_F, _N_s, 1, _phi,  _N_G   , U, _N_s, 0.0, _Uinteg, _N_G);
    _timers.stop_timer(8);
#endif
    //   printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //3

    //Matrix version: gradphi@integration_points * U_solution_weights = gradU@integration_points (wrt reference coordinates)
    //This is Marc's procedure for calculating naive gradient in each element,
    //may be necessary for source term in some systems. I prefer to go without it for full Navier-Stokes.
#ifdef RADSINGLEFLUID
    //Do the full C-parameter sweep right here, get communication out of the way.
    _timers.start_timer(9); //Get naive solution gradient at element quadrature points:
    blasGemm('N','N', _N_G*D, _N_E*N_F, _N_s, 1, _dphi, _N_G*D, U, _N_s, 0.0, _dUinteg, _N_G*D);
    _timers.stop_timer(9);
    
    //Use naive solution gradient to get the C parameters in all flesh elements;
    //then, need to calculate some global maxima
    _timers.start_timer(43);
    LGrabCparameters(_N_G, _N_E, _Uinteg, _dUinteg, _invJac, _DivMax, _LamMax, _CsMax, _elemCmax);
    LGrabCmax_Ghosts(_N_G, _N_E, _Ne_AUG, _Uinteg, _elemCmax);
    _timers.stop_timer(43);
    /*
      printf("pr=%d, passed Cparam routines\n", PROC);
      for (int a = 0; a < D; a++)
	{
	  printf("|dir=%d: LamMax=%f, DivMax=%f, CsMax=%f| \n",a,_LamMax[a], _DivMax[a],_CsMax[a]);
	}
      */
      //Have to include the global Lam,Div,Cs communication here;
      //using values from one substep back causes failure.
      //Calculate max C parameters across all partitions
      _timers.start_timer(53);
      //Implementation C
#ifdef USE_MPI
      //This is a bottleneck. Jugging by how the AD fluxes
      //are calculated, I don't think all this communication is
      //necessary. Communication penalty could be
      //substantially reduced by altering my
      //use of global maxima in the AD fluxes as well.

      //  MPI_Barrier(MPI_COMM_WORLD); // wait until every process gets here
      //For 4 processors, code works the same with and without the barrier
      MPI_Allreduce(_DivMax, GLO_DivMax, D, MPI_SCALAR, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(_LamMax, GLO_LamMax, D, MPI_SCALAR, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(_CsMax , GLO_CsMax, D, MPI_SCALAR, MPI_MAX, MPI_COMM_WORLD);
      //Also, while the code is stopped up for Allreduce, go ahead and communicate the sensor
      //array.
      //Communcation of sensor array is necessary to achieve consistent
      //calculation of the interface AD flux for partition-shared interfaces.
      _communicator.CommunicateSensor(SensorTag); //populates ghost element sensor value
      for (int a = 0; a < D; a++)
	{
	  _DivMax[a] = GLO_DivMax[a];
	  _CsMax[a] = GLO_CsMax[a];
	  _LamMax[a] = GLO_LamMax[a];
	}
#endif //end if for MPI
      _timers.stop_timer(53);
#endif //end for RADSINGLEFLUID



    // map U onto UF: requires Map, Ustar, UF and some integers for sizes, etc
    _timers.start_timer(6); //Populate lower-dimensional projection of solution along each interface, both elements
    //  printf("Preparing LMapToFace call:\n"); fflush(stdout);
    LmapToFace(_M_s, _M_T, _N_s, _map, U, _UF); //Takes U solution vector and populates the DOF along each local interface
    _timers.stop_timer(6);

    //  printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //1
    // Apply special boundary conditions; these procedures modify ghost state variables, not the BC themselves
    
    _timers.start_timer(7);
    LrflctiveBoundary(_M_s, _rflctiveIdx, _boundaryMap, _normals, 0,          _UF); //reflective (slip wall) boundary
    LnoslipBoundary(_M_s,   _noslpcount,  _boundaryMap, _normals, _noslpstrt, _UF); //the no-slip (viscous wall) boundary
    LAnflwBoundary(_M_s,    _Anflwcount,  _boundaryMap, _normals, _Anflwstrt, _UF); //A inflow (simple supersonic) inflow
    LKJetBoundary(_M_s,    _KJetcount,  _boundaryMap, _normals, _KJetstrt, _UF); //KJet inflow 
    LHomoBoundary(_M_s,    _Homocount,  _boundaryMap, _normals, _Homostrt, _UF); //Homogeneous (U=2) BC 
    LSubOutBoundary(_M_s,    _SubOutcount,  _boundaryMap, _normals, _SubOutstrt, _UF); //Subsonic Outflow BC 
    _timers.stop_timer(7);
    
    /*
    printf("pr=%d: UF distribution along boundary interfaces:\n",PROC);
    for (int tB = 0; tB < _M_B; tB++)
      {
	int t = _boundaryMap[tB];
	printf("\tpr=%d: boundary face %d, global Interface %d:\n",PROC, tB, t);
	    for (int k = 0; k < _M_s; k++)
	      {
#ifdef SCALARAD
		scalar diff = fabs(_UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k]-_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k]);
		if (fabs(diff) < 0.000001)
		  {
		    printf("\t\t face node %d: rhoL=%f, rhoR=%f\n", k, _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k]);
		  }
		else
		  {
		    printf("\t\t face node %d: rhoL=%f, rhoR=%f\t\t\tDISCREPANCY=%f\n", k, _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k],diff);
		  }
		
#endif
#ifdef RADSINGLEFLUID
		printf("\t\t face node %d: (rho,mox,Et,C)_L=(%f,%f,%f,%f), (rhoR,mox,Et,C)_R=(%f,%f,%f,%f)\n", k,
		       _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 1*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 2*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 3*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 1*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 2*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 3*2*_M_s + 1*_M_s + k]);
#endif
	      } 
      }
    */

 
    // printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //4
	
#ifdef SIMPLE
    // collocationUF: requires psi, UF, UintegF and some sizes
    //Matrix Version (ignore non-supported Lagrange functions): Take phi@face_integration_points * U_face_solution_weights = U@face_integration_points
    _timers.start_timer(10); //get solution at interface quadrature points:
    blasGemm('N','N', _M_G, _M_T*N_F*2, _M_s, 1, _psi, _M_G, _UF, _M_s, 0.0, _UintegF, _M_G);
    _timers.stop_timer(10);
    /*
    printf("Boundary solution at quadrature points, SIMPLE\n");
    for (int tB = 0; tB < _M_B; tB++)
      {
	int t = _boundaryMap[tB];
	printf("\tBoundary face %d, global face index %d:\n",tB,t);
	for (int g = 0; g < _M_G; g++)
	  {
	    scalar diff = fabs(_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g]-_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
	    if (diff < 0.000001)
	      {
		printf("\t\t quadrature node %d: rhoL=%f, rhoR=%f\n", g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
	      }
	    else
	      {
		printf("\t\t quadrature node %d: rhoL=%f, rhoR=%f\t\t\tDISCREPANCY=%f\n", g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g],diff);
	      }
	  }
      }
    */
#endif
    //printf("PREPARE FOR AWESOME\n");
#ifdef ICBN //counterpart to SIMPLE.
    //Perform ICB reconstruction for interface solution to Riemann solver. Store ICB solution as UintegF.
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //5

    //LUicbHalf_to_Uicb(_M_T, _M_G, _UicbHalf, _UintegF);
    //This has to include ghost elements regardless of whether I am using RADSINGLEFLUID,
    //because partition-boundary interface needs info from elements on both sides.

    /*
      //STARSCREAM implementation:
    _timers.start_timer(41);
    LUhat_to_UicbHalf(_N_E, _Ne_AUG, _N_N, _M_T, _M_G, _N_s, _RecoPair, _BinarySideAddress, _PSIxR_biased_elemwise, U, _UicbHalf);
    _timers.stop_timer(41);
    
    _timers.start_timer(42);
    LUicbHalf_to_Uicb(_M_T, _M_G, _UicbHalf, _UintegF);
    _timers.stop_timer(42);
    */

    //The new approach: Do ICB approximation in 1 step
    _timers.start_timer(41);
    LUhat_to_UicbDirect(_M_T, _M_G, _N_s, _BR2_Map, _PSIxR_biased_Global, U, _UintegF);
    _timers.stop_timer(41);

    //Also, populate flux properly along boundary interfaces. Uicb routine does it incorrectly, I remedy error here.
    LBuildBC_UintegF(_M_B, _M_G, _M_s, _boundaryMap, _psi, _UF, _UintegF);
    //blasGemm('N','N', _M_G, _M_T*N_F*2, _M_s, 1, _psi, _M_G, _UF, _M_s, 0.0, _UintegF, _M_G);
    /*
    printf("Just executed BuildBC_UintegF. Here is the new distribution along boundary faces:\n");
    for (int tB = 0; tB < _M_B; tB++)
      {
	int t = _boundaryMap[tB];
	printf("\tBoundary face %d, global face index %d:\n",tB,t);
	for (int g = 0; g < _M_G; g++)
	  {
	    scalar diff = fabs(_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g]-_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
	    if (diff < 0.000001)
	      {
		printf("\t\t quadrature node %d: rhoL=%f, rhoR=%f\n", g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
	      }
	    else
	      {
		printf("\t\t quadrature node %d: rhoL=%f, rhoR=%f\t\t\tDISCREPANCY=%f\n", g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g],diff);
	      }
	  }
      }
    */
#endif
    //printf("Exited the Uicb population procedures, back in dg_solver.h\n");
    
    //return 0;
    /*
    if (ruk == 0)
      {
	printf("pr=%d UF distribution:\n",PROC);
	for (int t = 0; t < _M_T; t++)
	  {
	    printf("\tpr=%d, Interface %d:\n",PROC,t);
	    for (int k = 0; k < _M_s; k++)
	      {
#ifdef SCALARAD
		scalar diff = fabs(_UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k]-_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k]);
		if (fabs(diff) < 0.000001)
		  {
		    printf("\t\t pr=%d, face node %d: rhoL=%f, rhoR=%f\n", PROC, k, _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k]);
		  }
		else
		  {
		    printf("\t\t pr=%d, face node %d: rhoL=%f, rhoR=%f\t\t\tDISCREPANCY=%f\n", PROC, k, _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],_UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k],diff);
		  }
		
#endif
#ifdef RADSINGLEFLUID
		printf("\t\t face node %d: (rho,mox,Et,C)_L=(%f,%f,%f,%f), (rhoR,mox,Et,C)_R=(%f,%f,%f,%f)\n", k,
		       _UF[t*N_F*2*_M_s + 0*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 1*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 2*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 3*2*_M_s + 0*_M_s + k],
		       _UF[t*N_F*2*_M_s + 0*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 1*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 2*2*_M_s + 1*_M_s + k],
		       _UF[t*N_F*2*_M_s + 3*2*_M_s + 1*_M_s + k]);
#endif
	      } 
	  }
	printf("pr=%d, _UintegF distribution:\n",PROC);
	for (int t = 0; t < _M_T; t++)
	  {
	    printf("\tpr=%d, Interface %d:\n",PROC,t);
	    for (int g = 0; g < _M_G; g++)
	      {
		scalar diff = fabs(_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g]-_UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
		if (diff < 0.000001)
		  {
		    printf("\t\t pr=%d, quadrature node %d: rhoL=%f, rhoR=%f\n", PROC, g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g]);
		  }
		else
		  {
		    printf("\t\t pr=%d, quadrature node %d: rhoL=%f, rhoR=%f\t\t\tDISCREPANCY=%f\n", PROC, g, _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 0*_M_G + g], _UintegF[t*N_F*2*_M_G + 0*2*_M_G + 1*_M_G + g],diff);
		  }
	      }
	  }
      }
    */

    // Physics
    //Send the U solution and gradient at interior integration points, receive the source and interior flux terms
    _timers.start_timer(11);
    Levaluate_sf(_N_G, _N_E, _s, _f, _Uinteg, _dUinteg, _invJac);//, _xyz);
    _timers.stop_timer(11);
    // printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //7
	
    //Send the U trace from each side of each quad point, each interface, receive F (dot) n on each side
    _timers.start_timer(12);
    Levaluate_q(_M_G, _M_T, _q, _UintegF, _normals);//, _xyzf);
    _timers.stop_timer(12);
    // printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; //8
	
    //CGR + Viscous Physics Subroutine:
    //BR2 Edit 01/06/2016: declaring _viscous to be true(1)/false(0) here. In future, specify it with deck
    //PEJ 05/18/2017: Adjusting how gradients are computed, but keeping
    //viscous flux functions from my earlier BR2 attempt.
    //  int _viscous = 1;
#ifdef NOVIS
    // _viscous = 0;
#endif
#ifndef NOVIS //novis is NOT defined, that means we need viscous routine
    {
	//	printf("Entered viscous loop\n");

	//Old approach, which may or may not have been an accurate implementation of BR2:
	/*
	//Obtain the solution jump at each interface quadrature point, integrate for BR2 residual
	Jump_Residual(_M_T, _N_s,  _M_G, _N_N, _BR2_Map, _UintegF, _phigF, _normals, _Jac_trace, _JF, _weight_trace,  _BR2_resi, ruk);
	
	//Solve BR2 residual, using inverse mass matrix, to get the correction DOF
	Lift_Solve(_N_E,  _N_N, _N_s, _M_inv, _BR2_resi, _BR2_poly, ruk); 
        
	//Use DG DOF alongside correction DOF to populate gradient over each element interior
	populate_CellGrad(_N_E, _N_G, _N_s, _N_N, U, _invJac, _dphi, _phi, _BR2_poly, _dUinteg_phys, ruk); 
        
	//USe DG DOF alongside correction DOF to populate gradient over each interface; this is where BR1/BR3 modification can be applied
	populate_FaceGrad(_M_T,  _M_G, _N_s, _N_N, _BR2_Map, U, _invJac_trace, _phigF, _dphigF, _BR2_poly, _dUintegF_phys, ruk);
	*/
	//05/18/2017 New approach:
	//Step 1: Populate the common solution value tilde{U} on each interface quadr point.
	//This needs to be done directly from DG DOF.
	//Specifics: Uhat=U is in column major form, NE*NF columns, K rows, with NF columns per element. UhCommon is a NF columns, MT*MG rows; each interface gets MG rows.

	/*
	if (ruk == 0)
	  {
	    printf("Global PSIxR:_M_T=%d, _M_G=%d\n",_M_T,_M_G);
	    for (int t = 0; t < _M_T; t++)
	      {
		printf("Interface %d:\n",t);
		for (int row = 0; row < _M_G; row++)
		  {
		    printf("Quadrature node %d: ",row);
		    for (int col = 0; col < 2*_N_s; col++)
		      {
			printf("%f, ",_PSIxR_Global[t*_M_G*2*_N_s + row*2*_N_s + col]);
		      }
		    printf("\n");
		  }
	      }
	  }
	*/

	/*
	Cool fact: inside this routine, calling printf causes segfault
	if you're trying to print an array value and using GPU.
	printf("The DOF vector:\n");
	for (int e = 0; e < _N_E; e++)
	  {
	    for (int k = 0; k < _N_s; k++)
	      {
		for (int f = 0; f < N_F; f++)
		  {
		    printf("DOF(e=%d,f=%d,k=%d) = %f\n", e, f, k, U[(e*N_F+f)*_N_s+k]);
		  }
	      }
	  }
	*/
	/*
	  LUhat_to_UhCommon populates each interface's common solution value
	  from the DOF of 2 surrounding elements. LUhat_to_CUGUF goes
	  one step further and grabs the interface gradient also
	printf("\nEntered DG Solver, executing LUhat_to_UhCommon\n");
	LUhat_to_UhCommon(_M_T, _M_G, _N_s,_BR2_Map, _PSIxR_Global, U, _UhCommon); 
	printf("Finished with LUhat_to_UhCommon\n");
	  */
	/*
	printf("RecoPair:\n");
	for (int e = 0; e < _N_E; e++)
	  {
	    for (int s = 0; s < D+_N_N; s++)
	      {
		printf("RecoPair(%d,%d) = %d\n", e, s, _RecoPair[e*(D+_N_N)+s]);
	      }
	  }
	*/
	//zero UhCommonHalf
	/*
	for (int j = 0; j < _M_T*_M_G*N_F; j++){
	  _UhCommonHalf[j] = 0; }
	*/

      /*
      //Previous working approach, hopefully dies with STARSCREAM
	_timers.start_timer(33);
	LUhat_to_UhCommonHalf(_N_E, _Ne_AUG, _N_N, _M_T, _M_G, _N_s, _RecoPair, _BinarySideAddress, _PSIxR_elemwise, U, _UhCommonHalf);
	_timers.stop_timer(33);

	_timers.start_timer(34);
	LUhCommonHalf_to_UhCommon(_M_T, _M_G, _UhCommonHalf, _UhCommon);
	_timers.stop_timer(34);
      */
      /*
      printf("SigFace_from_Uhat in DGsolver:\n");
      for (int t = 0; t < _M_T; t++){
	for (int g = 0; g < _M_G; g++){
	  for (int a = 0; a < D; a++){
	    for (int k = 0; k < 2*_N_s; k++)
	      {
		printf("\tentry(t=%d,g=%d,a=%d,k=%d)=%f\n", t, g, a, k, _SigFace_from_Uhat[t*(_M_G*D*2*_N_s) + g*(D*2*_N_s) + a*(2*_N_s) + k]);
	      }}}}
      */
      //The combined UhCommon and GradCommon population:
      _timers.start_timer(33);
      LUhat_to_CUGUF(_M_T, _M_G, _N_s, _BR2_Map, _PSIxR_Global, _SigFace_from_Uhat, U, _UhCommon, _gradCommon);
      _timers.stop_timer(33);

	//Next: Look at boundary distribution of UintegF and use it to correct UhCommon along boundaries
      LCorrectUhCommonBC(_M_T, _M_G, _M_B, _boundaryMap, _UintegF, _UhCommon);


	//LUhat_to_CUGUF(_M_T, _M_G, _N_s, _BR2_Map, _PSIxR_Global, _serial_SigFace_from_DOF, U, _UhCommon, _gradCommon);
	/*
	if (ruk == 0){
	  printf("processor %d: The Uh_common Distribution:\n",PROC);
	  for (int t = 0; t < _M_T; t++)
	    {
	      printf("pr=%d, tGlo=%d:\n",PROC,t);
	      for (int g = 0; g < _M_G; g++)
		{
		  printf("g=%d: N_F entries are    (",g);
		  for (int f = 0; f < N_F; f++)
		    {
		      printf("%f, ",_UhCommon[t*_M_G*N_F + g*N_F + f]);
		    }
		  printf(")\n");
		}
	      printf("\n");
	    }
	}
	*/
	/*
	  //Previous working approach, hopefully dies with STARTSCREAM
	_timers.start_timer(35);
	LUhat_to_GradCommonHalf(_N_E, _Ne_AUG, _N_N, _M_T, _M_G, _N_s, _RecoPair, _BinarySideAddress, _SigFace_from_Uhat_elemwise, U, _gradCommonHalf);
	_timers.stop_timer(35);

	_timers.start_timer(36);
	LGradCommonHalf_to_GradCommon(_M_T, _M_G, _gradCommonHalf, _gradCommon);
	_timers.stop_timer(36);
	*/

	//Get the gradient along boundary interfaces (use sigma trace)
	//	scalar* _serial_Uhat2GradBC = new scalar[1];
	//	scalar* _serial_UhCommon2GradBC = new scalar[1]; 
	LUhat_to_GradCommonBC(_M_B, _M_G, _N_s, _N_N, _boundaryMap, _BR2_Map, _serial_Uhat2GradBC, U, _RecoPair, _Alt_FaceFromElem, _UhCommon, _serial_UhCommon2GradBC, _gradCommon);

	//LUhat_to_GradCommon_v2(_N_E, _N_N, _M_T, _M_G, _N_s, _RecoPair, _SigFace_from_Uhat_elemwise, U, _gradCommon);
	//LUhat_to_CUGUF(_M_T, _M_G, _N_s, _BR2_Map, _PSIxR_Global, _serial_SigFace_from_DOF, U, _UhCommon, _gradCommon);

	/*
	if (ruk == 0){
	  printf("processor=%d: The grad_common Distribution:\n",PROC);
	  for (int t = 0; t < _M_T; t++)
	    {
	      printf("pr=%d, tGlo=%d, omA=%d, omB=%d:\n",PROC,t,_BR2_Map[t*4+0], _BR2_Map[t*4+2]);
	      // if (_BR2_Map[t*4+2] >= _N_E)
		{
	      for (int g = 0; g < _M_G; g++)
		{
		  for (int a = 0; a < D; a++)
		    {
		      printf("g=%d,a=%d: N_F entries are    (",g,a);
		      for (int f = 0; f < N_F; f++)
			{
			  printf("%f, ",_gradCommon[t*_M_G*N_F*D + g*N_F*D + f*D + a]);
			}
		      printf(")\n");
		    }
		}
		}
	      printf("\n");
	    }
	}
	*/

	//return 0;
	//	printf("\nEntered DF solver, executing LUhat_to_CUGUF\n");
	//LUhat_to_CUGUF(_M_T, _M_G, _N_s, _BR2_Map, _PSIxR_Global, _serial_SigFace_from_DOF, U, _UhCommon, _gradCommon);
	//	printf("Finished with LUhat_to_CUGUF\n");
	/*
	if (ruk == 0){
	  printf("The Uh_common Distribution:\n");
	  for (int t = 0; t < _M_T; t++)
	    {
	      printf("tGlo=%d:\n",t);
	      for (int g = 0; g < _M_G; g++)
		{
		  printf("g=%d: N_F entries are \n (",g);
		  for (int f = 0; f < N_F; f++)
		    {
		      printf("%f, ",_UhCommon[t*_M_G*N_F + g*N_F + f]);
		    }
		  printf(")\n");
		}
	      printf("\n");
	    }
	*/
	/*
	  printf("The gradCommon Distribution:\n");
	  for (int t = 0; t < _M_T; t++)
	    {
	      printf("tGlo=%d:\n",t);
	      for (int g = 0; g < _M_G; g++)
		{
		  printf("g=%d: drho_dx entry is \n (",g);
		  for (int f = 0; f < 1; f++)
		    {
		      printf("%f, ",_gradCommon[t*_M_G*N_F*D + g*N_F*D + f*D + 0]);
		    }
		  printf(")\n");
		}
	      printf("\n");
	    }
	}
	*/



	//Step 2: Solve weak form for element interior gradients (sigma), store as _dUinteg_phys
	//Inputs to this step are DG DOF and the interface tilde{U} values
	//Step 2a: The volume contribution
	//	printf("About to call LUHat_to_Sigma\n");
	_timers.start_timer(37);
	LUhat_to_Sigma(_N_E, _N_G, _N_s, _serial_MSigVol, U, _sigma);
	_timers.stop_timer(37);
	//	printf("Finished with LUhat_to_Sigma\n");

	//	printf("About to callLUhCommon_to_Sigma\n");
	_timers.start_timer(38);
	LUhCommon_to_Sigma(_N_E, _N_G, _N_N, _M_G, _serial_MSigSurf, _Alt_FaceFromElem, _UhCommon, _sigma);
	_timers.stop_timer(38);
	//	printf("Finished with LUhCommon_to_Sigma\n");
	/*
	if (ruk == 0)
	  {
	    printf("processor=%d, The sigma(rho) approximation, after volume and surface contribution\n",PROC);
	    for (int e = 0; e < _N_E; e++)
	      {
		printf("\tpr=%d, element %d:\n", PROC, e);
		for (int a = 0; a < D; a++)
		  {
		    for (int g = 0; g < _N_G; g++)
		      {
			printf("\t\tnode %d: drho_dX(a=%d) = %f\n", g, a, _sigma[e*D*_N_G*N_F + a*_N_G*N_F + g*N_F + 0]);
		      }
		  }
		printf("\n");
	      }
	  }
	*/
	
	/*
	  gradients look good for sine wave in x test. Next, need to deal with the physics package
	  and put in scalar advection-diffusion so I can be sure that things are working properly.
	 */

	//Step 3: Solve for interface gradients, store as _dUintegF_phys.
	//Inputs here are just the DG DOF, I must grab interface gradients efficiently


	//Apply Naumann condition, normal derivative=0 along boundary interfaces
	//For the no-gradient boundaries, the tangential gradient must be supplied by the boundary element;
	// this supply is automatically attended to in Populate_FaceGrad
	//This function will modify dUintegF_phys, so it cannot be in same spot as other boundary routines
	//The necessary values of UintegF are populated by marc's routine, require no special treatment
	//This function sets the solution gradient zero normal to the boundary, an act for which I do not know
	//and physical justification. It is outlawed for now.

	//LnogradBoundary(_M_s,   _nogradcount,   _boundaryMap, _normals, _nogradstrt, _dUintegF_phys); 
      
	//Viscous flux subroutines: requure _dUinteg_phys= sigma for interior gradient 
	//Viscous flux subroutines: requre _UintegF and _dUintegF_phys for interface gradientes


        //Apply viscous flux component to cell interior fluxes

	_timers.start_timer(39);
	//	Levaluate_sf_vis(_N_G, _N_E, _Mu, _Pr, _KTemp, _f, _Uinteg, _dUinteg_phys); //Arbitrary viscous flux subroutine
	Levaluate_sf_vis(_N_G, _N_E, PROC, _f, _Uinteg, _sigma); //Arbitrary viscous flux subroutine
	_timers.stop_timer(39);
	
        //Average viscous flux component along each interface, then apply the reuslt to interface fluxes
	_timers.start_timer(40);
	//Levaluate_q_vis(_M_G, _M_T, _Mu, _Pr, _KTemp,  _q, _UintegF, _dUintegF_phys, _normals);
	Levaluate_q_vis(_M_G, _M_T, PROC, _q, _UhCommon, _gradCommon, _normals);
	_timers.stop_timer(40);

      } //end case structure for viscous
#endif

      
    //END BR2+Navier-Stokes Section

    //PEJ 10/15/2017: Reisner AD subroutines.
#ifdef RADSINGLEFLUID
        {
	  //	  printf("pr=%d, entering RADSINGLEGFLUID section of dgsolver\n", PROC);
      //Identify some global maxima: (looking for DivMax, LamMax, and CsMax) and also max C per element.
      //Currently using naive gradient for this estimation subroutine because I don't think gradient accuracy is important here
     


      /*
       printf("pr=%d, calculated global Cparam routines\n", PROC);
       for (int a = 0; a < D; a++)
	{
	  printf("|dir=%d: GLO_LamMax=%f, GLO_DivMax=%f, GLO_CsMax=%f| \n",a,GLO_LamMax[a], GLO_DivMax[a],GLO_CsMax[a]);
	}
      */
      /*
      printf("pr=%d:The calculated Cparam values on flesh elements:\n",PROC);
      for (int e = 0; e < _N_E; e++)
	{
	  printf("pr=%d, e=%d: elemCmax=%f\n ", PROC, e, _elemCmax[e]);
	}
      printf("pr=%d: The calculated elemCmax in ghosts:\n",PROC);
      for (int e = _N_E; e < _Ne_AUG; e++)
	{
	  printf("pr=%d, e=%d: elemCmax=%f\n", PROC, e, _elemCmax[e]);
	}
      */
      //Next: Get common interface solution on all AD interfaces
      _timers.start_timer(48);
      //BEGIN YOUR CONQUEST HERE
      LUhat_to_UhCommonHalf_AD(_N_E, _Ne_AUG, _N_N, _M_T, _M_G, _N_s, _Cthresh, _elemCmax, SensorTag, _RecoPair, _BinarySideAddress, _PSIxR_elemwise, U, _UhCommonHalf); //done
      _timers.stop_timer(48);
      
      _timers.start_timer(49);
      LUhCommonHalf_to_UhCommon_AD(_M_T, _M_G, _Cthresh, _elemCmax, SensorTag,  _BR2_Map,_UhCommonHalf, _UhCommon); //no need to change
      _timers.stop_timer(49);

      //Some boundary considerations: No sensor input because the boundary interface
      //population is small compared to full interface set, and the sensor
      //is probably triggered anyway at the boundary
      //Next: Look at boundary distribution of UintegF and use it to correct UhCommon along boundaries
      LCorrectUhCommonBC(_M_T, _M_G, _M_B, _boundaryMap, _UintegF, _UhCommon); //no change
      
      
      //Get the gradient along boundary interfaces (use sigma trace)
      //	scalar* _serial_Uhat2GradBC = new scalar[1];
      //	scalar* _serial_UhCommon2GradBC = new scalar[1]; 
      //	LUhat_to_GradCommonBC(_M_B, _M_G, _N_s, _N_N, _boundaryMap, _BR2_Map, _serial_Uhat2GradBC, U, _RecoPair, _Alt_FaceFromElem, _UhCommon, _serial_UhCommon2GradBC, _NaivegradCommon);
      
      //Calculate mixed-form gradient over element interiors, wherever AD operation is needed
      _timers.start_timer(50);
      LUhat_to_Sigma_AD(_N_E, _N_G, _N_s, _Cthresh, _elemCmax, SensorTag, _serial_MSigVol, U, _sigma); //no change, only relevant over flesh elements
      _timers.stop_timer(50);
      
      _timers.start_timer(51);
      LUhCommon_to_Sigma_AD(_N_E, _N_G, _N_N, _M_G, _Cthresh, _elemCmax, SensorTag, _BR2_Map, _serial_MSigSurf, _Alt_FaceFromElem, _UhCommon, _sigma); //no change here
      _timers.stop_timer(51);

      //Calculate interface gradient wherever AD operation is needed.
      //NOt actually naive, that is a mislabel. It becomes naive if you take Chi_AD=0 in main.
      _timers.start_timer(44);
      LUhat_to_GradCommonHalf_AD(_N_E, _Ne_AUG, _N_N, _M_T, _M_G, _N_s, _Cthresh, _elemCmax, SensorTag, _RecoPair, _BinarySideAddress, _SigFace_from_Uhat_elemwise, U, _NaiveGradCommonHalf); //done
      _timers.stop_timer(44);
      
      _timers.start_timer(45);
      LGradCommonHalf_to_GradCommon_AD(_M_T, _M_G, _Cthresh, _elemCmax, SensorTag,  _BR2_Map, _NaiveGradCommonHalf, _NaiveGradCommon); //no change
      _timers.stop_timer(45);
      
      //Get the gradient along boundary interfaces (use sigma trace)
      LUhat_to_GradCommonBC(_M_B, _M_G, _N_s, _N_N, _boundaryMap, _BR2_Map, _serial_Uhat2GradBC, U, _RecoPair, _Alt_FaceFromElem, _UhCommon, _serial_UhCommon2GradBC, _NaiveGradCommon); //no change
      /*
      printf("AD The Uh_common Distribution:\n");
      for (int t = 0; t < _M_T; t++)
	{
	  printf("tGlo=%d:\n",t);
	  for (int g = 0; g < _M_G; g++)
	    {
	      printf("g=%d: N_F entries are    (",g);
	      for (int f = 0; f < N_F; f++)
		{
		  printf("%f, ",_UhCommon[t*_M_G*N_F + g*N_F + f]);
		}
	      printf(")\n");
	    }
	  printf("\n");
	}
      printf("The gradCommon Distribution:\n");
	  for (int t = 0; t < _M_T; t++)
	    {
	      printf("tGlo=%d:\n",t);
	      for (int g = 0; g < _M_G; g++)
		{
		  printf("g=%d: drho_dx entry is \n (",g);
		  for (int f = 0; f < N_F; f++)
		    {
		      printf("%f, ",_NaiveGradCommon[t*_M_G*N_F*D + g*N_F*D + f*D + 0]);
		    }
		  printf(")\n");
		}
	      printf("\n");
	    }
	   printf("The sigma_x approximation, after volume and surface contribution\n");
	    for (int e = 0; e < _N_E; e++)
	      {
		printf("\telement %d:\n", e);
		for (int a = 0; a < 1; a++)
		  {
		    for (int g = 0; g < _N_G; g++)
		      {
			printf("\t\tnode %d: ", g);
			for (int fc = 0; fc < N_F; fc++)
			  {
			    printf("f(U%d)_dx=%f, ",fc, _sigma[e*D*_N_G*N_F + a*_N_G*N_F + g*N_F + fc]);
			  }
			printf("\n");
			    
		      }
		  }
		printf("\n");
	      }
      */
      int _order = 0; //not actually used for anything right now, but fed as an argument, might use later.

      //gradient is now known at all quadrature points, evaluate AD flux.
      //CRITICAL: Levaluate_sf_rad must import proper gradient; sigma is wrt phys coordinates, dUinteg is wrt ref coordinates/
      //printf("Calling Levaluate_sf_rad\n"); fflush(stdout);
      _timers.start_timer(46);
      //Levaluate_sf_rad(_N_G, _N_E, _order, _DivMax, _Beta_S, _Mew_S, _CsMax, _LamMax, _epsGen, _ADeps, _Cthresh, _elemCmax, SensorTag, _invJac, _Uinteg, _dUinteg, _s, _f);
      Levaluate_sf_rad(_N_G, _N_E, _order, _DivMax, _Beta_S, _Mew_S, _CsMax, _LamMax, _epsGen, _ADeps, _Cthresh, _elemCmax, SensorTag, _invJac, _Uinteg, _sigma, _s, _f); //no change
      _timers.stop_timer(46);
      
      //printf("Calling Levaluate_q_rad\n"); fflush(stdout);
      _timers.start_timer(47);
      Levaluate_q_rad(_M_G, _M_T, _q, _UintegF, _NaiveGradCommon, _DivMax, _Beta_S, _Mew_S, _CsMax, _LamMax, _epsGen, _ADepsF, _Cthresh, _elemCmax, SensorTag, _BR2_Map, _normals); //no change
      _timers.stop_timer(47);
      

    }
#endif
    //exit(1);
    // redistribute_sf: requires J, invJac, s, f, phi_w, dphi_w, sJ, fJ, S, F
    //Take source and interior flux output, multiply by Jacobian to prepare for integration
    //This function also multiplies the interior flux by inverse Jacobian matrix, such that differentiation
    //of phi can be based on reference element (compromise between memory usage and flops)
    _timers.start_timer(13);
    //Lredistribute_sf(_N_G, _N_E, _sJ, _fJ, _s, _f, _J, _invJac);
    Lredistribute_sf(_N_G, _N_E, _sJ, _fJ, _s, _f, _J, _detJ_full, _invJac);
    _timers.stop_timer(13);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    // matrix-matrix multiply for sf
    //Perform integration of source distribution and populate Source residual
    _timers.start_timer(14);
    blasGemm('T','N', _N_s, _N_E*N_F, _N_G   , 1, _phi_w , _N_G   , _sJ, _N_G  , 0.0, _S, _N_s);
    _timers.stop_timer(14);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    //Perform integration of interiorflux * gradphi(ref elements) to get interior residual
    _timers.start_timer(15);
    blasGemm('T','N', _N_s, _N_E*N_F, _N_G*D, 1, _dphi_w, _N_G*D, _fJ, _N_G*D, 0.0, _F, _N_s);
    _timers.stop_timer(15);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    // redistribute_q: requires JF, q, qJ, psi_w, Qtcj,
    //At each interface, multiply surface flux outflow by JF, which I think is the Jacobian determinant on each side of interface
    
    _timers.start_timer(16);
    Lredistribute_q(_M_G, _M_T, _qJ, _q, _JF);
    _timers.stop_timer(16);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    // matrix-matrix multiply for q
    //Take integration of interface flux against supported test functions for Qtcj, which contains
    //surface residual, but in interface-oriented format
    _timers.start_timer(17);
    blasGemm('T','N', _M_s, _M_T*N_F*2, _M_G, 1, _psi_w , _M_G, _qJ, _M_G, 0.0, _Qtcj, _M_s);
    _timers.stop_timer(17);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    // map_q: requires map, Qtcj, Q (might want to do this in the previous step)
    //Takes Qtcj residual contributions, moves them to element orientation
    //I think that for a single testing basis function on the element,
    //there can be multiple interface contributions (think corner nodes)
    _timers.start_timer(18);
    LmapToElement(_N_s, _N_E, _M_s, _N_N, _invmap, _Q, _Qtcj);
    _timers.stop_timer(18);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 
	
    // Make f_rk = S+F+Q
    _timers.start_timer(19);
    LaddSFQ(_N_s, _N_E, f_rk, _S, _F, _Q);
    _timers.stop_timer(19);
    //printf("DGSOLVER CHECK %d\n",autocheck); fflush(stdout); autocheck++; 

    _timers.stop_timer(5);
  }; // end solver function


  void conservation(scalar* U, double time){
    /*!
      \brief Function to calculate and output conservation of certain quantities
      \param[in] U solution to evaluate (Lagrange nodal), on the host (CPU)
      \param[in] time time step
    */
    
    // Collocate the solution to the integration points
    hostblasGemm('N','N', _N_G, _N_E*N_F, _N_s, 1, _phiC, _N_G, U, _N_s, 0.0, _UgC, _N_G);
    
    // Take the cell average of the solution
    makeZero(_I, N_F);
    for(int fc = 0; fc < N_F; fc++){
      for(int e = 0; e < _N_E; e++){
    	for(int g = 0; g < _N_G; g++){
    	  _I[fc] += _UgC[(e*N_F+fc)*_N_G+g]*_JC[e]*_weight[g];
    	}
      }
    }
    // write to file
    fprintf(consf,"%20.16E\t", time); for(int fc = 0; fc < N_F; fc++) fprintf(consf,"%20.16E\t", _I[fc]); fprintf(consf,"\n");

  };//end conservation function

  scalar* getPhiW()const {/*! Return phi*/return _phi_w;}

  //I need a couple routines to send C equation information to the rk routine
  scalar* get_ADepsMax() {/*! Return ADeps*/ return _ADepsMax;}
  scalar* get_DivMax(){/*! return max super-divergence*/ return _DivMax;}
  scalar get_Beta_S() {/* Return AD strength*/ return _Beta_S;}
  scalar get_Mew_S() {/* Return AD spread factor*/ return _Mew_S;}

  //More routines to send C equation information to rk routine
  scalar* getLamMax(){/*! return max wavespeed (each direction)*/ return _LamMax;}
  scalar* getCsMax(){/*!\brief return max Cs value (each direction)*/ return _CsMax;}

  void BulkyAD_RelayParams(scalar* rkUhat, scalar* rkLamMax, scalar* rkDivMax, scalar* rkCsMax)
  {
    /*!\brief calculate global maxima for AD approach starting from Uhat
      \param[in] rkUhat the DOF vector stored in RK routine, relayed in
      \param[out] rkLamMax partition's max wavespeed in each direction
      \param[out] rkDivMax partition's max divergence component in each direction
      \param[out] rkCsMax partition's max Cs value in each direction
     */
    //This routine only needs to run on flesh elements, as ghost
    //elements will be properly accounted for by their home partitions.
    //Step 1: Calculate solution at quadrature points
    blasGemm('N','N', _N_G   , _N_E*N_F, _N_s, 1, _phi,  _N_G   , rkUhat, _N_s, 0.0, _Uinteg, _N_G);
    //Step 2: Calculate gradient wrt ref cords at quadrature points
    blasGemm('N','N', _N_G*D, _N_E*N_F, _N_s, 1, _dphi, _N_G*D, rkUhat, _N_s, 0.0, _dUinteg, _N_G*D);
    //Use the classic GrabCparameters function to get all the items that rk routine needs
    LGrabCparameters(_N_G, _N_E, _Uinteg, _dUinteg, _invJac, rkDivMax, rkLamMax, rkCsMax, _elemCmax);
    /*
    //Steps 3,4,5,6 happen inside the e-g loop.
    for (int a = 0; a < D; a++){
      rkLamMax[a] = 0.0;
      rkDivMax[a] = 0.0;
      rkCsMax[a] = 0.0; }
    
    for (int e = 0; e < _N_E; e++){
      for (int g = 0; g < _N_G; g++){

	//Step 3: Use invJac to get gradient wrt phys cords
	//Step 4: get max wavespeeds
	//Step 5: get max divergence
	//Step 6: Identify max Cs value
      }}
    */
  }

  void SleekAD_RelayParams(scalar* rkLamMax, scalar* rkDivMax, scalar* rkCsMax)
  {
    /*!\brief send calculated global maxima for AD approach 
     */
    //To use this routine, the Max params must have already been calculated with GrabCparameters
    //during the dg residual calculation. This program just relays the known values.
    for (int a = 0; a < D; a++)
      {
	//rk_StoredVariable = DG_StoredVariable (recall that rk CALLS the dgsolver routine)
	rkLamMax[a] = _LamMax[a];
	rkDivMax[a] = _DivMax[a];
	rkCsMax[a]  = _CsMax[a];
      }
  }

#ifdef TGVSTATS
  void TGV_statistics(scalar* U, double time)
  {
    /*!
      \brief Function to gather time-dependent turbulence statistics for HiOCFD5 Taylor-Green vortex problem
      \param[in] U the DG solution
      \param[in] time current solution time, necessary for post-processing temporal derivatives
     */
    /*
      This function integrates both the kinetic energy and the enstrophy
      over the spatial domain. Enstrophy requires velocity derivatives,
      which gets a bit nasty. Modeling usage of InvJac after the 
      evaluate_sf routine in physics/physics.cu
     */
    //PHASE 1: The kinetic energy
    scalar EK = 0.0;
    //zero the elemental EK contributions (stored in dgsolver so I don't have to reallocate every time step)
    for (int j = 0; j < _N_E; j++){_EKlocal[j] = 0.0;}
    //Collocate the solution to element quadrature points.
    blasGemm('N','N', _N_G   , _N_E*N_F, _N_s, 1, _phi,  _N_G   , U, _N_s, 0.0, _Uinteg, _N_G);
    //For each element: calucate local EK contribution, relay to global storage
    scalar piLocal = 3.14159265358979;
    scalar L_tgv = 1.0;
    scalar volume = pow(2*piLocal*L_tgv, 3); //volume of TGV domain (2 pi L cube);
    scalar rho0 = 1.0; //reference density
    for (int e = 0; e < _N_E; e++)
      {
	scalar rho;
	scalar mox;
	scalar moy; 
	scalar moz;
	scalar u_sq;
	scalar v_sq;
	scalar w_sq;
	scalar inv_Rhosq;
	scalar QLocal; //local integrated quantity, stored at quadrature points
	for (int g = 0; g < _N_G; g++)
	  {
	    //gather the conserved variables:
	    rho = _Uinteg[(e*N_F + 0)*_N_G + g];
	    mox = _Uinteg[(e*N_F + 1)*_N_G + g];
	    moy = _Uinteg[(e*N_F + 2)*_N_G + g];
	    moz = _Uinteg[(e*N_F + 3)*_N_G + g];
	    
	    //calculate the squares of velocities:
	    inv_Rhosq = 1.0 / (rho*rho);
	    u_sq = mox*mox * inv_Rhosq;
	    v_sq = moy*moy * inv_Rhosq;
	    w_sq = moz*moz * inv_Rhosq;
	    
	    //Calculate kinetic energy
	    QLocal = 0.5 * rho * (u_sq + v_sq + w_sq);

	    //Multiply Q_local by weights*detJ to get integral contribution
	    //detJxW_full means each quadrature point has distint detJ value,
	    //and that value is multiplied by appropriate quadrature weight
	    _EKlocal[e] += QLocal * _detJxW_full[e*_N_G + g];
	  }
      }
    //Combine the local integrals to get full EK value
    for (int e = 0; e < _N_E; e++)
      {
	EK += _EKlocal[e];
      }
    EK = EK / (rho0*volume);

    //PHASE 2: The Enstrophy integral
    //Going to reuse EKlocal storage
    scalar TGVeps = 0.0;
    //zero the EKlocal array
    for (int j = 0; j < _N_E; j++){_EKlocal[j] = 0.0;}
    //We already have solution at quadrature points, also need gradient wrt ref coords.
    blasGemm('N','N', _N_G*D, _N_E*N_F, _N_s, 1, _dphi, _N_G*D, U, _N_s, 0.0, _dUinteg, _N_G*D);
    for (int e = 0; e < _N_E; e++)
      {
	scalar rho;
	scalar mox;
	scalar moy; 
	scalar moz;
	scalar GradRho_Ref[D];
	scalar GradMox_Ref[D];
	scalar GradMoy_Ref[D];
	scalar GradMoz_Ref[D];
	scalar GradRho_Phys[D];
	scalar GradMox_Phys[D];
	scalar GradMoy_Phys[D];
	scalar GradMoz_Phys[D];
	scalar du_dx;
	scalar du_dy;
	scalar du_dz;
	scalar dv_dx;
	scalar dv_dy;
	scalar dv_dz;
	scalar dw_dx;
	scalar dw_dy;
	scalar dw_dz;
	scalar vorticity[D];
	scalar inv_Rhosq;
	scalar QLocal; //local integrated quantity, stored at quadrature points
	for (int g = 0; g < _N_G; g++)
	  {
	    //gather the conserved variables
	    rho = _Uinteg[(e*N_F + 0)*_N_G + g];
	    mox = _Uinteg[(e*N_F + 1)*_N_G + g];
	    moy = _Uinteg[(e*N_F + 2)*_N_G + g];
	    moz = _Uinteg[(e*N_F + 3)*_N_G + g];

	    //gather the conserved gradients wrt reference coordinates
	    for (int a = 0; a < D; a++)
	      {
		GradRho_Ref[a] = _dUinteg[((e*N_F+0)*_N_G+g)*D+a];
		GradMox_Ref[a] = _dUinteg[((e*N_F+1)*_N_G+g)*D+a];
		GradMoy_Ref[a] = _dUinteg[((e*N_F+2)*_N_G+g)*D+a];
		GradMoz_Ref[a] = _dUinteg[((e*N_F+3)*_N_G+g)*D+a];
	      }

	    //calculate conserved variable gradients wrt physical coordinates
	    for (int a = 0; a < D; a++) //physical coordinate
	      {
		GradRho_Phys[a] = 0.0;
		GradMox_Phys[a] = 0.0;
		GradMoy_Phys[a] = 0.0;
		GradMoz_Phys[a] = 0.0;
		//After some thought, realized the approach below
		//mathches with how I calculate gradient for Cparams;
		//this is good.
		for (int alpha = 0; alpha < D; alpha++) //reference coordinate
		  {
		    GradRho_Phys[a] += GradRho_Ref[alpha] * _invJac[((e*_N_G+g)*D+a)*D+alpha];
		    GradMox_Phys[a] += GradMox_Ref[alpha] * _invJac[((e*_N_G+g)*D+a)*D+alpha];
		    GradMoy_Phys[a] += GradMoy_Ref[alpha] * _invJac[((e*_N_G+g)*D+a)*D+alpha];
		    GradMoz_Phys[a] += GradMoz_Ref[alpha] * _invJac[((e*_N_G+g)*D+a)*D+alpha];
		  }
		
	      }
	    //Calculate most of the velocity derivatives (don't need all components)
	    inv_Rhosq = 1.0 / (rho*rho);
	    //du_dx = inv_Rhosq * (rho*GradMox_Phys[0] - GradRho_Phys[0]*mox);
	    du_dy = inv_Rhosq * (rho*GradMox_Phys[1] - GradRho_Phys[1]*mox);
	    du_dz = inv_Rhosq * (rho*GradMox_Phys[2] - GradRho_Phys[2]*mox);

	    dv_dx = inv_Rhosq * (rho*GradMoy_Phys[0] - GradRho_Phys[0]*moy);
	    //dv_dy = inv_Rhosq * (rho*GradMoy_Phys[1] - GradRho_Phys[1]*moy);
	    dv_dz = inv_Rhosq * (rho*GradMoy_Phys[2] - GradRho_Phys[2]*moy);

	    dw_dx = inv_Rhosq * (rho*GradMoz_Phys[0] - GradRho_Phys[0]*moz);
	    dw_dy = inv_Rhosq * (rho*GradMoz_Phys[1] - GradRho_Phys[1]*moz);
	    //dw_dz = inv_Rhosq * (rho*GradMoz_Phys[2] - GradRho_Phys[2]*moz);

	    //Get the local vorticity
	    vorticity[0] = dw_dy - dv_dz;
	    vorticity[1] = du_dz - dw_dx;
	    vorticity[2] = dv_dx - du_dy;
	    

	    //calculate local enstrophy
	    QLocal = 0.5 * rho * (vorticity[0]*vorticity[0] + vorticity[1]*vorticity[1] + vorticity[2]*vorticity[2]);

	    //Multiply Q_local by weights*detJ to get integral contribution
	    //detJxW_full means each quadrature point has distint detJ value,
	    //and that value is multiplied by appropriate quadrature weight
	    _EKlocal[e] += QLocal * _detJxW_full[e*_N_G + g];
	  }
      } //end element loop to get enstrophy integral contribution
    
    //Combine the local integrals to get full enstrophy value
    for (int e = 0; e < _N_E; e++)
      {
       TGVeps += _EKlocal[e];
      }
    TGVeps = TGVeps / (rho0*volume);
    //printf("time=%f, EK=%f, TGVeps=%f\n", time, EK, TGVeps); 
    //Send the timestep's enstrphy and kinetic energy to appropriate file (using csv because I known how to post-process effectively)
    fprintf(TGVstatsf,"%20.16f,", time); 
    fprintf(TGVstatsf,"%20.16f,", EK); 
    fprintf(TGVstatsf,"%20.16f,", TGVeps); 
    fprintf(TGVstatsf,"\n");
  } //end of TGVstatisits subroutine
#endif //end if for TGVSTATS

};

#endif // DG_SOLVER_H
