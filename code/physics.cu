#include <physics.h>
#include <basic_fluxes.h>
#include <oned_passive_fluxes.h>
#include <twod_passive_fluxes.h>
#include <oned_multifluid_fluxes.h>
#include <twod_multifluid_fluxes.h>

//==========================================================================
//
// Kernel definitions
//
//==========================================================================

//==========================================================================
arch_global void evaluate_sf_1D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
#elif USE_GPU
      int e = blockIdx.x;
      int g = threadIdx.x;
#endif

      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar Et    = Ug[(e*N_F+2)*N_G+g];
#ifdef MULTIFLUID
#ifdef GAMCONS
      gamma=1+rho/Ug[(e*N_F+3)*N_G+g];
#elif  GAMNCON
      gamma=1+1.0/Ug[(e*N_F+3)*N_G+g];
#endif
#endif
      scalar p = (gamma-1)*(Et - 0.5*rho*u*u);
      scalar EtplusP = Et + p;

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);       
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);
      f[((e*N_F+2)*N_G+g)*D+0] = flux_ab(EtplusP,u);

#ifdef PASSIVE
      s[(e*N_F+3)*N_G+g] = 0;
      s[(e*N_F+4)*N_G+g] = -u*dUg[(e*N_F+4)*N_G+g]*invJac[e*N_G+g];// = -u*dphincdx;
      f[((e*N_F+3)*N_G+g)*D+0] = flux_abc(rho,u,Ug[(e*N_F+3)*N_G+g]/rho); // rho*u*phic
      f[((e*N_F+4)*N_G+g)*D+0] = 0;
#elif MULTIFLUID
#ifdef GAMCONS
      s[(e*N_F+3)*N_G+g] = 0;
      f[((e*N_F+3)*N_G+g)*D+0] = flux_abc(rho,u,1/(gamma-1));
#elif  GAMNCON
      s[(e*N_F+3)*N_G+g] = -u*dUg[(e*N_F+3)*N_G+g]*invJac[e*N_G+g]; // = -u*dalphadx
      f[((e*N_F+3)*N_G+g)*D+0] = 0;
#endif
#endif
      
#ifdef USE_CPU
    }
  }
#endif
}

//==========================================================================
arch_global void evaluate_sf_2D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
#elif USE_GPU
      int e = blockIdx.x;
      int g = threadIdx.x;
#endif

      scalar rho   = Ug[(e*N_F+0)*N_G+g];   
      scalar u     = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v     = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Et    = Ug[(e*N_F+3)*N_G+g];
      scalar vdotv = u*u+v*v;
#ifdef MULTIFLUID
#ifdef GAMCONS
      gamma=1+rho/Ug[(e*N_F+4)*N_G+g];
#elif  GAMNCON
      gamma=1+1.0/Ug[(e*N_F+4)*N_G+g];
#endif
#endif
      scalar p = (gamma-1)*(Et - 0.5*rho*vdotv);
      scalar EtplusP = Et + p;

      // Source term
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
      s[(e*N_F+3)*N_G+g] = 0;

      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = flux_ab(rho,u);      // rho*u     
      f[((e*N_F+1)*N_G+g)*D+0] = flux_ab2pc(rho,u,p);    // rho*u*u + p
      f[((e*N_F+2)*N_G+g)*D+0] = flux_abc(rho,u,v);    // rho*u*v
      f[((e*N_F+3)*N_G+g)*D+0] = flux_ab(EtplusP,u);  // u(E+p)

      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = flux_ab(rho,v);      // rho*v     
      f[((e*N_F+1)*N_G+g)*D+1] = flux_abc(rho,u,v);    // rho*u*v
      f[((e*N_F+2)*N_G+g)*D+1] = flux_ab2pc(rho,v,p);    // rho*v*v + p
      f[((e*N_F+3)*N_G+g)*D+1] = flux_ab(EtplusP,v);  // v(E+p)

#ifdef PASSIVE
      scalar phic = Ug[(e*N_F+4)*N_G+g]/rho; 
      s[(e*N_F+4)*N_G+g] = 0;
      f[((e*N_F+4)*N_G+g)*D+0] = flux_abc(rho,u,phic); // flux wrt x: rho*u*phic
      f[((e*N_F+4)*N_G+g)*D+1] = flux_abc(rho,v,phic); // flux wrt y: rho*v*phic
      
      scalar vdotgradphinc = 0; // = -u*dphincdx-v*dphincdy
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradphinc += 
	  -u*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+5)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
      }
      s[(e*N_F+5)*N_G+g] = vdotgradphinc;
      f[((e*N_F+5)*N_G+g)*D+0] = 0;// flux wrt x
      f[((e*N_F+5)*N_G+g)*D+1] = 0;// flux wrt y
#elif MULTIFLUID
#ifdef GAMCONS
      s[(e*N_F+4)*N_G+g] = 0;
      f[((e*N_F+4)*N_G+g)*D+0] = flux_abc(rho,u,1/(gamma-1));// flux wrt x
      f[((e*N_F+4)*N_G+g)*D+1] = flux_abc(rho,v,1/(gamma-1));// flux wrt y
#elif  GAMNCON
      scalar vdotgradalpha = 0; // = -u*dalphadx-v*dalphady
      for(int alpha = 0; alpha < D; alpha++){
	vdotgradalpha += 
	  -u*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] // dphidx = dphidxi*dxidx + dphideta*detadx
	  -v*dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];// dphidy = dphidxi*dxidy + dphideta*detady
      }
      s[(e*N_F+4)*N_G+g] = vdotgradalpha; 
      f[((e*N_F+4)*N_G+g)*D+0] = 0;// flux wrt x
      f[((e*N_F+4)*N_G+g)*D+1] = 0;// flux wrt y
#endif
#endif
      
#ifdef USE_CPU
    }
  }
#endif
}

//==========================================================================
arch_global void evaluate_q_1D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals){
  
#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    scalar* vap = new scalar[4*(4+4)];
#elif USE_GPU
    int t = blockIdx.x;
    extern __shared__ scalar vap[];
#endif

    scalar nx = normals[t];

    scalar rhoL= UgF[(t*N_F+0)*2+0];
    scalar rhoR= UgF[(t*N_F+0)*2+1];
    scalar uL  = UgF[(t*N_F+1)*2+0]/rhoL;
    scalar uR  = UgF[(t*N_F+1)*2+1]/rhoR;
    scalar EtL = UgF[(t*N_F+2)*2+0];
    scalar EtR = UgF[(t*N_F+2)*2+1];
    scalar phicL  = 0; scalar phicR  = 0;
    scalar phincL = 0; scalar phincR = 0;
    scalar alphaL = 0; scalar alphaR = 0;

#ifdef PASSIVE
    phicL  = UgF[(t*N_F+3)*2+0]/rhoL;
    phicR  = UgF[(t*N_F+3)*2+1]/rhoR;
    phincL = UgF[(t*N_F+4)*2+0];
    phincR = UgF[(t*N_F+4)*2+1];
    alphaL = 1/(gamma-1);
    alphaR = 1/(gamma-1);
#elif MULTIFLUID
#ifdef GAMCONS
    alphaL = UgF[(t*N_F+3)*2+0]/rhoL;
    alphaR = UgF[(t*N_F+3)*2+1]/rhoR;
#elif  GAMNCON
    alphaL = UgF[(t*N_F+3)*2+0];
    alphaR = UgF[(t*N_F+3)*2+1];
#endif
#endif

    scalar gammaL = 1.0+1.0/alphaL;
    scalar gammaR = 1.0+1.0/alphaR;
    scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*uL*uL);
    scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*uR*uR);
    scalar EtPL = EtL+pL;
    scalar EtPR = EtR+pR;

    // Evaluate the right and left eigenvalues
    int sizevap = 4;
    scalar aL = sqrt((gammaL*pL)/rhoL);
    scalar aR = sqrt((gammaR*pR)/rhoR);

    vap[0*sizevap+0] = fabs(uL) - aL;
    vap[0*sizevap+1] = fabs(uL);
    vap[0*sizevap+2] = fabs(uL) + aL;
    vap[0*sizevap+3] = fabs(uL);

    vap[1*sizevap+0] = fabs(uR) - aR;
    vap[1*sizevap+1] = fabs(uR);
    vap[1*sizevap+2] = fabs(uR) + aR;
    vap[1*sizevap+3] = fabs(uR);

    scalar maxvap = 0;
    for (int k = 0; k < 2*sizevap; k++){
      if (maxvap<vap[k]) maxvap = vap[k];
    }
   
    //
    // Evaluate the fluxes on the right and left
    //

    // Rusanov flux
#ifdef RUS
    //first: fx = rho*u; 
    scalar qL = -0.5*((flux_ab(rhoL,uL) + flux_ab(rhoR,uR))*nx
		      -maxvap*(rhoR-rhoL));
    q[(t*N_F+0)*2+0] = qL;
    q[(t*N_F+0)*2+1] = -qL;
      
    //second: fx = rho*u*u+Bx*Bx+Pbar; 
    qL = -0.5*((flux_ab2pc(rhoL,uL,pL)  + flux_ab2pc(rhoR,uR,pR))*nx
	       -maxvap*(rhoR*uR-rhoL*uL));
    q[(t*N_F+1)*2+0] = qL;
    q[(t*N_F+1)*2+1] = -qL;

    //third: fx = EtplusP*u; 
    qL = -0.5*((flux_ab(EtPL,uL) + flux_ab(EtPR,uR))*nx
	       -maxvap*(EtR-EtL));
    q[(t*N_F+2)*2+0] = qL; 
    q[(t*N_F+2)*2+1] = -qL;

    scalar ncterm = 0;
#ifdef PASSIVE
    //fourth: fx = rho*u*phi
    qL = -0.5*((flux_abc(rhoL,uL,phicL) + flux_abc(rhoR,uR,phicR))*nx
	       -maxvap*(rhoR*phicR-rhoL*phicL));
    q[(t*N_F+3)*2+0] = qL; 
    q[(t*N_F+3)*2+1] = -qL;
	
    //fifth:
    qL = -0.5*(-maxvap*(phincR-phincL));
    ncterm = 0.5*0.5*(uL+uR)*(phincL-phincR)*nx;
    q[(t*N_F+4)*2+0] = qL + ncterm; 
    q[(t*N_F+4)*2+1] = -qL+ ncterm;
#elif MULTIFLUID
#ifdef GAMCONS
    qL = -0.5*((flux_abc(rhoL,uL,alphaL) + flux_abc(rhoR,uR,alphaR))*nx
	       -maxvap*(rhoR*alphaR-rhoL*alphaL));
#elif  GAMNCON
    qL = -0.5*(-maxvap*(alphaR-alphaL));
    ncterm = -0.5*0.5*(uL+uR)*(alphaR-alphaL)*nx;
#endif
    q[(t*N_F+3)*2+0] = qL + ncterm; 
    q[(t*N_F+3)*2+1] = -qL+ ncterm;
#endif

    // Non-conservative flux
#elif HLL

    // Wave estimates
    scalar SL = 0; scalar SR = 0;
    int estimate = 0;
    switch (estimate){
    case 1:{ // Toro 10.49
      scalar rhoRoe = sqrt(rhoL*rhoR);
      scalar uRoe = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar HL = (EtL + pL)/rhoL;
      scalar HR = (EtR + pR)/rhoR;
      scalar HRoe = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar alphaRoe = (sqrt(rhoL)*alphaL+sqrt(rhoR)*alphaR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar gammaRoe = 1+1.0/alphaRoe;
      scalar aRoe = sqrt((gammaRoe-1)*(HRoe-0.5*uRoe*uRoe));
      SL = uRoe*nx-aRoe;
      SR = uRoe*nx+aRoe;
      break;}
    case 2: { // Toro 10.52
      scalar eta2 = 0.5*sqrt(rhoL*rhoR)/((sqrt(rhoL)+sqrt(rhoR))*(sqrt(rhoL)+sqrt(rhoR)));
      scalar dbar2 = (sqrt(rhoL)*aL*aL + sqrt(rhoR)*aR*aR)/(sqrt(rhoL)+sqrt(rhoR)) + eta2*(uR-uL)*(uR-uL);
      scalar ubar = 0.5*(uR+uL);
      SL = ubar*nx - sqrt(dbar2);
      SR = ubar*nx + sqrt(dbar2);
      break;}
    default: // Toro 10.48
      SL = MIN((uL*nx-aL),(uR*nx-aR));
      SR = MAX((uL*nx+aL),(uR*nx+aR));
    }
    
    scalar pnc1=0, pnc2=0, pnc3=0, pnc4=0, pnc5=0;
    scalar vnc = 0, vncabs = 0;
#ifdef PASSIVE
    vnc    = -0.5*    (uL*nx+uR*nx)*(phincL-phincR);
    vncabs = -0.5*fabs(uL*nx+uR*nx)*(phincL-phincR);
#elif MULTIFLUID
    vnc    = -0.5*    (uL*nx+uR*nx)*(alphaL-alphaR);
    vncabs = -0.5*fabs(uL*nx+uR*nx)*(alphaL-alphaR);
#endif

    // define the flux
    if (SL > 0){
      pnc1 = flux_ab(rhoL,uL)*nx;
      pnc2 = flux_ab2pc(rhoL,uL,pL)*nx;
      pnc3 = flux_ab(EtPL,uL)*nx;
#ifdef PASSIVE
      pnc4 = flux_abc(rhoL,uL,phicL)*nx;
      pnc5 = -0.5*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
      pnc4 = flux_abc(rhoL,uL,alphaL)*nx;
#elif  GAMNCON
      pnc4 = -0.5*vncabs;
#endif
#endif
    }
    else if ((SL < 0)&&(SR > 0)){
      pnc1 = fhll(rhoL, SL, flux_ab(rhoL,uL)*nx,   rhoR, SR, flux_ab(rhoR,uR)*nx);
      pnc2 = fhll(  uL, SL, flux_ab2pc(rhoL,uL,pL)*nx,  uR, SR, flux_ab2pc(rhoR,uR,pR)*nx);
      pnc3 = fhll( EtL, SL, flux_ab(EtPL,uL)*nx,    EtR, SR, flux_ab(EtPR,uR)*nx);
#ifdef PASSIVE
      pnc4 = fhll(phicL, SL, flux_abc(rhoL,uL,phicL)*nx, phicR, SR, flux_abc(rhoR,uR,phicR)*nx);
      pnc5 = fhll(phincL, SL, 0, phincR, SR, 0) - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
      pnc4 = fhll(alphaL, SL, flux_abc(rhoL,uL,alphaL)*nx, alphaR, SR, flux_abc(rhoR,uR,alphaR)*nx);
#elif  GAMNCON
      pnc4 = fhll(alphaL, SL, 0, alphaR, SR, 0) - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
#endif
#endif
    }
    else if (SR < 0){
      pnc1 = flux_ab(rhoR,uR)*nx;
      pnc2 = flux_ab2pc(rhoR,uR,pR)*nx;
      pnc3 = flux_ab(EtPR,uR)*nx;
#ifdef PASSIVE
      pnc4 = flux_abc(rhoR,uR,phicR)*nx;
      pnc5 = 0.5*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
      pnc4 = flux_abc(rhoR,uR,alphaR)*nx;
#elif  GAMNCON
      pnc4 = 0.5*vncabs;
#endif
#endif
    }

    // first
    scalar qL = -pnc1;
    q[(t*N_F+0)*2+0] = qL;
    q[(t*N_F+0)*2+1] = -qL;
    
    //second: fx = rho*u*u+Bx*Bx+Pbar; 
    qL = -pnc2;
    q[(t*N_F+1)*2+0] = qL;
    q[(t*N_F+1)*2+1] = -qL;
    
    //third: fx = EtplusP*u; 
    qL = -pnc3;
    q[(t*N_F+2)*2+0] = qL; 
    q[(t*N_F+2)*2+1] = -qL;
      
    scalar ncterm = 0;
#ifdef PASSIVE
    //fourth: fx = rho*u*phic
    qL = -pnc4;
    q[(t*N_F+3)*2+0] = qL;
    q[(t*N_F+3)*2+1] = -qL;
      
    //fifth:
    qL = -pnc5;
    ncterm = - 0.5*vnc;
    q[(t*N_F+4)*2+0] = qL  + ncterm; 
    q[(t*N_F+4)*2+1] = -qL + ncterm;
#elif MULTIFLUID
#ifdef GAMNCON
    ncterm = -0.5*vnc;
#endif
    qL = -pnc4;
    q[(t*N_F+3)*2+0] = qL  + ncterm; 
    q[(t*N_F+3)*2+1] = -qL + ncterm;
#endif

    // Non-conservative Roe flux
#elif ROE
    scalar rhoRoe = sqrt(rhoL*rhoR);
    scalar uRoe = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
    scalar HL = (EtL + pL)/rhoL;
    scalar HR = (EtR + pR)/rhoR;
    scalar HRoe = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
    scalar alphaRoe = (sqrt(rhoL)*alphaL+sqrt(rhoR)*alphaR)/(sqrt(rhoL)+sqrt(rhoR));
    scalar gammaRoe = 1+1.0/alphaRoe;
    scalar aRoe = sqrt((gammaRoe-1)*(HRoe-0.5*uRoe*uRoe));
    scalar iL = pL*alphaL;
    scalar iR = pR*alphaR;
    scalar iRoe = (sqrt(rhoL)*iL+sqrt(rhoR)*iR)/(sqrt(rhoL)+sqrt(rhoR));
    scalar DpRoe= (gammaRoe-1)*(gammaRoe-1)*(alphaRoe*(iR-iL) - iRoe*(alphaR-alphaL));
    scalar pRoe = (gammaRoe-1)*iRoe;

    // Roe eigenvalues
    vap[2*sizevap+0] = uRoe-aRoe;
    vap[2*sizevap+1] = uRoe;
    vap[2*sizevap+2] = uRoe+aRoe;
    vap[2*sizevap+3] = uRoe;

    // Roe waves strengths
    vap[3*sizevap+1] = (rhoR-rhoL) - DpRoe/(aRoe*aRoe);
    vap[3*sizevap+2] = (DpRoe + rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
    vap[3*sizevap+0] = (DpRoe - rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
    vap[3*sizevap+3] = alphaR-alphaL;
    // aiRoe[1]= (gamma-1)/(aRoe*aRoe)*((rhoR-rhoL)*(HRoe-uRoe*uRoe)+uRoe*(rhoR*uR-rhoL*uL)-(EtR-EtL));
    // aiRoe[0]=1/(2*aRoe)*((rhoR-rhoL)*(uRoe+aRoe)-(rhoR*uR-rhoL*uL)-aRoe*aiRoe[1]);
    // aiRoe[2]=(rhoR-rhoL)-(aiRoe[0]+aiRoe[1]);

    // Roe eigenvectors
    vap[(4+0)*sizevap+0] = 1;
    vap[(4+0)*sizevap+1] = uRoe-aRoe;
    vap[(4+0)*sizevap+2] = HRoe-uRoe*aRoe;
    vap[(4+0)*sizevap+3] = 0;

    vap[(4+1)*sizevap+0] = 1;
    vap[(4+1)*sizevap+1] = uRoe;
    vap[(4+1)*sizevap+2] = 0.5*uRoe*uRoe;
    vap[(4+1)*sizevap+3] = 0;

    vap[(4+2)*sizevap+0] = 1;
    vap[(4+2)*sizevap+1] = uRoe+aRoe;
    vap[(4+2)*sizevap+2] = HRoe+uRoe*aRoe;
    vap[(4+2)*sizevap+3] = 0;
      
    vap[(4+3)*sizevap+0] = 0;
    vap[(4+3)*sizevap+1] = 0;
    vap[(4+3)*sizevap+2] = pRoe;
    vap[(4+3)*sizevap+3] = 1;

    //first: fx = rho*u;
    // if      (uRoe>0)  qL = flux_ab(rhoL,uL) + aiRoe[0]*vapRoe[0]*vep[0*3+0];
    // else if (uRoe<=0) qL = flux_ab(rhoR,uR) - aiRoe[2]*vapRoe[2]*vep[2*3+0];
    scalar qL = 0;
    qL = 0.5*(flux_ab(rhoL,uL) + flux_ab(rhoR,uR))*nx;
    for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
    q[(t*N_F+0)*2+0] = -qL;
    q[(t*N_F+0)*2+1] = qL;
      
    //second: fx = rho*u*u+Bx*Bx+Pbar; 
    // if      (uRoe>0)  qL = flux_ab2pc(rhoL,uL,pL) + aiRoe[0]*vapRoe[0]*vep[0*3+1];
    // else if (uRoe<=0) qL = flux_ab2pc(rhoR,uR,pR) - aiRoe[2]*vapRoe[2]*vep[2*3+1];
    qL = 0.5*(flux_ab2pc(rhoL,uL,pL)  + flux_ab2pc(rhoR,uR,pR))*nx;
    for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+1];
    q[(t*N_F+1)*2+0] = -qL;
    q[(t*N_F+1)*2+1] = qL;

    //third: fx = EtplusP*u; 
    // if      (uRoe>0)  qL = flux_ab(EtPL,uL) + aiRoe[0]*vapRoe[0]*vep[0*3+2];
    // else if (uRoe<=0) qL = flux_ab(EtPR,uR) - aiRoe[2]*vapRoe[2]*vep[2*3+2];
    qL = 0.5*(flux_ab(EtPL,uL) + flux_ab(EtPR,uR))*nx;
    for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+2];
    q[(t*N_F+2)*2+0] = -qL; 
    q[(t*N_F+2)*2+1] = qL;

    scalar ncterm = 0;
#ifdef PASSIVE
    //fourth: fx = rho*u*phic
    // if      (uRoe>0)  qL = flux_abc(rhoL,uL,phicL) + aiRoe[0]*vapRoe[0]*vep[0*3+1];
    // else if (uRoe<=0) qL = flux_abc(rhoR,uR,phicR) - aiRoe[2]*vapRoe[2]*vep[2*3+1];
    qL = 0.5*(flux_abc(rhoL,uL,phicL) + flux_abc(rhoR,uR,phicR))*nx;
    for(int k=0;k<3;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
    q[(t*N_F+3)*2+0] = -qL; 
    q[(t*N_F+3)*2+1] = qL;
      
    //fifth:
    qL = 0.0;
    for(int k=0;k<3;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
    q[(t*N_F+4)*2+0] = -qL - 0.5*uRoe*(phincR-phincL)*nx;
    q[(t*N_F+4)*2+1] = qL  - 0.5*uRoe*(phincR-phincL)*nx;
#elif MULTIFLUID
#ifdef GAMCONS
    qL = 0.5*(flux_abc(rhoL,uL,alphaL) + flux_abc(rhoR,uR,alphaR))*nx;
#elif  GAMNCON
    qL = 0;
    ncterm = -0.5*uRoe*(alphaR-alphaL)*nx;
#endif
    for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+3];
    q[(t*N_F+3)*2+0] = -qL + ncterm;
    q[(t*N_F+3)*2+1] = qL  + ncterm;
#endif
#endif // end flux ifs'

#ifdef USE_CPU
    delete[] vap;
  }
#endif
}

//==========================================================================
arch_global void evaluate_q_2D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals){
  
#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    //scalar* vap = new scalar[M_G*2*4*(4+4)];
    scalar* vap = new scalar[M_G*2*4];
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
      int t = blockIdx.x;
      int g = threadIdx.x;
      extern __shared__ scalar vap[];
#endif

      scalar nx = normals[t*2+0];
      scalar ny = normals[t*2+1];
      scalar rhoL= UgF[((t*N_F+0)*2+0)*M_G+g];
      scalar rhoR= UgF[((t*N_F+0)*2+1)*M_G+g];
      scalar uL  = UgF[((t*N_F+1)*2+0)*M_G+g]/rhoL;
      scalar uR  = UgF[((t*N_F+1)*2+1)*M_G+g]/rhoR;
      scalar vL  = UgF[((t*N_F+2)*2+0)*M_G+g]/rhoL;
      scalar vR  = UgF[((t*N_F+2)*2+1)*M_G+g]/rhoR;
      scalar EtL = UgF[((t*N_F+3)*2+0)*M_G+g];
      scalar EtR = UgF[((t*N_F+3)*2+1)*M_G+g]; 
      scalar vdotvL = uL*uL+vL*vL;
      scalar vdotvR = uR*uR+vR*vR;
      scalar vdotnL = uL*nx+vL*ny;
      scalar vdotnR = uR*nx+vR*ny;
      scalar phicL  = 0; scalar phicR  = 0;
      scalar phincL = 0; scalar phincR = 0;
      scalar alphaL = 0; scalar alphaR = 0;

#ifdef PASSIVE
      phicL  = UgF[((t*N_F+4)*2+0)*M_G+g]/rhoL;
      phicR  = UgF[((t*N_F+4)*2+1)*M_G+g]/rhoR;
      phincL = UgF[((t*N_F+5)*2+0)*M_G+g];
      phincR = UgF[((t*N_F+5)*2+1)*M_G+g];
      alphaL = 1/(gamma-1);
      alphaR = 1/(gamma-1);
#elif MULTIFLUID
#ifdef GAMCONS
      alphaL = UgF[((t*N_F+4)*2+0)*M_G+g]/rhoL;
      alphaR = UgF[((t*N_F+4)*2+1)*M_G+g]/rhoR;
#elif  GAMNCON
      alphaL = UgF[((t*N_F+4)*2+0)*M_G+g];
      alphaR = UgF[((t*N_F+4)*2+1)*M_G+g];
#endif
#endif

      scalar gammaL = 1.0+1.0/alphaL;
      scalar gammaR = 1.0+1.0/alphaR;
      scalar pL = (gammaL-1)*(EtL - 0.5*rhoL*vdotvL);
      scalar pR = (gammaR-1)*(EtR - 0.5*rhoR*vdotvR);
      scalar EtPL = EtL+pL;
      scalar EtPR = EtR+pR;

      // Evaluate the right and left eigenvalues
      int sizevap = 4;
      scalar aL = sqrt((gammaL*pL)/rhoL);
      scalar aR = sqrt((gammaR*pR)/rhoR);

      vap[(g*2+0)*sizevap+0] = fabs(vdotnL) - aL;
      vap[(g*2+0)*sizevap+1] = fabs(vdotnL);
      vap[(g*2+0)*sizevap+2] = fabs(vdotnL) + aL;
      vap[(g*2+0)*sizevap+3] = fabs(vdotnL);

      vap[(g*2+1)*sizevap+0] = fabs(vdotnR) - aR;
      vap[(g*2+1)*sizevap+1] = fabs(vdotnR);
      vap[(g*2+1)*sizevap+2] = fabs(vdotnR) + aR;
      vap[(g*2+1)*sizevap+3] = fabs(vdotnR);

      scalar maxvap = 0;
      for (int k = 0; k < 2*sizevap; k++){
	if (maxvap<vap[g*2*sizevap+k]) maxvap = vap[g*2*sizevap+k];
      }
      
      //
      // Evaluate the fluxes on the right and left
      //

      // Rusanov flux
#ifdef RUS
      //first: fx = rho*u; fy = rho*v
      scalar qL = -0.5*((flux_ab(rhoL,uL) + flux_ab(rhoR,uR))*nx+
			(flux_ab(rhoL,vL) + flux_ab(rhoR,vR))*ny
			-maxvap*(rhoR-rhoL));
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      
      //second: fx = rho*u*u+p; fy = rho*u*v
      qL = -0.5*((flux_ab2pc(rhoL,uL,pL)  + flux_ab2pc(rhoR,uR,pR))*nx+
		 (flux_abc(rhoL,uL,vL)  + flux_abc(rhoR,uR,vR))*ny
		 -maxvap*(rhoR*uR-rhoL*uL));
      q[((t*N_F+1)*2+0)*M_G+g] = qL;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL;

      //third: fx = rho*u*v; fy = rho*v*v+p
      qL = -0.5*((flux_abc(rhoL,uL,vL)  + flux_abc(rhoR,uR,vR))*nx+
		 (flux_ab2pc(rhoL,vL,pL)  + flux_ab2pc(rhoR,vR,pR))*ny
		 -maxvap*(rhoR*vR-rhoL*vL));
      q[((t*N_F+2)*2+0)*M_G+g] = qL;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL;

      //fourth: fx = EtplusP*u; fy = EtplusP*v;
      qL = -0.5*((flux_ab(EtPL,uL) + flux_ab(EtPR,uR))*nx+
		 (flux_ab(EtPL,vL) + flux_ab(EtPR,vR))*ny
		 -maxvap*(EtR-EtL));
      q[((t*N_F+3)*2+0)*M_G+g] = qL;
      q[((t*N_F+3)*2+1)*M_G+g] = -qL;

      scalar ncterm = 0;
#ifdef PASSIVE
      //fourth: fx = rho*u*phi; fy = rho*v*phi
      qL = -0.5*((flux_abc(rhoL,uL,phicL) + flux_abc(rhoR,uR,phicR))*nx+
		 (flux_abc(rhoL,vL,phicL) + flux_abc(rhoR,vR,phicR))*ny
		 -maxvap*(rhoR*phicR-rhoL*phicL));
      q[((t*N_F+4)*2+0)*M_G+g] = qL;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL;
	
      //fifth:
      qL = -0.5*(-maxvap*(phincR-phincL));
      ncterm = 0.5*0.5*(vdotnL+vdotnR)*(phincL-phincR);
      q[((t*N_F+5)*2+0)*M_G+g] = qL + ncterm;
      q[((t*N_F+5)*2+1)*M_G+g] = -qL+ ncterm;

#elif MULTIFLUID
#ifdef GAMCONS
      qL = -0.5*((flux_abc(rhoL,uL,alphaL) + flux_abc(rhoR,uR,alphaR))*nx+
		 (flux_abc(rhoL,vL,alphaL) + flux_abc(rhoR,vR,alphaR))*ny
		 -maxvap*(rhoR*alphaR-rhoL*alphaL));
#elif  GAMNCON
      qL = -0.5*(-maxvap*(alphaR-alphaL));
      ncterm = -0.5*0.5*(vdotnL+vdotnR)*(alphaR-alphaL);
#endif
      q[((t*N_F+4)*2+0)*M_G+g] = qL + ncterm;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL+ ncterm;
#endif

      // Non-conservative flux
#elif HLL

      // Wave estimates
      scalar SL = 0; scalar SR = 0;
      int estimate = 0;
      switch (estimate){
      case 1:{ // Toro 10.49
	scalar rhoRoe = sqrt(rhoL*rhoR);
	scalar VRoe = (sqrt(rhoL)*vdotnL+sqrt(rhoR)*vdotnR)/(sqrt(rhoL)+sqrt(rhoR));
	scalar HL = (EtL + pL)/rhoL;
	scalar HR = (EtR + pR)/rhoR;
	scalar HRoe = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
	scalar alphaRoe = (sqrt(rhoL)*alphaL+sqrt(rhoR)*alphaR)/(sqrt(rhoL)+sqrt(rhoR));
	scalar gammaRoe = 1+1.0/alphaRoe;
	scalar aRoe = sqrt((gammaRoe-1)*(HRoe-0.5*VRoe*VRoe));
	SL = VRoe-aRoe;
	SR = VRoe+aRoe;
	break;}
      case 2: { // Toro 10.52
	scalar eta2 = 0.5*sqrt(rhoL*rhoR)/((sqrt(rhoL)+sqrt(rhoR))*(sqrt(rhoL)+sqrt(rhoR)));
	scalar dbar2 = (sqrt(rhoL)*aL*aL + sqrt(rhoR)*aR*aR)/(sqrt(rhoL)+sqrt(rhoR)) + eta2*(vdotnR-vdotnL)*(vdotnR-vdotnL);
	scalar Vbar = 0.5*(vdotnR+vdotnL);
	SL = Vbar - sqrt(dbar2);
	SR = Vbar + sqrt(dbar2);
	break;}
      default: // Toro 10.48
	SL = MIN((vdotnL-aL),(vdotnR-aR));
	SR = MAX((vdotnL+aL),(vdotnR+aR));
      }
    
      scalar pnc1=0, pnc2=0, pnc3=0, pnc4=0, pnc5=0, pnc6 = 0;
      scalar vnc = 0, vncabs = 0;
#ifdef PASSIVE
      vnc    = -0.5*    (vdotnL+vdotnR)*(phincL-phincR);
      vncabs = -0.5*fabs(vdotnL+vdotnR)*(phincL-phincR);
#elif MULTIFLUID
      vnc    = -0.5*    (vdotnL+vdotnR)*(alphaL-alphaR);
      vncabs = -0.5*fabs(vdotnL+vdotnR)*(alphaL-alphaR);
#endif

      // define the flux
      if (SL > 0){
	pnc1 = flux_ab(rhoL,uL)*nx       + flux_ab(rhoL,vL)*ny;
	pnc2 = flux_ab2pc(rhoL,uL,pL)*nx    + flux_abc(rhoL,uL,vL)*ny;
	pnc3 = flux_abc(rhoL,uL,vL)*nx    + flux_ab2pc(rhoL,vL,pL)*ny;
	pnc4 = flux_ab(EtPL,uL)*nx       + flux_ab(EtPL,vL)*ny;
#ifdef PASSIVE
	pnc5 = flux_abc(rhoL,uL,phicL)*nx + flux_abc(rhoL,vL,phicL)*ny;
	pnc6 = -0.5*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
	pnc5 = flux_abc(rhoL,uL,alphaL)*nx+flux_abc(rhoL,vL,alphaL)*ny;
#elif  GAMNCON
	pnc5 = -0.5*vncabs;
#endif
#endif
      }
      else if ((SL < 0)&&(SR > 0)){
	pnc1 = fhll(rhoL, SL, flux_ab(rhoL,uL)   *nx + flux_ab(rhoL,vL)   *ny,rhoR, SR, flux_ab(rhoR,uR)   *nx + flux_ab(rhoR,vR)*ny);
	pnc2 = fhll(  uL, SL, flux_ab2pc(rhoL,uL,pL)*nx + flux_abc(rhoL,uL,vL)*ny,  uR, SR, flux_ab2pc(rhoR,uR,pR)*nx + flux_abc(rhoR,uR,vR)*ny);
	pnc3 = fhll(  vL, SL, flux_abc(rhoL,uL,vL)*nx + flux_ab2pc(rhoL,vL,pL)*ny,  vR, SR, flux_ab2pc(rhoR,uR,pR)*nx + flux_ab2pc(rhoR,vR,pR)*ny);
	pnc4 = fhll( EtL, SL, flux_ab(EtPL,uL)   *nx + flux_ab(EtPL,vL)   *ny, EtR, SR, flux_ab(EtPR,uR)   *nx + flux_ab(EtPR,vR)*ny);
#ifdef PASSIVE
	pnc5 = fhll(phicL, SL, flux_abc(rhoL,uL,phicL)*nx + flux_abc(rhoL,vL,phicL)*ny, phicR, SR, flux_abc(rhoR,uR,phicR)*nx + flux_abc(rhoR,vR,phicR)*ny);
	pnc6 = fhll(phincL, SL, 0, phincR, SR, 0) - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
	pnc5 = fhll(alphaL, SL, flux_abc(rhoL,uL,alphaL)*nx+flux_abc(rhoL,vL,alphaL)*ny, alphaR, SR, flux_abc(rhoR,uR,alphaR)*nx+flux_abc(rhoR,vR,alphaR)*ny);
#elif  GAMNCON
	pnc5 = fhll(alphaL, SL, 0, alphaR, SR, 0) - 0.5*fabs(SR+SL)/fabs(SR-SL)*vncabs;
#endif
#endif
      }
      else if (SR < 0){
	pnc1 = flux_ab(rhoR,uR)*nx        + flux_ab(rhoR,vR)*ny;
	pnc2 = flux_ab2pc(rhoR,uR,pR)*nx     + flux_abc(rhoR,uR,vR)*ny;
	pnc3 = flux_abc(rhoR,uR,vR)*nx     + flux_ab2pc(rhoR,vR,pR)*ny;
	pnc4 = flux_ab(EtPR,uR)*nx        + flux_ab(EtPR,vR)*ny;
#ifdef PASSIVE
	pnc5 = flux_abc(rhoR,uR,phicR)*nx  + flux_abc(rhoR,vR,phicR)*ny;
	pnc6 = 0.5*vncabs;
#elif MULTIFLUID
#ifdef GAMCONS
	pnc5 = flux_abc(rhoR,uR,alphaR)*nx+ flux_abc(rhoR,vR,alphaR)*ny;
#elif  GAMNCON
	pnc5 = 0.5*vncabs;
#endif
#endif
      }

      // first
      scalar qL = -pnc1;
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
    
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      qL = -pnc2;
      q[((t*N_F+1)*2+0)*M_G+g] = qL ;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL;

      //third: fx = rho*v*v+Bx*Bx+Pbar; 
      qL = -pnc3;
      q[((t*N_F+2)*2+0)*M_G+g] = qL ;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL;
      
      //fourth: fx = EtplusP*u; 
      qL = -pnc4;
      q[((t*N_F+3)*2+0)*M_G+g] = qL ;
      q[((t*N_F+3)*2+1)*M_G+g] = -qL;

      scalar ncterm = 0;
#ifdef PASSIVE
      //fifth: fx = rho*u*phic
      qL = -pnc5;
      q[((t*N_F+4)*2+0)*M_G+g] = qL ;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL;
      
      //sixth:
      qL = -pnc6;
      ncterm = - 0.5*vnc;
      q[((t*N_F+5)*2+0)*M_G+g] = qL + ncterm;
      q[((t*N_F+5)*2+1)*M_G+g] = -qL+ ncterm;
#elif MULTIFLUID
#ifdef GAMNCON
      ncterm = -0.5*vnc;
#endif
      qL = -pnc5;
      q[((t*N_F+4)*2+0)*M_G+g] = qL + ncterm;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL+ ncterm;
#endif

      // Non-conservative Roe flux
#elif ROE
      scalar rhoRoe = sqrt(rhoL*rhoR);
      scalar uRoe = (sqrt(rhoL)*vdotnL+sqrt(rhoR)*vdotnR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar HL = (EtL + pL)/rhoL;
      scalar HR = (EtR + pR)/rhoR;
      scalar HRoe = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar alphaRoe = (sqrt(rhoL)*alphaL+sqrt(rhoR)*alphaR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar gammaRoe = 1+1.0/alphaRoe;
      scalar aRoe = sqrt((gammaRoe-1)*(HRoe-0.5*uRoe*uRoe));
      scalar iL = pL*alphaL;
      scalar iR = pR*alphaR;
      scalar iRoe = (sqrt(rhoL)*iL+sqrt(rhoR)*iR)/(sqrt(rhoL)+sqrt(rhoR));
      scalar DpRoe= (gammaRoe-1)*(gammaRoe-1)*(alphaRoe*(iR-iL) - iRoe*(alphaR-alphaL));
      scalar pRoe = (gammaRoe-1)*iRoe;

      // Roe eigenvalues
      vap[2*sizevap+0] = uRoe-aRoe;
      vap[2*sizevap+1] = uRoe;
      vap[2*sizevap+2] = uRoe+aRoe;
      vap[2*sizevap+3] = uRoe;

      // Roe waves strengths
      vap[3*sizevap+1] = (rhoR-rhoL) - DpRoe/(aRoe*aRoe);
      vap[3*sizevap+2] = (DpRoe + rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
      vap[3*sizevap+0] = (DpRoe - rhoRoe*aRoe*(uR-uL))/(2*aRoe*aRoe);
      vap[3*sizevap+3] = alphaR-alphaL;
      // aiRoe[1]= (gamma-1)/(aRoe*aRoe)*((rhoR-rhoL)*(HRoe-uRoe*uRoe)+uRoe*(rhoR*uR-rhoL*uL)-(EtR-EtL));
      // aiRoe[0]=1/(2*aRoe)*((rhoR-rhoL)*(uRoe+aRoe)-(rhoR*uR-rhoL*uL)-aRoe*aiRoe[1]);
      // aiRoe[2]=(rhoR-rhoL)-(aiRoe[0]+aiRoe[1]);

      // Roe eigenvectors
      vap[(4+0)*sizevap+0] = 1;
      vap[(4+0)*sizevap+1] = uRoe-aRoe;
      vap[(4+0)*sizevap+2] = HRoe-uRoe*aRoe;
      vap[(4+0)*sizevap+3] = 0;

      vap[(4+1)*sizevap+0] = 1;
      vap[(4+1)*sizevap+1] = uRoe;
      vap[(4+1)*sizevap+2] = 0.5*uRoe*uRoe;
      vap[(4+1)*sizevap+3] = 0;

      vap[(4+2)*sizevap+0] = 1;
      vap[(4+2)*sizevap+1] = uRoe+aRoe;
      vap[(4+2)*sizevap+2] = HRoe+uRoe*aRoe;
      vap[(4+2)*sizevap+3] = 0;
      
      vap[(4+3)*sizevap+0] = 0;
      vap[(4+3)*sizevap+1] = 0;
      vap[(4+3)*sizevap+2] = pRoe;
      vap[(4+3)*sizevap+3] = 1;

      //first: fx = rho*u;
      // if      (uRoe>0)  qL = flux_ab(rhoL,uL) + aiRoe[0]*vapRoe[0]*vep[0*3+0];
      // else if (uRoe<=0) qL = flux_ab(rhoR,uR) - aiRoe[2]*vapRoe[2]*vep[2*3+0];
      scalar qL = 0;
      qL = 0.5*(flux_ab(rhoL,uL) + flux_ab(rhoR,uR))*nx;
      for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
      q[(t*N_F+0)*2+0] = -qL;
      q[(t*N_F+0)*2+1] = qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; 
      // if      (uRoe>0)  qL = flux_ab2pc(rhoL,uL,pL) + aiRoe[0]*vapRoe[0]*vep[0*3+1];
      // else if (uRoe<=0) qL = flux_ab2pc(rhoR,uR,pR) - aiRoe[2]*vapRoe[2]*vep[2*3+1];
      qL = 0.5*(flux_ab2pc(rhoL,uL,pL)  + flux_ab2pc(rhoR,uR,pR))*nx;
      for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+1];
      q[(t*N_F+1)*2+0] = -qL;
      q[(t*N_F+1)*2+1] = qL;

      //third: fx = EtplusP*u; 
      // if      (uRoe>0)  qL = flux_ab(EtPL,uL) + aiRoe[0]*vapRoe[0]*vep[0*3+2];
      // else if (uRoe<=0) qL = flux_ab(EtPR,uR) - aiRoe[2]*vapRoe[2]*vep[2*3+2];
      qL = 0.5*(flux_ab(EtPL,uL) + flux_ab(EtPR,uR))*nx;
      for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+2];
      q[(t*N_F+2)*2+0] = -qL; 
      q[(t*N_F+2)*2+1] = qL;

      scalar ncterm = 0;
#ifdef PASSIVE
      //fourth: fx = rho*u*phic
      // if      (uRoe>0)  qL = flux_abc(rhoL,uL,phicL) + aiRoe[0]*vapRoe[0]*vep[0*3+1];
      // else if (uRoe<=0) qL = flux_abc(rhoR,uR,phicR) - aiRoe[2]*vapRoe[2]*vep[2*3+1];
      qL = 0.5*(flux_abc(rhoL,uL,phicL) + flux_abc(rhoR,uR,phicR))*nx;
      for(int k=0;k<3;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
      q[(t*N_F+3)*2+0] = -qL; 
      q[(t*N_F+3)*2+1] = qL;
      
      //fifth:
      qL = 0.0;
      for(int k=0;k<3;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+0];
      q[(t*N_F+4)*2+0] = -qL - 0.5*uRoe*(phincR-phincL)*nx;
      q[(t*N_F+4)*2+1] = qL  - 0.5*uRoe*(phincR-phincL)*nx;
#elif MULTIFLUID
#ifdef GAMCONS
      qL = 0.5*(flux_abc(rhoL,uL,alphaL) + flux_abc(rhoR,uR,alphaR))*nx;
#elif  GAMNCON
      qL = 0;
      ncterm = -0.5*uRoe*(alphaR-alphaL)*nx;
#endif
      for(int k=0;k<4;k++) qL += -0.5*vap[3*sizevap+k]*fabs(vap[2*sizevap+k])*vap[(4+k)*sizevap+3];
      q[(t*N_F+3)*2+0] = -qL + ncterm;
      q[(t*N_F+3)*2+1] = qL  + ncterm;
#endif
#endif // end flux ifs'

#ifdef USE_CPU
    } // end loop on g
    delete[] vap;
  } // end loop on t
#endif
}

//==========================================================================
arch_global void evaluate_q(int M_G, int M_T, int N_F, int D, scalar* q, scalar* UgF, scalar* normals){


  scalar* uL = new scalar[N_F]; // conserved variables
  scalar* uR = new scalar[N_F]; 
  scalar* n  = new scalar[D];   // normals
  scalar* F  = new scalar[N_F];   // fluxes
  scalar* ncterm = new scalar[N_F];
  
#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    //scalar* vap = new scalar[M_G*2*4*(4+4)];
    scalar* vap = new scalar[M_G*2*4];
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
      int t = blockIdx.x;
      int g = threadIdx.x;
      extern __shared__ scalar vap[];
#endif

      // Initialize these variables
      for(int fc = 0; fc < N_F; fc++){
	uL[fc] = UgF[((t*N_F+fc)*2+0)*M_G+g];
	uR[fc] = UgF[((t*N_F+fc)*2+1)*M_G+g];
	F[fc]  = 0;
	ncterm[fc] = 0;
      }
      for(int alpha = 0; alpha < D; alpha++) n[alpha] = normals[t*D+alpha];
      
      // Send the data to the Riemann solvers
#ifdef ONED

#ifdef PASSIVE

#ifdef RUS
      oned_passive_rusanov(uL,uR,n,F,ncterm);
#elif HLL
      oned_passive_hll(uL,uR,n,F,ncterm);
#elif ROE
      oned_passive_roe(uL,uR,n,F,ncterm);
#endif // flux if

#elif MULTIFLUID

#ifdef RUS
      oned_multifluid_rusanov(uL,uR,n,F,ncterm);
#elif HLL
      oned_multifluid_hll(uL,uR,n,F,ncterm);
#elif ROE
      oned_multifluid_roe(uL,uR,n,F,ncterm);
#endif // flux if

#endif // physics if

#elif TWOD

#ifdef PASSIVE

#ifdef RUS
      twod_passive_rusanov(uL,uR,n,F,ncterm);
#elif HLL
      twod_passive_hll(uL,uR,n,F,ncterm);
#elif ROE
      twod_passive_roe(uL,uR,n,F,ncterm);
#endif // flux if

#elif MULTIFLUID

#ifdef RUS
      twod_multifluid_rusanov(uL,uR,n,F,ncterm);
#elif HLL
      //twod_multifluid_hll(uL,uR,n,F,ncterm);
#elif ROE
      //twod_multifluid_roe(uL,uR,n,F,ncterm);
#endif // flux if

#endif // physics if

#endif // dimension if

      // Apply the fluxes
      for(int fc = 0; fc < N_F; fc++){
	q[((t*N_F+fc)*2+0)*M_G+g] =-F[fc] + ncterm[fc];
	q[((t*N_F+fc)*2+1)*M_G+g] = F[fc] + ncterm[fc];
      }

      
#ifdef USE_CPU
    } // end loop on g
    delete[] vap;
  } // end loop on t
#endif

  // Free the vectors
  delete[] uL;
  delete[] uR;
  delete[] n;
  delete[] F;
  
}


      
//==========================================================================
//
//  Host C functions
//
//==========================================================================
extern "C" 
void Levaluate_sf_1D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

#ifdef USE_GPU
  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
#endif

  evaluate_sf_1D arch_args (D, N_G, N_E, N_F, gamma, s, f, Ug, dUg, invJac);
}

extern "C" 
void Levaluate_q_1D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals){

#ifdef USE_GPU
  dim3 dimBlock(1,1,1);
  dim3 dimGrid(M_T,1);
#endif

  evaluate_q_1D arch_args_array(4*(4+4)*sizeof(scalar)) (M_G, M_T, N_F, gamma, q, UgF, normals);
}

extern "C" 
void Levaluate_sf_2D(int D, int N_G, int N_E, int N_F, scalar gamma, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac){

#ifdef USE_GPU
  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
#endif

  evaluate_sf_2D arch_args (D, N_G, N_E, N_F, gamma, s, f, Ug, dUg, invJac);
}

extern "C" 
void Levaluate_q_2D(int M_G, int M_T, int N_F, scalar gamma, scalar* q, scalar* UgF, scalar* normals){

#ifdef USE_GPU
  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
#endif

  evaluate_q_2D arch_args_array(M_G*4*(4+4)*sizeof(scalar)) (M_G, M_T, N_F, gamma, q, UgF, normals);
}

extern "C" 
void Levaluate_q(int M_G, int M_T, int N_F, int D, scalar* q, scalar* UgF, scalar* normals){

#ifdef USE_GPU
  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
#endif

  evaluate_q arch_args_array(M_G*4*(4+4)*sizeof(scalar)) (M_G, M_T, N_F, D, q, UgF, normals);
}



//
// Possibly broken
// 
//==========================================================================
arch_global void cpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
#elif USE_GPU
      int e = blockIdx.x;
      int g = threadIdx.x;
#endif
  
      scalar eta  =  Ug[(e*N_F+0)*N_G+g];
  
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;
      s[(e*N_F+2)*N_G+g] = 0;
  
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = H0*Ug[(e*N_F+1)*N_G+g]; // u_x
      f[((e*N_F+1)*N_G+g)*D+0] = G0*eta; // eta 
      f[((e*N_F+2)*N_G+g)*D+0] = 0; 
  
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = H0*Ug[(e*N_F+2)*N_G+g]; // u_y
      f[((e*N_F+1)*N_G+g)*D+1] = 0;
      f[((e*N_F+2)*N_G+g)*D+1] = G0*eta; // eta

#ifdef USE_CPU
    }
  }
#endif
}

arch_device scalar cpu_flux1_mhd(scalar rho, scalar u){return rho*u;}                                        // for f0X, f0X, f0X
arch_device scalar cpu_flux2_mhd(scalar rho, scalar u, scalar Bx, scalar pbar){return rho*u*u-Bx*Bx+pbar;}   // for f1X, f2Y, f3Z
arch_device scalar cpu_flux3_mhd(scalar rho, scalar u, scalar v, scalar Bx, scalar By){return rho*u*v-Bx*By;}// for f1Y, f1Z, f2X, f2Z, f3X, f3Y
arch_device scalar cpu_flux4_mhd(scalar u, scalar By, scalar v, scalar Bx){return u*By-v*Bx;}                // for f4Y, f4Z, f5X, f5Z, f6X, f6Y
arch_device scalar cpu_flux5_mhd(scalar EtplusPbar, scalar u, scalar vdotB, scalar Bx) {return EtplusPbar*u - vdotB*Bx;} // for f7X, f7Y, f7Z

//==========================================================================
arch_global void cpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

#ifdef USE_CPU
  for(int e = 0; e < N_E; e++){
    for(int g = 0; g < N_G; g++){
#elif USE_GPU
      int e = blockIdx.x;
      int g = threadIdx.x;
#endif

      scalar rho = Ug[(e*N_F+0)*N_G+g];   
      scalar u   = Ug[(e*N_F+1)*N_G+g]/rho;  // (rho u / rho) = u
      scalar v   = Ug[(e*N_F+2)*N_G+g]/rho;  // (rho v / rho) = v
      scalar Bx  = Ug[(e*N_F+3)*N_G+g];
      scalar By  = Ug[(e*N_F+4)*N_G+g];
      scalar Et  = Ug[(e*N_F+5)*N_G+g];
      scalar vdotv  = u*u+v*v;       
      scalar BdotB  = Bx*Bx+By*By;
      scalar vdotB  = u*Bx+v*By;
      scalar divB   = 0;
      for(int alpha = 0; alpha < D; alpha++){
	divB += dUg[((e*N_F+3)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+0)*D+alpha] + dUg[((e*N_F+4)*N_G+g)*D+alpha]*invJac[((e*N_G+g)*D+1)*D+alpha];
      }
      scalar p = (gamma-1)*(Et - 0.5*(rho*vdotv+BdotB));
      scalar Pbar = p+0.5*BdotB;
      scalar EtPbar = Et+Pbar;

      // Powel source
      s[(e*N_F+0)*N_G+g] = 0;
      s[(e*N_F+1)*N_G+g] = 0;//-divB*Bx;
      s[(e*N_F+2)*N_G+g] = 0;//-divB*By;
      s[(e*N_F+3)*N_G+g] = 0;//-divB*u;
      s[(e*N_F+4)*N_G+g] = 0;//-divB*v;
      s[(e*N_F+5)*N_G+g] = 0;//-divB*vdotB;
            
      // Flux derive par rapport a x
      f[((e*N_F+0)*N_G+g)*D+0] = cpu_flux1_mhd(rho,u);                //rho*u; 
      f[((e*N_F+1)*N_G+g)*D+0] = cpu_flux2_mhd(rho,u,Bx,Pbar);        //rho*u*u-Bx*Bx+Pbar; 
      f[((e*N_F+2)*N_G+g)*D+0] = cpu_flux3_mhd(rho,u,v,Bx,By);        //rho*u*v-Bx*By; 
      f[((e*N_F+3)*N_G+g)*D+0] = 0;                                   //0;
      f[((e*N_F+4)*N_G+g)*D+0] = cpu_flux4_mhd(u,By,v,Bx);            //u*By-v*Bx;
      f[((e*N_F+5)*N_G+g)*D+0] = cpu_flux5_mhd(EtPbar,u,vdotB,Bx);    //EtplusPbar*u-vdotB*Bx;
      
      // Flux derive par rapport a y
      f[((e*N_F+0)*N_G+g)*D+1] = cpu_flux1_mhd(rho,v);                //rho*v;
      f[((e*N_F+1)*N_G+g)*D+1] = cpu_flux3_mhd(rho,v,u,By,Bx);        //rho*v*u-By*Bx;
      f[((e*N_F+2)*N_G+g)*D+1] = cpu_flux2_mhd(rho,v,By,Pbar);        //rho*v*v-By*By+Pbar;
      f[((e*N_F+3)*N_G+g)*D+1] = cpu_flux4_mhd(v,Bx,u,By);            //v*Bx-u*By;
      f[((e*N_F+4)*N_G+g)*D+1] = 0;                                   //0;
      f[((e*N_F+5)*N_G+g)*D+1] = cpu_flux5_mhd(EtPbar,v,vdotB,By);    //EtplusPbar*v-vdotB*By;

#ifdef USE_CPU
    }
  }
#endif
}

//==========================================================================
arch_global void cpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
      int t = blockIdx.x;
      int g = threadIdx.x;
#endif

      scalar nx = normals[t*2+0];
      scalar ny = normals[t*2+1];
      scalar etaL= UgF[((t*N_F+0)*2+0)*M_G+g];
      scalar etaR= UgF[((t*N_F+0)*2+1)*M_G+g];
      scalar uLn = UgF[((t*N_F+1)*2+0)*M_G+g] * nx + UgF[((t*N_F+2)*2+0)*M_G+g] * ny;
      scalar uRn = UgF[((t*N_F+1)*2+1)*M_G+g] * nx + UgF[((t*N_F+2)*2+1)*M_G+g] * ny;
      
      scalar h0 = H0;
      scalar g0 = G0;
      
      // first equation
      scalar qL = -0.5*h0*(uLn + uRn + sqrt(g0/h0)*(etaL-etaR)); // Left
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      // second
      qL = -0.5*g0*nx*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
      q[((t*N_F+1)*2+0)*M_G+g] = qL;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL;
      // third
      qL = -0.5*g0*ny*(etaL+etaR+sqrt(h0/g0)*(uLn-uRn)); // Left
      q[((t*N_F+2)*2+0)*M_G+g] = qL;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL;

#ifdef USE_CPU
    }
  }
#endif
}


//==========================================================================
arch_global void cpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){
  
#ifdef USE_CPU
  for(int t = 0; t < M_T; t++){
    scalar* vap = new scalar[M_G*2*8];
    for(int g = 0; g < M_G; g++){
#elif USE_GPU
      int t = blockIdx.x;
      int g = threadIdx.x;
      extern __shared__ scalar vap[];
#endif

      scalar nx = normals[t*2+0];
      scalar ny = normals[t*2+1];
      scalar rhoL= UgF[((t*N_F+0)*2+0)*M_G+g];
      scalar rhoR= UgF[((t*N_F+0)*2+1)*M_G+g];
      scalar uL  = UgF[((t*N_F+1)*2+0)*M_G+g]/rhoL;
      scalar uR  = UgF[((t*N_F+1)*2+1)*M_G+g]/rhoR;
      scalar vL  = UgF[((t*N_F+2)*2+0)*M_G+g]/rhoL;
      scalar vR  = UgF[((t*N_F+2)*2+1)*M_G+g]/rhoR;
      scalar BxL = UgF[((t*N_F+3)*2+0)*M_G+g];
      scalar BxR = UgF[((t*N_F+3)*2+1)*M_G+g];
      scalar ByL = UgF[((t*N_F+4)*2+0)*M_G+g];
      scalar ByR = UgF[((t*N_F+4)*2+1)*M_G+g];
      scalar EtL = UgF[((t*N_F+5)*2+0)*M_G+g];
      scalar EtR = UgF[((t*N_F+5)*2+1)*M_G+g];
      scalar vdotvL = uL*uL+vL*vL;
      scalar vdotvR = uR*uR+vR*vR;
      scalar BdotBL = BxL*BxL+ByL*ByL;
      scalar BdotBR = BxR*BxR+ByR*ByR;
      scalar vdotBL = uL*BxL+vL*ByL;
      scalar vdotBR = uR*BxR+vR*ByR;
      scalar vdotnL = uL*nx+vL*ny;
      scalar vdotnR = uR*nx+vR*ny;
      scalar BdotnL = BxL*nx+ByL*ny;
      scalar BdotnR = BxR*nx+ByR*ny;
      //scalar aveBdotn = 0.5*(BdotnL + BdotnR);      
      scalar pL = (gamma-1)*(EtL - 0.5*(rhoL*vdotvL+BdotBL));
      scalar pR = (gamma-1)*(EtR - 0.5*(rhoR*vdotvR+BdotBR));
      scalar PbarL = pL+0.5*BdotBL;
      scalar PbarR = pR+0.5*BdotBR;
      scalar EtPbarL = EtL+PbarL;
      scalar EtPbarR = EtR+PbarR;

      // Evaluate the right and left eigenvalues
      int sizevap = 8;
      scalar alfenL  = BdotnL/sqrt(rhoL);
      scalar alfenR  = BdotnR/sqrt(rhoR);
      scalar a2L = (gamma*pL+BdotBL)/rhoL;
      scalar a2R = (gamma*pR+BdotBR)/rhoR;
      scalar cfL = sqrt(0.5*(a2L + sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
      scalar cfR = sqrt(0.5*(a2R + sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));
      scalar csL = sqrt(0.5*(a2L - sqrt( a2L*a2L - 4.0*gamma*pL*BdotnL*BdotnL/(rhoL*rhoL))));
      scalar csR = sqrt(0.5*(a2R - sqrt( a2R*a2R - 4.0*gamma*pR*BdotnR*BdotnR/(rhoR*rhoR))));

      vap[(g*2+0)*sizevap+0] = fabs(vdotnL);
      vap[(g*2+0)*sizevap+1] = fabs(vdotnL) + alfenL;
      vap[(g*2+0)*sizevap+2] = fabs(vdotnL) - alfenL;
      vap[(g*2+0)*sizevap+3] = fabs(vdotnL) + cfL;
      vap[(g*2+0)*sizevap+4] = fabs(vdotnL) - cfL;
      vap[(g*2+0)*sizevap+5] = fabs(vdotnL) + csL;
      vap[(g*2+0)*sizevap+6] = fabs(vdotnL) - csL;
      vap[(g*2+0)*sizevap+7] = fabs(vdotnL);
      
      vap[(g*2+1)*sizevap+0] = fabs(vdotnR);
      vap[(g*2+1)*sizevap+1] = fabs(vdotnR) + alfenR;
      vap[(g*2+1)*sizevap+2] = fabs(vdotnR) - alfenR;
      vap[(g*2+1)*sizevap+3] = fabs(vdotnR) + cfR;
      vap[(g*2+1)*sizevap+4] = fabs(vdotnR) - cfR;
      vap[(g*2+1)*sizevap+5] = fabs(vdotnR) + csR;
      vap[(g*2+1)*sizevap+6] = fabs(vdotnR) - csR;
      vap[(g*2+1)*sizevap+7] = fabs(vdotnR);

      scalar maxvap = 0;
      for (int k = 0; k < 2*sizevap; k++){
	if (maxvap<vap[g*16+k]) maxvap = vap[g*16+k];
      }

      // Upwinding on the source term
      // scalar upBdotnL = aveBdotn;
      // scalar upBdotnR = aveBdotn;
      scalar upBdotnL = 0.0;
      scalar upBdotnR = 0.0;
      if      (vdotnL >= 0) upBdotnL = BdotnL;
      else if (vdotnL <  0) upBdotnL = BdotnR;
      if      (vdotnR >= 0) upBdotnR = BdotnL;
      else if (vdotnR <  0) upBdotnR = BdotnR;  
      
      //
      // Evaluate the fluxes on the right and left
      //

      //first: fx = rho*u; fy = rho*v; fz = rho*w; 
      scalar qL = -0.5*((cpu_flux1_mhd(rhoL,uL) + cpu_flux1_mhd(rhoR,uR))*nx +
			(cpu_flux1_mhd(rhoL,vL) + cpu_flux1_mhd(rhoR,vR))*ny 
			-maxvap*(rhoR-rhoL));
      q[((t*N_F+0)*2+0)*M_G+g] = qL;
      q[((t*N_F+0)*2+1)*M_G+g] = -qL;
      
      //second: fx = rho*u*u+Bx*Bx+Pbar; fy = rho*v*u-By*Bx; fz = rho*w*u-Bz*Bx;
      qL = -0.5*((cpu_flux2_mhd(rhoL,uL,BxL,PbarL)  + cpu_flux2_mhd(rhoR,uR,BxR,PbarR) )*nx +
		 (cpu_flux3_mhd(rhoL,vL,uL,ByL,BxL) + cpu_flux3_mhd(rhoR,vR,uR,ByR,BxR))*ny
		 -maxvap*(rhoR*uR-rhoL*uL));
      q[((t*N_F+1)*2+0)*M_G+g] = qL  - BxL*upBdotnL;
      q[((t*N_F+1)*2+1)*M_G+g] = -qL + BxR*upBdotnR;

      //third: fx = rho*u*v-Bx*By; fy = rho*v*v-By*By+Pbar; fz = rho*w*v-Bz*By;
      qL = -0.5*((cpu_flux3_mhd(rhoL,uL,vL,BxL,ByL) + cpu_flux3_mhd(rhoR,uR,vR,BxR,ByR))*nx +
		 (cpu_flux2_mhd(rhoL,vL,ByL,PbarL)  + cpu_flux2_mhd(rhoR,vR,ByR,PbarR) )*ny
		 -maxvap*(rhoR*vR-rhoL*vL));
      q[((t*N_F+2)*2+0)*M_G+g] = qL  - ByL*upBdotnL;
      q[((t*N_F+2)*2+1)*M_G+g] = -qL + ByR*upBdotnR;
      
      //fourth: fx = 0;  fy = v*Bx-u*By; fz = w*Bx-u*Bz;
      qL = -0.5*((0                            + 0                           )*nx +
		 (cpu_flux4_mhd(vL,BxL,uL,ByL) + cpu_flux4_mhd(vR,BxR,uR,ByR))*ny
		 -maxvap*(BxR-BxL));
      q[((t*N_F+3)*2+0)*M_G+g] = qL  - uL*upBdotnL;
      q[((t*N_F+3)*2+1)*M_G+g] = -qL + uR*upBdotnR;

      //fifth: fx = u*By-v*Bx;  fy = 0; fz = w*By-v*Bz;
      qL = -0.5*((cpu_flux4_mhd(uL,ByL,vL,BxL) + cpu_flux4_mhd(uR,ByR,vR,BxR))*nx +
		 (0                            + 0                           )*ny 
		 -maxvap*(ByR-ByL));
      q[((t*N_F+4)*2+0)*M_G+g] = qL  - vL*upBdotnL;
      q[((t*N_F+4)*2+1)*M_G+g] = -qL + vR*upBdotnR;

      //sixth: fx = EtplusPbar*u-vdotB*Bx; fy = EtplusPbar*v-vdotB*By; fz = EtplusPbar*w-vdotB*Bz;
      qL = -0.5*((cpu_flux5_mhd(EtPbarL,uL,vdotBL,BxL) + cpu_flux5_mhd(EtPbarR,uR,vdotBR,BxR))*nx +
		 (cpu_flux5_mhd(EtPbarL,vL,vdotBL,ByL) + cpu_flux5_mhd(EtPbarR,vR,vdotBR,ByR))*ny 
		 -maxvap*(EtR-EtL));
      q[((t*N_F+5)*2+0)*M_G+g] = qL  - vdotBL*upBdotnL;
      q[((t*N_F+5)*2+1)*M_G+g] = -qL + vdotBR*upBdotnR;

#ifdef USE_CPU
    }
    delete[] vap;
  }
#endif
}

extern "C" 
void Lcpu_evaluate_sf_shallow(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar H0, scalar G0){

#ifdef USE_GPU
  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
#endif

  cpu_evaluate_sf_shallow arch_args (D, N_G, N_E, N_F, s, f, Ug, H0, G0);
}

extern "C" 
void Lcpu_evaluate_sf_mhd(int D, int N_G, int N_E, int N_F, scalar* s, scalar* f, scalar* Ug, scalar* dUg, scalar* invJac, scalar gamma){

#ifdef USE_GPU
  dim3 dimBlock(N_G,1,1);
  dim3 dimGrid(N_E,1);
#endif

  cpu_evaluate_sf_mhd arch_args (D, N_G, N_E, N_F, s, f, Ug, dUg, invJac, gamma);
}

extern "C" 
void Lcpu_evaluate_q_shallow(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar H0, scalar G0, scalar* normals){

#ifdef USE_GPU
  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
#endif

  cpu_evaluate_q_shallow arch_args (M_G, M_T, N_F, q, UgF, H0, G0, normals);
}

extern "C" 
void Lcpu_evaluate_q_mhd(int M_G, int M_T, int N_F, scalar* q, scalar* UgF, scalar gamma, scalar* normals){

#ifdef USE_GPU
  dim3 dimBlock(M_G,1,1);
  dim3 dimGrid(M_T,1);
#endif

  cpu_evaluate_q_mhd arch_args_array(M_G*2*8*sizeof(scalar)) (M_G, M_T, N_F, q, UgF, gamma, normals);
}