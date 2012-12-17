#include <print_sol.h>

//===========================================
// Output Shallow Water solutions
//===========================================
void print_dg_shallow(const int N_s, const int N_E, const int N_F, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append){

 
  fullMatrix<scalar> H (N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> Uy(N_s, N_E);
  
  // separate the fields
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){
      H (i,e) = U(i,e*N_F+0);
      Ux(i,e) = U(i,e*N_F+1);
      Uy(i,e) = U(i,e*N_F+2);
    }
  }
  // print to the output file
  m.writeSolution(H , msh_tri,  "h.pos", "H" , step, time, append);
  m.writeSolution(Ux, msh_tri, "ux.pos", "Ux", step, time, append);
  m.writeSolution(Uy, msh_tri, "uy.pos", "Uy", step, time, append);
}

void print_dg_shallow(const int N_s, const int N_E, const int N_F, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append){

 
  fullMatrix<scalar> H (N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> Uy(N_s, N_E);
  
  // separate the fields
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){
      H (i,e) = U[(e*N_F+0)*N_s+i];
      Ux(i,e) = U[(e*N_F+1)*N_s+i];
      Uy(i,e) = U[(e*N_F+2)*N_s+i];
    }
  }
  // print to the output file
  m.writeSolution(H , msh_tri,  "h.pos", "H" , step, time, append);
  m.writeSolution(Ux, msh_tri, "ux.pos", "Ux", step, time, append);
  m.writeSolution(Uy, msh_tri, "uy.pos", "Uy", step, time, append);
}

//===========================================
// Output MHD solutions
//===========================================
void print_dg_mhd(const int N_s, const int N_E, const int N_F, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma){

 
  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Uy ;
  fullMatrix<scalar> Bx ;
  fullMatrix<scalar> By ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))Rho.resize(N_s, N_E);
  if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
  if((all==-1)||(all==2)) Uy.resize(N_s, N_E);
  if((all==-1)||(all==3)) Bx.resize(N_s, N_E);
  if((all==-1)||(all==4)) By.resize(N_s, N_E);
  if((all==-1)||(all==5)) Et.resize(N_s, N_E);
  if((all==-1)||(all==6))  P.resize(N_s, N_E);
  
  // separate the fields
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){
      if((all==-1)||(all==0)) Rho(i,e) = U(i,e*N_F+0);
      if((all==-1)||(all==1)) Ux (i,e) = U(i,e*N_F+1)/Rho(i,e);
      if((all==-1)||(all==2)) Uy (i,e) = U(i,e*N_F+2)/Rho(i,e);
      if((all==-1)||(all==3)) Bx (i,e) = U(i,e*N_F+3);
      if((all==-1)||(all==4)) By (i,e) = U(i,e*N_F+4);
      if((all==-1)||(all==5)) Et (i,e) = U(i,e*N_F+5);
      if((all==-1)||(all==6)) P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*(Rho(i,e)*(Ux(i,e)*Ux(i,e)+Uy(i,e)*Uy(i,e))
								   + Bx(i,e)*Bx(i,e)+By(i,e)*By(i,e)));
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(Rho, msh_tri,"rho.pos", "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_tri, "ux.pos",  "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(Uy,  msh_tri, "uy.pos",  "Uy", step, time, append);
  if((all==-1)||(all==3)) m.writeSolution(Bx,  msh_tri, "bx.pos",  "Bx", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution(By,  msh_tri, "by.pos",  "By", step, time, append);
  if((all==-1)||(all==5)) m.writeSolution(Et,  msh_tri, "et.pos",  "Et", step, time, append);
  if((all==-1)||(all==6)) m.writeSolution( P,  msh_tri,  "p.pos",   "P", step, time, append);
}

void print_dg_mhd(const int N_s, const int N_E, const int N_F, scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma){

  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Uy ;
  fullMatrix<scalar> Bx ;
  fullMatrix<scalar> By ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))Rho.resize(N_s, N_E);
  if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
  if((all==-1)||(all==2)) Uy.resize(N_s, N_E);
  if((all==-1)||(all==3)) Bx.resize(N_s, N_E);
  if((all==-1)||(all==4)) By.resize(N_s, N_E);
  if((all==-1)||(all==5)) Et.resize(N_s, N_E);
  if((all==-1)||(all==6))  P.resize(N_s, N_E);
  
  // separate the fields
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){
      if((all==-1)||(all==0)) Rho(i,e) = U[(e*N_F+0)*N_s+i];
      if((all==-1)||(all==1)) Ux (i,e) = U[(e*N_F+1)*N_s+i]/Rho(i,e);
      if((all==-1)||(all==2)) Uy (i,e) = U[(e*N_F+2)*N_s+i]/Rho(i,e);
      if((all==-1)||(all==3)) Bx (i,e) = U[(e*N_F+3)*N_s+i];
      if((all==-1)||(all==4)) By (i,e) = U[(e*N_F+4)*N_s+i];
      if((all==-1)||(all==5)) Et (i,e) = U[(e*N_F+5)*N_s+i];
      if((all==-1)||(all==6)) P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*(Rho(i,e)*(Ux(i,e)*Ux(i,e)+Uy(i,e)*Uy(i,e))
								   + Bx(i,e)*Bx(i,e)+By(i,e)*By(i,e)));
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(Rho, msh_tri,"rho.pos", "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_tri, "ux.pos",  "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(Uy,  msh_tri, "uy.pos",  "Uy", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution(Bx,  msh_tri, "bx.pos",  "Bx", step, time, append);
  if((all==-1)||(all==5)) m.writeSolution(By,  msh_tri, "by.pos",  "By", step, time, append);
  if((all==-1)||(all==7)) m.writeSolution(Et,  msh_tri, "et.pos",  "Et", step, time, append);
  if((all==-1)||(all==8)) m.writeSolution( P,  msh_tri,  "p.pos",   "P", step, time, append);
}

//===========================================
// Output generic solutions
//===========================================
void print_dg(const int N_s, const int N_E, const int N_F, scalar gamma, scalar* U, const simpleMesh m, const int elem_type, const int step, const double time, const int append){

  fullMatrix<scalar> Rho(N_s, N_E);
  fullMatrix<scalar> Ux(N_s, N_E);
  fullMatrix<scalar> Et(N_s, N_E);
  fullMatrix<scalar> P(N_s, N_E);
#ifdef PASSIVE
  fullMatrix<scalar> PhiC(N_s, N_E)  ;
  fullMatrix<scalar> PhiNC(N_s, N_E)  ;
#elif  MULTIFLUID  
  fullMatrix<scalar> G(N_s, N_E)  ;
#endif
  
#ifdef ONED
      int Uxind = 1;
      int Eind  = 2;
      int Phicind = 3;
      int Phincind = 4;
      int Gind = 3;
#elif TWOD
      int Uxind = 1;
      int Uyind = 2; fullMatrix<scalar> Uy(N_s, N_E);
      int Eind  = 3;
      int Phicind = 4;
      int Phincind = 5;
      int Gind = 4;
#endif      
  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U[(e*N_F+0)*N_s+i];
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	//exit(1);
      }

      Rho(i,e) = rho;
      Ux (i,e) = U[(e*N_F+Uxind)*N_s+i]/rho;
#ifdef TWOD
      Uy (i,e) = U[(e*N_F+Uyind)*N_s+i]/rho;
#endif	
      Et (i,e) = U[(e*N_F+Eind)*N_s+i];
#ifdef PASSIVE
      PhiC  (i,e) = U[(e*N_F+Phicind)*N_s+i]/rho;
      PhiNC (i,e) = U[(e*N_F+Phincind)*N_s+i];
      P     (i,e) = (gamma-1)*(Et(i,e) - 0.5*U[(e*N_F+Uxind)*N_s+i]*U[(e*N_F+Uxind)*N_s+i]/rho);
#ifdef TWOD
      P     (i,e) = P(i,e) - (gamma-1)*0.5*U[(e*N_F+Uyind)*N_s+i]*U[(e*N_F+Uyind)*N_s+i]/rho;
#endif
#elif MULTIFLUID
#ifdef GAMCONS
      G  (i,e) = 1+rho/U[(e*N_F+Gind)*N_s+i];
#elif  GAMNCON
      G  (i,e) = 1+1.0/U[(e*N_F+Gind)*N_s+i];
#endif
      P  (i,e) = (G(i,e)-1)*(Et(i,e) - 0.5*U[(e*N_F+Uxind)*N_s+i]*U[(e*N_F+Uxind)*N_s+i]/rho);
#ifdef TWOD
      P  (i,e) = P(i,e) - (G(i,e)-1)*0.5*U[(e*N_F+Uyind)*N_s+i]*U[(e*N_F+Uyind)*N_s+i]/rho;
#endif
#endif
    }
  }

  // print to the output file
  m.writeSolution(  Rho,  elem_type,  "rho.pos",  "Rho", step, time, append);
  m.writeSolution(   Ux,  elem_type,   "ux.pos",   "Ux", step, time, append);
#ifdef TWOD
  m.writeSolution(   Uy,  elem_type,   "uy.pos",   "Uy", step, time, append);
#endif
  m.writeSolution(   Et,  elem_type,   "et.pos",   "Et", step, time, append);
#ifdef PASSIVE
  m.writeSolution( PhiC,  elem_type, "phic.pos", "PhiC", step, time, append);
  m.writeSolution(PhiNC,  elem_type,"phinc.pos","PhiNC", step, time, append);
  m.writeSolution(    P,  elem_type,    "p.pos",    "P", step, time, append);
#elif MULTIFLUID
  m.writeSolution( G,  elem_type,  "g.pos",   "G", step, time, append);
  m.writeSolution( P,  elem_type,  "p.pos",   "P", step, time, append);
#endif
}
