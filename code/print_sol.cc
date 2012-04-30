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
// Output multifluid solutions
//===========================================
void print_dg_multifluid(const int N_s, const int N_E, const int N_F, const int model, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all){

 
  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> G  ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))Rho.resize(N_s, N_E);
  if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
  if((all==-1)||(all==2)) Et.resize(N_s, N_E);
  if((all==-1)||(all==3))  G.resize(N_s, N_E);
  if((all==-1)||(all==4))  P.resize(N_s, N_E);

  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U(i,e*N_F+0);
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0)) Rho(i,e) = rho;
      if((all==-1)||(all==1)) Ux (i,e) = U(i,e*N_F+1)/rho;
      if((all==-1)||(all==2)) Et (i,e) = U(i,e*N_F+2);
      if((all==-1)||(all==3)){
	if      (model==0) G  (i,e) = U(i,e*N_F+3)/rho;
	else if (model==1) G  (i,e) = 1+1.0/U(i,e*N_F+3);}
      if((all==-1)||(all==4)) P  (i,e) = (G(i,e)-1)*(Et(i,e) - 0.5*U(i,e*N_F+1)*U(i,e*N_F+1)/rho);
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(Rho, msh_lin,"rho.pos", "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_lin, "ux.pos",  "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(Et,  msh_lin, "et.pos",  "Et", step, time, append);
  if((all==-1)||(all==3)) m.writeSolution( G,  msh_lin,  "g.pos",   "G", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution( P,  msh_lin,  "p.pos",   "P", step, time, append);
}

void print_dg_multifluid(const int N_s, const int N_E, const int N_F, const int model, scalar* U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all){

  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> G  ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))Rho.resize(N_s, N_E);
  if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
  if((all==-1)||(all==2)) Et.resize(N_s, N_E);
  if((all==-1)||(all==3))  G.resize(N_s, N_E);
  if((all==-1)||(all==4))  P.resize(N_s, N_E);
    
  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U[(e*N_F+0)*N_s+i];
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0)) Rho(i,e) = rho;
      if((all==-1)||(all==1)) Ux (i,e) = U[(e*N_F+1)*N_s+i]/rho;
      if((all==-1)||(all==2)) Et (i,e) = U[(e*N_F+2)*N_s+i];
      if((all==-1)||(all==3)){
	if      (model==0) G  (i,e) = U[(e*N_F+3)*N_s+i]/rho;
	else if (model==1) G  (i,e) = 1+1.0/U[(e*N_F+3)*N_s+i];}
      if((all==-1)||(all==4)) P  (i,e) = (G(i,e)-1)*(Et(i,e) - 0.5*U[(e*N_F+1)*N_s+i]*U[(e*N_F+1)*N_s+i]/rho);
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(Rho, msh_lin,"rho.pos", "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_lin, "ux.pos",  "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(Et,  msh_lin, "et.pos",  "Et", step, time, append);
  if((all==-1)||(all==3)) m.writeSolution( G,  msh_lin,  "g.pos",   "G", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution( P,  msh_lin,  "p.pos",   "P", step, time, append);

  // std::string rhoF = "rho.txt"; 
  // FILE *f = fopen(rhoF.c_str(),"w");
  // for (int e=0; e<N_E; e++){
  //   for (int fc=0; fc<N_F; fc++){
  //     for (int i=0; i<N_s; i++){
  // 	fprintf(f,"%i %12.7f %i %i ", step, time, 0, e, ); for(int fc = 0; fc < N_F; fc++) fprintf(f,"%20.16E\t", h_Err1[fc]/N_E);          fprintf(f,"\n");
  //     }
  //   }
  // }
  
}

void print_dg_multifluid_err(const int N_s, const int N_E, const int N_F, const int model, scalar* U, const simpleMesh m, const int msh_lin, const int all){

  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> G  ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))Rho.resize(N_s, N_E);
  if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
  if((all==-1)||(all==2)) Et.resize(N_s, N_E);
  if((all==-1)||(all==3))  G.resize(N_s, N_E);
  if((all==-1)||(all==4))  P.resize(N_s, N_E);
    
  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U[(e*N_F+0)*N_s+i];
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0)) Rho(i,e) = rho;
      if((all==-1)||(all==1)) Ux (i,e) = U[(e*N_F+1)*N_s+i];
      if((all==-1)||(all==2)) Et (i,e) = U[(e*N_F+2)*N_s+i];
      if((all==-1)||(all==3)) G  (i,e) = U[(e*N_F+3)*N_s+i];
      if((all==-1)||(all==4)) P  (i,e) = (G(i,e)-1)*(Et(i,e) - 0.5*U[(e*N_F+1)*N_s+i]*U[(e*N_F+1)*N_s+i]);
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(  Rho,  msh_lin,  "rho_err.pos",  "ErrRho", 0, 0, 0);
  if((all==-1)||(all==1)) m.writeSolution(   Ux,  msh_lin,   "ux_err.pos",   "ErrUx", 0, 0, 0);
  if((all==-1)||(all==2)) m.writeSolution(   Et,  msh_lin,   "et_err.pos",   "ErrEt", 0, 0, 0);
  if((all==-1)||(all==3)) m.writeSolution(    G,  msh_lin,    "g_err.pos",    "ErrG", 0, 0, 0);
  if((all==-1)||(all==5)) m.writeSolution(    P,  msh_lin,    "p_err.pos",    "ErrP", 0, 0, 0);
}


//===========================================
// Output passive solutions
//===========================================
void print_dg_passive(const int N_s, const int N_E, const int N_F, scalar gamma, const fullMatrix<scalar> &U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all){

 
  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> PhiC;
  fullMatrix<scalar> PhiNC;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))  Rho.resize(N_s, N_E);
  if((all==-1)||(all==1))   Ux.resize(N_s, N_E);
  if((all==-1)||(all==2))   Et.resize(N_s, N_E);
  if((all==-1)||(all==3)) PhiC.resize(N_s, N_E);
  if((all==-1)||(all==4))PhiNC.resize(N_s, N_E);
  if((all==-1)||(all==5))    P.resize(N_s, N_E);

  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U(i,e*N_F+0);
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0)) Rho(i,e) = rho;
      if((all==-1)||(all==1)) Ux (i,e) = U(i,e*N_F+1)/rho;
      if((all==-1)||(all==2)) Et (i,e) = U(i,e*N_F+2);
      if((all==-1)||(all==3)) PhiC(i,e)= U(i,e*N_F+3)/rho;
      if((all==-1)||(all==4))PhiNC(i,e)= U(i,e*N_F+4);
      if((all==-1)||(all==5)) P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*U(i,e*N_F+1)*U(i,e*N_F+1)/rho);
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(Rho  , msh_lin,   "rho.pos",  "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(Ux   ,  msh_lin,   "ux.pos",   "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(Et   ,  msh_lin,   "et.pos",   "Et", step, time, append);
  if((all==-1)||(all==3)) m.writeSolution(PhiC ,  msh_lin, "phic.pos", "PhiC", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution(PhiNC,  msh_lin,"phinc.pos","PhiNC", step, time, append);
  if((all==-1)||(all==5)) m.writeSolution( P   ,  msh_lin,    "p.pos",    "P", step, time, append);
}

void print_dg_passive(const int N_s, const int N_E, const int N_F, scalar gamma, scalar* U, const simpleMesh m, const int msh_lin, const int step, const double time, const int append, const int all){

  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> PhiC  ;
  fullMatrix<scalar> PhiNC  ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))  Rho.resize(N_s, N_E);
  if((all==-1)||(all==1))   Ux.resize(N_s, N_E);
  if((all==-1)||(all==2))   Et.resize(N_s, N_E);
  if((all==-1)||(all==3)) PhiC.resize(N_s, N_E);
  if((all==-1)||(all==4))PhiNC.resize(N_s, N_E);
  if((all==-1)||(all==5))    P.resize(N_s, N_E);
    
  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U[(e*N_F+0)*N_s+i];
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0))    Rho(i,e) = rho;
      if((all==-1)||(all==1))    Ux (i,e) = U[(e*N_F+1)*N_s+i]/rho;
      if((all==-1)||(all==2))    Et (i,e) = U[(e*N_F+2)*N_s+i];
      if((all==-1)||(all==3))  PhiC (i,e) = U[(e*N_F+3)*N_s+i]/rho;
      if((all==-1)||(all==4)) PhiNC (i,e) = U[(e*N_F+4)*N_s+i];
      if((all==-1)||(all==5))    P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*U[(e*N_F+1)*N_s+i]*U[(e*N_F+1)*N_s+i]/rho);
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(  Rho,  msh_lin,  "rho.pos",  "Rho", step, time, append);
  if((all==-1)||(all==1)) m.writeSolution(   Ux,  msh_lin,   "ux.pos",   "Ux", step, time, append);
  if((all==-1)||(all==2)) m.writeSolution(   Et,  msh_lin,   "et.pos",   "Et", step, time, append);
  if((all==-1)||(all==3)) m.writeSolution( PhiC,  msh_lin, "phic.pos", "PhiC", step, time, append);
  if((all==-1)||(all==4)) m.writeSolution(PhiNC,  msh_lin,"phinc.pos","PhiNC", step, time, append);
  if((all==-1)||(all==5)) m.writeSolution(    P,  msh_lin,    "p.pos",    "P", step, time, append);
}

void print_dg_passive_err(const int N_s, const int N_E, const int N_F, scalar gamma, scalar* U, const simpleMesh m, const int msh_lin, const int all){

  fullMatrix<scalar> Rho;
  fullMatrix<scalar> Ux ;
  fullMatrix<scalar> Et ;
  fullMatrix<scalar> PhiC  ;
  fullMatrix<scalar> PhiNC  ;
  fullMatrix<scalar> P  ;

  if((all==-1)||(all==0))  Rho.resize(N_s, N_E);
  if((all==-1)||(all==1))   Ux.resize(N_s, N_E);
  if((all==-1)||(all==2))   Et.resize(N_s, N_E);
  if((all==-1)||(all==3)) PhiC.resize(N_s, N_E);
  if((all==-1)||(all==4))PhiNC.resize(N_s, N_E);
  if((all==-1)||(all==5))    P.resize(N_s, N_E);
    
  // separate the fields
  scalar rho = 0;
  for (int e = 0; e < N_E; e++){
    for (int i = 0; i < N_s; i++){

      // Check for NaN error
      rho = U[(e*N_F+0)*N_s+i];
      if(rho != rho){
	printf("NaN error. Code crashed... bummer.\n");
	exit(1);
      }

      if((all==-1)||(all==0))    Rho(i,e) = rho;
      if((all==-1)||(all==1))    Ux (i,e) = U[(e*N_F+1)*N_s+i];
      if((all==-1)||(all==2))    Et (i,e) = U[(e*N_F+2)*N_s+i];
      if((all==-1)||(all==3))  PhiC (i,e) = U[(e*N_F+3)*N_s+i];
      if((all==-1)||(all==4)) PhiNC (i,e) = U[(e*N_F+4)*N_s+i];
      if((all==-1)||(all==5))    P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*Rho(i,e)*Ux(i,e)*Ux(i,e));
    }
  }
  // print to the output file
  if((all==-1)||(all==0)) m.writeSolution(  Rho,  msh_lin,  "rho_err.pos",  "ErrRho", 0, 0, 0);
  if((all==-1)||(all==1)) m.writeSolution(   Ux,  msh_lin,   "ux_err.pos",   "ErrUx", 0, 0, 0);
  if((all==-1)||(all==2)) m.writeSolution(   Et,  msh_lin,   "et_err.pos",   "ErrEt", 0, 0, 0);
  if((all==-1)||(all==3)) m.writeSolution( PhiC,  msh_lin, "phic_err.pos", "ErrPhiC", 0, 0, 0);
  if((all==-1)||(all==4)) m.writeSolution(PhiNC,  msh_lin,"phinc_err.pos","ErrPhiNC", 0, 0, 0);
  if((all==-1)||(all==5)) m.writeSolution(    P,  msh_lin,    "p_err.pos",    "ErrP", 0, 0, 0);
}
