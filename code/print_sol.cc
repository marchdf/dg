// #include <print_sol.h>

// //===========================================
// // Output Shallow Water solutions
// //===========================================
// void print_dg_shallow(const int N_s, const int N_E,  const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append){

 
//   fullMatrix<scalar> H (N_s, N_E);
//   fullMatrix<scalar> Ux(N_s, N_E);
//   fullMatrix<scalar> Uy(N_s, N_E);
  
//   // separate the fields
//   for (int e = 0; e < N_E; e++){
//     for (int i = 0; i < N_s; i++){
//       H (i,e) = U(i,e*N_F+0);
//       Ux(i,e) = U(i,e*N_F+1);
//       Uy(i,e) = U(i,e*N_F+2);
//     }
//   }
//   // print to the output file
//   m.writeSolution(H , msh_tri,  "h.pos", "H" , step, time, append);
//   m.writeSolution(Ux, msh_tri, "ux.pos", "Ux", step, time, append);
//   m.writeSolution(Uy, msh_tri, "uy.pos", "Uy", step, time, append);
// }

// void print_dg_shallow(const int N_s, const int N_E,  scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append){

 
//   fullMatrix<scalar> H (N_s, N_E);
//   fullMatrix<scalar> Ux(N_s, N_E);
//   fullMatrix<scalar> Uy(N_s, N_E);
  
//   // separate the fields
//   for (int e = 0; e < N_E; e++){
//     for (int i = 0; i < N_s; i++){
//       H (i,e) = U[(e*N_F+0)*N_s+i];
//       Ux(i,e) = U[(e*N_F+1)*N_s+i];
//       Uy(i,e) = U[(e*N_F+2)*N_s+i];
//     }
//   }
//   // print to the output file
//   m.writeSolution(H , msh_tri,  "h.pos", "H" , step, time, append);
//   m.writeSolution(Ux, msh_tri, "ux.pos", "Ux", step, time, append);
//   m.writeSolution(Uy, msh_tri, "uy.pos", "Uy", step, time, append);
// }

// //===========================================
// // Output MHD solutions
// //===========================================
// void print_dg_mhd(const int N_s, const int N_E,  const fullMatrix<scalar> &U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma){

 
//   fullMatrix<scalar> Rho;
//   fullMatrix<scalar> Ux ;
//   fullMatrix<scalar> Uy ;
//   fullMatrix<scalar> Bx ;
//   fullMatrix<scalar> By ;
//   fullMatrix<scalar> Et ;
//   fullMatrix<scalar> P  ;

//   if((all==-1)||(all==0))Rho.resize(N_s, N_E);
//   if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
//   if((all==-1)||(all==2)) Uy.resize(N_s, N_E);
//   if((all==-1)||(all==3)) Bx.resize(N_s, N_E);
//   if((all==-1)||(all==4)) By.resize(N_s, N_E);
//   if((all==-1)||(all==5)) Et.resize(N_s, N_E);
//   if((all==-1)||(all==6))  P.resize(N_s, N_E);
  
//   // separate the fields
//   for (int e = 0; e < N_E; e++){
//     for (int i = 0; i < N_s; i++){
//       if((all==-1)||(all==0)) Rho(i,e) = U(i,e*N_F+0);
//       if((all==-1)||(all==1)) Ux (i,e) = U(i,e*N_F+1)/Rho(i,e);
//       if((all==-1)||(all==2)) Uy (i,e) = U(i,e*N_F+2)/Rho(i,e);
//       if((all==-1)||(all==3)) Bx (i,e) = U(i,e*N_F+3);
//       if((all==-1)||(all==4)) By (i,e) = U(i,e*N_F+4);
//       if((all==-1)||(all==5)) Et (i,e) = U(i,e*N_F+5);
//       if((all==-1)||(all==6)) P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*(Rho(i,e)*(Ux(i,e)*Ux(i,e)+Uy(i,e)*Uy(i,e))
// 								   + Bx(i,e)*Bx(i,e)+By(i,e)*By(i,e)));
//     }
//   }
//   // print to the output file
//   if((all==-1)||(all==0)) m.writeSolution(Rho, msh_tri,"rho.pos", "Rho", step, time, append);
//   if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_tri, "ux.pos",  "Ux", step, time, append);
//   if((all==-1)||(all==2)) m.writeSolution(Uy,  msh_tri, "uy.pos",  "Uy", step, time, append);
//   if((all==-1)||(all==3)) m.writeSolution(Bx,  msh_tri, "bx.pos",  "Bx", step, time, append);
//   if((all==-1)||(all==4)) m.writeSolution(By,  msh_tri, "by.pos",  "By", step, time, append);
//   if((all==-1)||(all==5)) m.writeSolution(Et,  msh_tri, "et.pos",  "Et", step, time, append);
//   if((all==-1)||(all==6)) m.writeSolution( P,  msh_tri,  "p.pos",   "P", step, time, append);
// }

// void print_dg_mhd(const int N_s, const int N_E,  scalar* U, const simpleMesh m, const int msh_tri, const int step, const double time, const int append, const int all, const scalar gamma){

//   fullMatrix<scalar> Rho;
//   fullMatrix<scalar> Ux ;
//   fullMatrix<scalar> Uy ;
//   fullMatrix<scalar> Bx ;
//   fullMatrix<scalar> By ;
//   fullMatrix<scalar> Et ;
//   fullMatrix<scalar> P  ;

//   if((all==-1)||(all==0))Rho.resize(N_s, N_E);
//   if((all==-1)||(all==1)) Ux.resize(N_s, N_E);
//   if((all==-1)||(all==2)) Uy.resize(N_s, N_E);
//   if((all==-1)||(all==3)) Bx.resize(N_s, N_E);
//   if((all==-1)||(all==4)) By.resize(N_s, N_E);
//   if((all==-1)||(all==5)) Et.resize(N_s, N_E);
//   if((all==-1)||(all==6))  P.resize(N_s, N_E);
  
//   // separate the fields
//   for (int e = 0; e < N_E; e++){
//     for (int i = 0; i < N_s; i++){
//       if((all==-1)||(all==0)) Rho(i,e) = U[(e*N_F+0)*N_s+i];
//       if((all==-1)||(all==1)) Ux (i,e) = U[(e*N_F+1)*N_s+i]/Rho(i,e);
//       if((all==-1)||(all==2)) Uy (i,e) = U[(e*N_F+2)*N_s+i]/Rho(i,e);
//       if((all==-1)||(all==3)) Bx (i,e) = U[(e*N_F+3)*N_s+i];
//       if((all==-1)||(all==4)) By (i,e) = U[(e*N_F+4)*N_s+i];
//       if((all==-1)||(all==5)) Et (i,e) = U[(e*N_F+5)*N_s+i];
//       if((all==-1)||(all==6)) P  (i,e) = (gamma-1)*(Et(i,e) - 0.5*(Rho(i,e)*(Ux(i,e)*Ux(i,e)+Uy(i,e)*Uy(i,e))
// 								   + Bx(i,e)*Bx(i,e)+By(i,e)*By(i,e)));
//     }
//   }
//   // print to the output file
//   if((all==-1)||(all==0)) m.writeSolution(Rho, msh_tri,"rho.pos", "Rho", step, time, append);
//   if((all==-1)||(all==1)) m.writeSolution(Ux,  msh_tri, "ux.pos",  "Ux", step, time, append);
//   if((all==-1)||(all==2)) m.writeSolution(Uy,  msh_tri, "uy.pos",  "Uy", step, time, append);
//   if((all==-1)||(all==4)) m.writeSolution(Bx,  msh_tri, "bx.pos",  "Bx", step, time, append);
//   if((all==-1)||(all==5)) m.writeSolution(By,  msh_tri, "by.pos",  "By", step, time, append);
//   if((all==-1)||(all==7)) m.writeSolution(Et,  msh_tri, "et.pos",  "Et", step, time, append);
//   if((all==-1)||(all==8)) m.writeSolution( P,  msh_tri,  "p.pos",   "P", step, time, append);
// }

