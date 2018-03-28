//08/12/2016: gmsh data translator for 3D calculations.
//Currently set up for p1 elements, I don't know if it will work past that

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <iostream>
#include<vector>
#include <fstream>
#include <sstream>


const double pi = 3.141592653589793238;
const double e = 2.71828182845904523536;

const int N = 2; //DG discretization order (like, p1,p2,...)
const int ND = 2;
const int Nfield = 4;
const int plus = 1;
const double PERIOD = 50.0;

const int climit = 15;
const int interval_out = 100;

const int RK1_yes        = 0; //if 1, switch from ideal RK scheme to RK1

//basis choice: 10 is okay lagrange basis
//11 is normalized smart monomial basis (my personal creation, I hope it works)
const int basis_choice_ICB = 1; //1=normalized bi-unit Legendre
const int basis_choice_DG = 2; 
/*
  basis_choice_DG:
  0 is shifted Legendre, needs to be used if limiting is applied
  1 is normalized Legendre
  2 is the Lagrange basis, MUST BE ORGANIZED LIKE GMSH ELEMENT

  10: Actual Lagrange basis, except p2 is for bi-unit domain
  11: Hierarchical basis, bad cond. for icb on side 1
  12: Similar to 11, with a fix for icb conditioning
/*
//Solution space paramters if working on triangles:
const int Nsolut = (N+1)*(N+2)/2;
const int delN = N+1;//N+0; //this describes growth from DG to icb solution
*/
//Solution space parameters for tensor product basis on quads
const int Nsolut = pow(N+1,ND);

const int N_N = 2*ND; //sides per element. Should make easier transition
const char* mesh_file = "mesh.msh";
const char* ux_file =   "ux0000000023.pos";
const char* uy_file =   "uy0000000023.pos";
const char* rho_file = "rho0000000023.pos";
const char* p_file =     "p0000000023.pos";
//Mesh tolerance parameters for periodicity fx, High-Order workshop
//const double msh_tol = 0.3*pow(10,-2); //for grid0
const double msh_tol = 1.0*pow(10,-3); //for grid1. Also grid2
//const double msh_tol = 0.5*pow(10,-3); //for grid3
//const double msh_tol = 0.2*pow(10,-4) //for grid4


//Mesh resolution is coming directly from gmsh file
//See get_mesh_tri and get_c2n_tri for how I set these
int RES_nodes; //total vertex count
int RES_cells; //total element count

/*
int RES_nodes = 289;
int RES_cells = 512; //total count of computational cells
int Mx = 2*sqrt(RES_cells/2); //Mx,My just for data processing usage
int My = sqrt(RES_cells/2);
*/
int t_count; //used for surface flux looping

//Matrix Inversion Parameters
const int master_max = 1;
double cond_mult = 1.0;
double mat_tol = pow(10,-15);



//Gaussian Quadrature Resolution. N+2 for recovery on Cartesian elements
//in this case, GQvol is going to be the full 2d count of quad nodes
const int GQvol = (N+1)*(N+1); //quad resolution over element volume
const int GQsurf = N+1; //this is quad resolution alomg a line (like, edge of an element)
const int GQinit = GQvol; //this is for error quantification and matrix formation.
//output resolution
const int GQout = GQvol;
const double eps = pow(10,-13.0);
double weights_vol[GQvol];
double weights_surf[GQsurf];
double weights_init[GQinit];
class timer 
{
public:
  clock_t start;
  std::string name;
  double elapsed;
  
  timer(std::string n) : name(n), elapsed(0.0f) {}
  ~timer(){}
  
  void startBlock() 
  {
    start = clock();
  }
  
  void endBlock() 
  {
    elapsed += (double)(clock() - start)/CLOCKS_PER_SEC;
  }
  
  void print() {
    printf("%s total time: %fs\n", name.c_str(), elapsed);
  }
};
std::vector<timer*> timers;

timer subTimer_quadSurf(" SURFACE QUADRATURE");

timer subTimer_quadVol ("  VOLUME QUADRATURE");


typedef struct
{
  double loc;
  double weight;
} gwn;
typedef struct //x,y components of position vector
{
  double rx;
  double ry;
  double rz;
} cord;
typedef struct //i,j vector components
{
  double ei;
  double ej;
} vec_dir;
typedef struct
{
  int ei;
  int ej;
} int_dir;
typedef struct
{
  double ma;
  double mx;
  double my;
  double en;
} vec4; //vector of conserved variables
typedef struct
{
  int ma;
  int mx;
  int my;
  int en;
} vec4_int; //integer vector, conserved vars
typedef struct 
{
  double dxx;
  double dx;
  double dy;
  double dyy;
  double dxy;
  double dt;
} grad6; //holds various derivatives of a scalar for source term calculation
typedef struct
{
  vec_dir ma;
  vec_dir mx;
  vec_dir my;
  vec_dir en;
} vec4_dir; //directional vector for each of the conserved variables; used for gradient
typedef struct //2 by 2 matrix
{
  double nw;
  double ne;
  double sw;
  double se;
} mat_22;
typedef struct
{
  cord z0a0; //zeta=0, ada=0
  cord z1a0;
  cord z1a1;
  cord z0a1;
} cell_verts; //stores 4 cell vertices for each element
typedef struct 
{
  double r0;
  double r1;
  double r2;
  double r3;
} mat_41; //4*1 vector, for Transformation vectors
typedef std::vector<double> D1_double; //array of doubles
typedef std::vector<std::vector<double > > D2_double; //double array of doubles
typedef std::vector<std::vector<std::vector<double > > > D3_double; //triple array of doubles
typedef std::vector<std::vector<std::vector<std::vector<double> > > > D4_double; //quadruple array of doubles
typedef std::vector<vec4 > D1_vec4; //single array for vec4
typedef std::vector<std::vector<vec4  > > D2_vec4; //double array for vec4
typedef std::vector<std::vector<std::vector<vec4 > > > D3_vec4; //triple array for vec4

typedef std::vector<cord> D1_cord;
typedef std::vector<std::vector<cord> > D2_cord;

typedef std::vector<std::vector<int > > D2_integer; //double array of integers
typedef std::vector<std::vector<std::vector<int > > > D3_integer; //triple array of integers

typedef std::vector<vec_dir> D1_vecdir;
typedef std::vector<std::vector<vec_dir > > D2_vecdir;
typedef std::vector<std::vector<std::vector<vec_dir> > > D3_vecdir;

typedef std::vector<std::vector<vec4_dir> > D2_vec4dir;
typedef std::vector<std::vector<std::vector<vec4_dir> > > D3_vec4dir;

//Data structure allpcation tools:
D1_double Allocate_1D_double(int d1)
{
  D1_double output;
  output.resize(d1);
  return output;
}
D2_double Allocate_2D_double(int d1, int d2)
{
  std::vector<std::vector<double> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
    }
  for (int A = 0; A < d1; A++)
    {
      for (int B = 0; B < d2; B++)
	{
	  output[A][B] = 0.0;
	}
    }
  return output;
}
D3_double Allocate_3D_double(int d1, int d2, int d3)
{
  D3_double output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	}
    }
  return output;
}
D4_double Allocate_4D_double(int d1, int d2, int d3, int d4)
{
  D4_double output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	  for (int C = 0; C < d3; C++)
	    {
	      output[A][B][C].resize(d4);
	    }
	}
    }
  return output;
}

D1_vec4 Allocate_1D_vec4(int d1)
{
  D1_vec4 output;
  output.resize(d1);
  return output;
}
D2_vec4 Allocate_2D_vec4(int d1, int d2)
{
  D2_vec4 output;
  output.resize(d1);
  for (int A = 0; A < d1; A++)
    {
      output[A].resize(d2);
    }
  return output;
}
D3_vec4 Allocate_3D_vec4(int d1, int d2, int d3)
{
  D3_vec4 output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	}
    }
  return output;
}

D2_vec4dir Allocate_2D_vec4dir(int d1, int d2)
{
  D2_vec4dir output;
  output.resize(d1);
  for (int A = 0; A < d1; A++)
    {
      output[A].resize(d2);
    }
  return output;
}
D3_vec4dir Allocate_3D_vec4dir(int d1, int d2, int d3)
{
  D3_vec4dir output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	}
    }
  return output;
}


D1_cord Allocate_1D_cord(int d1)
{
  D1_cord output;
  output.resize(d1);
  return output;
}
D2_cord Allocate_2D_cord(int d1, int d2)
{
  std::vector<std::vector<cord> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
    }
  return output;
}

D2_integer Allocate_2D_integer(int d1, int d2)
{
  std::vector<std::vector<int> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
    }
  for (int n1 = 0; n1 < d1; n1++)
    {
      for (int n2 = 0; n2 < d2; n2++)
	{
	  output[n1][n2] = 0.0;
	}
    }
  return output;
}
D3_integer Allocate_3D_integer(int d1, int d2, int d3)
{
  D3_integer output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B++)
	{
	  output[A][B].resize(d3);
	}
    }
  
  for (int n1 = 0; n1 < d1; n1++)
    {
      for (int n2 = 0; n2 < d2; n2++)
	{
	  for (int n3 = 0; n3 < d3; n3++)
	    {
	      output[n1][n2][n3] = 0.0;
	    }
	}
    }
  return output;
}

D1_vecdir Allocate_1D_vecdir(int d1)
{
  D1_vecdir output;
  output.resize(d1);
  return output;
}
D2_vecdir Allocate_2D_vecdir(int d1, int d2)
{
  D2_vecdir output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
    }
  return output;
}
D3_vecdir Allocate_3D_vecdir(int d1, int d2, int d3)
{
  D3_vecdir output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	}
    }
  return output;
}

//quadrature parameters
gwn gwn_return(int dim, int spot)
{
  gwn output;
   //Gaussian weights and nodes:
  double weight20[20];	//nodal weights
  double GQ20[20];	//nodal points
  
  double weight10[10];
  double GQ10[10];

  if (dim == 8)
    {
      printf("CATASTROPHE in gwn_return\n");
    }
  
  double weight7[7];
  double GQ7[7];

  double weight6[6];
  double GQ6[6];

  double weight5[5];
  double GQ5[5];
  
  double weight4[4];
  double GQ4[4];
  
  double weight3[3];
  double GQ_loc_3[3];
  
  double weight2[2];
  double GQ_loc_2[2];

  
  weight20[9] = 0.1527533871307258506980843;
  weight20[8] = 0.1491729864726037467878287;
  weight20[7] = 0.1420961093183820513292983;
  weight20[6] = 0.1316886384491766268984945;
  weight20[5] = 0.1181945319615184173123774;
  weight20[4] = 0.1019301198172404350367501;
  weight20[3] = 0.0832767415767047487247581;
  weight20[2] = 0.0626720483341090635695065;
  weight20[1] = 0.0406014298003869413310400;
  weight20[0] = 0.0176140071391521183118620;
  GQ20[9] = -0.0765265211334973337546404;
  GQ20[8] = -0.2277858511416450780804962;
  GQ20[7] = -0.3737060887154195606725482;
  GQ20[6] = -0.5108670019508270980043641;
  GQ20[5] = -0.6360536807265150254528367;
  GQ20[4] = -0.7463319064601507926143051;
  GQ20[3] = -0.8391169718222188233945291;
  GQ20[2] = -0.9122344282513259058677524;
  GQ20[1] = -0.9639719272779137912676661;
  GQ20[0] = -0.9931285991850949247861224;
  
  weight10[4] = 0.2955242247147528701738930;
  weight10[3] = 0.2692667193099963550912269;
  weight10[2] = 0.2190863625159820439955349;
  weight10[1] = 0.1494513491505805931457763;
  weight10[0] = 0.0666713443086881375935688;
  GQ10[4] = -0.1488743389816312108848260;
  GQ10[3] = -0.4333953941292471907992659;
  GQ10[2] = -0.6794095682990244062343274;
  GQ10[1] = -0.8650633666889845107320967;
  GQ10[0] = -0.9739065285171717200779640;
  
  weight6[2] = 0.4679139345726910473898703;
  weight6[1] = 0.3607615730481386075698335;
  weight6[0] = 0.1713244923791703450402961;
  GQ6[2] = -0.2386191860831969086305017;
  GQ6[1] = -0.6612093864662645136613996;
  GQ6[0] = -0.9324695142031520278123016;

  weight4[1] = 0.6521451548625461;
  weight4[0] = 0.3478548451374538;
  GQ4[1] = -0.3399810435848563;
  GQ4[0] = -0.8611363115940526;

  GQ_loc_2[0] = -0.5773502691896257;
  weight2[0] = 1.0000000000000000;

  weight3[2] = 0.5555555555555556;
  weight3[1] = 0.8888888888888888;
  weight3[0] = 0.5555555555555556;
  GQ_loc_3[2] = 0.7745966692414834;
  GQ_loc_3[1] = 0.0000000000000000;
  GQ_loc_3[0] = -0.7745966692414834;


  GQ7[0] = -0.9491079123427585;
  GQ7[1] = -0.7415311855993945;
  GQ7[2] = -0.4058451513773972;
  GQ7[3] = 0.0000000000000000;
  GQ7[4] = 0.4058451513773972;
  GQ7[5] = 0.7415311855993945;
  GQ7[6] = 0.9491079123427585;
  weight7[0] = 0.1294849661688697;
  weight7[1] = 0.2797053914892766;
  weight7[2] = 0.3818300505051189;
  weight7[3] = 0.4179591836734694;
  weight7[4] = 0.3818300505051189;
  weight7[5] = 0.2797053914892766;
  weight7[6] = 0.1294849661688697;
 

  GQ5[0] = -0.9061798459386640;
  GQ5[1] = -0.5384693101056831;
  GQ5[2] = 0.0000000000000000;
  GQ5[3] = 0.5384693101056831;
  GQ5[4] = 0.9061798459386640;
  weight5[0] = 0.2369268850561891;
  weight5[1] = 0.4786286704993665;
  weight5[2] = 0.5688888888888889;
  weight5[3] = 0.4786286704993665;
  weight5[4] = 0.2369268850561891;
  
  
  //More Gaussian weights and nodes:
  for (int i = 1; i <= 10; i = i + 1)
    {
      weight20[i + 9] = weight20[10 - i];
      GQ20[i + 9] = -1 * GQ20[10 - i];
    }
  for (int i = 1; i <= 5; i = i + 1)
    {
      weight10[i + 4] = weight10[5 - i];
      GQ10[i + 4] = -1 * GQ10[5 - i];
    }
  for (int i = 1; i <= 3; i = i + 1)
    {
      weight6[i + 2] = weight6[3 - i];
      GQ6[i + 2] = -1 * GQ6[3 - i];
    }

for (int i = 1; i <= 2; i = i + 1)
    {
      weight4[i + 1] = weight4[2 - i];
      GQ4[i + 1] = -1 * GQ4[2 - i];
    }

for (int i = 1; i <= 1; i = i + 1)
    {
      weight2[i + 0] = weight2[1 - i];
      GQ_loc_2[i + 0] = -1 * GQ_loc_2[1 - i];
    }
  switch (dim)
    {
    case 2: 
      {
	output.loc = GQ_loc_2[spot];
	output.weight = weight2[spot];
	return output;
	break;
      }
    case 3: 
      {
	output.loc = GQ_loc_3[spot];
	output.weight = weight3[spot];
	return output;
	break;
      }
    case 4: 
      {
	output.loc = GQ4[spot];
	output.weight = weight4[spot];
	return output;
	break;
      }
    case 5: 
      {
	output.loc = GQ5[spot];
	output.weight = weight5[spot];
	return output;
	break;
      }
    case 6: 
      {
	output.loc = GQ6[spot];
	output.weight = weight6[spot];
	return output;
	break;
      }
      case 7: 
      {
	output.loc = GQ7[spot];
	output.weight = weight7[spot];
	return output;
	break;
      }
    case 10:
      {
	output.loc = GQ10[spot];
	output.weight = weight10[spot];
	return output;
	break;
      }
    case 20:
      {
	output.loc = GQ20[spot];
	output.weight = weight20[spot];
	return output;
	break;
      }
    }
}

//Some matrix algebra tools:
D2_double initialize_Iden(int d1, int d2)
{
  D2_double output = Allocate_2D_double(d1,d2);
  for (int row = 0; row < d1; row = row + 1)
    {
      for (int col = 0; col < d2; col = col + 1)
	{
	  output[row][col] = 0.0;
	}
      output[row][row] = 1.0;
    }
  return output;
}
void print_matrix(D2_double input, int d1, int d2)
{
  printf("\n");
  for (int row = 0; row < d1; row = row + 1)
    {
      printf("row %d:  ", row);
      for (int col = 0; col < d2; col = col + 1)
	{
	  printf("%6.3f,",input[row][col]);
	}
      printf("\n");
    }
}
D2_double invert_matrix_v1(D2_double input, int n)
{
  D2_double output = Allocate_2D_double(n,n);
  //Populate Appropriately sized Identity
  D2_double Iden_local = Allocate_2D_double(n,n);
  D2_double shift_Iden = Allocate_2D_double(n,n); //records preemptive row swithces
  D2_double MAT_local = Allocate_2D_double(n,n);
  int grab;
  double holder;
  double L1_cond[n];
  D2_double P_local = Allocate_2D_double(n,n);
  D2_double P_inverse_local = Allocate_2D_double(n,n);
  double sum;
  D2_double PMAT_local = Allocate_2D_double(n,n);
  double PIc_local[n];
  D2_double sys_mod_local = Allocate_2D_double(n,n+1);
  D2_double sys_ori_local = Allocate_2D_double(n,n+1);
  double MAT_inverse_col_local[n];

  double factor;
  int pivot_finder;
  double val_pivot;
  double b_ori[n];
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  Iden_local[row][col] = 0.0;
	}
      Iden_local[row][row] = 1.0;
    }
  /*printf("\tINPUT MATRIX:\n");
  for (int row = 0; row < n; row = row + 1)
    {
      printf("Row %d:\t",row);
      for (int col = 0; col < n; col = col + 1)
	{
	  printf("%5.3f,  ", input[row][col]);
	}
      printf("\n");
    }
  printf("\tIDENTITY MATRIX:\n");
  for (int row = 0; row < n; row = row + 1)
    {
      printf("Row %d:\t",row);
      for (int col = 0; col < n; col = col + 1)
	{
	  printf("%5.3f,  ", Iden_local[row][col]);
	}
      printf("\n");
      }*/
  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      //We work in one column of the full inversion at a time
      for (int row = 0; row < n; row = row + 1)
	{
	  shift_Iden[row][oseer] = Iden_local[row][oseer];
	}
      //Reset MAT to input matrix for peace of mind
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      MAT_local[row][col] = input[row][col];
	    }
	}
      //Perform preemptive row switch to prevent zeros in diagonal
      for (int row = 0; row < n; row = row + 1)
	{
	  grab = row;
	  while (fabs(MAT_local[grab][row]) <= mat_tol && grab < n)
	    {
	      // printf("Performed preemtive row switch for Cmat row %d\n",row);
	      grab = grab + 1;
	    }
	  if (grab == n) //we've reached the end of the of the column without finding nonzero
	    {
	      grab = row-1;
	      while (fabs(MAT_local[grab][row]) <= mat_tol)
		{
		  grab = grab-1;
		}
	      //when the while loop terminates, we have found a nonzero in the column
	    }
	  
	  //Make appropriate row switch in MAT and Iden_col
	  if (grab >= row) //Normal case. We did forward track to make row switch
	    {
	      for (int col = 0; col < n; col = col + 1)
		{
		  holder = MAT_local[row][col];
		  MAT_local[row][col] = MAT_local[grab][col];
		  MAT_local[grab][col] = holder;
		}

	      holder = shift_Iden[row][oseer];
	      shift_Iden[row][oseer] = shift_Iden[grab][oseer];
	      shift_Iden[grab][oseer] = holder;
	    }
	  if (grab < row) //we had to backtrack. Perform row arithmetic to populate diagonal
	    {
	      for (int col = 0; col < n; col = col + 1)
		{
		  MAT_local[row][col] = MAT_local[row][col] + MAT_local[grab][col];
		}
	      
	      shift_Iden[row][oseer]  = shift_Iden[row][oseer] + shift_Iden[grab][oseer];
	    }
	}
    }
      
  /*
    FILE*file99;
    file99 = fopen("switched_f_LHS.csv","w");
    for (int row = 0; row < n; row = row + 1)
    {
    for (int col = 0; col < n; col = col + 1)
    {
    fprintf(file99, "%30.16f,",MAT_local[row][col]);
    }
    fprintf(file99, "\n");
    }
    fclose(file99);
  */
  
  //Now, set preconditioner P based on MAT------------------------------------
  //my_Precon(1); //P and P_inverse are now properly populated
  
  //Preconditioning: move PMAT closer to diagonal dominance than MAT
  for (int col = 0; col < n; col = col + 1)
    {
      L1_cond[col] = 0;
      for (int row = 0; row < n; row = row + 1)
	{
	  if (row != col)
	    {
	      L1_cond[col] = L1_cond[col] + fabs(MAT_local[row][col]);
	    }
	}
      L1_cond[col] = L1_cond[col] / fabs(MAT_local[col][col]);
    }
  //L1_cond communicates the ratio between a column's entries and its diagonal
  for (int row = 0; row < n; row = row + 1)
    {
      //Set diagonal value to be exaggerated in PA
      P_local[row][row] = cond_mult * L1_cond[row] / MAT_local[row][row];
      P_inverse_local[row][row] = MAT_local[row][row] / (cond_mult*L1_cond[row]);
      // P_local[row][row] = 1;
      //P_inverse_local[row][row] = 1;
      for (int col = 0; col < n; col = col + 1)
	{
	  if (row != col)
	    {
	      P_local[row][col] = 0;
	      P_inverse_local[row][col] = 0;
	    }
	}
    }
  //variables P and P_inverse have been set for the current problem-------------
  //Set PMAT and P(Iden_col) as system inputs, populate sys_mod--------
  //populate_sys(1);
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * MAT_local[k][col];
	    }
	  PMAT_local[row][col] = sum;
	}
    }
  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      for (int row = 0; row < n; row = row + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * shift_Iden[k][oseer];
	    }
	  PIc_local[row] = sum;
	}
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sys_mod_local[row][col] = PMAT_local[row][col];
	    }
	  sys_mod_local[row][n] = PIc_local[row];
	}
      //done-----------------------------------------------------------
      //function takes sys_mod, returns solution x for PMx = PIc------
      //vec_sys_solve(1);
      
      
      //The term ''pivot'' identifies the column number where the pivot occurs
      for (int pivot = 0; pivot < n; pivot = pivot + 1)
	{
	  //First, find row with highest magnitude in pivot column, switch it to pivot row
	  pivot_finder = pivot;
	  val_pivot = sys_mod_local[pivot][pivot];
	  for (int scroll = pivot+1; scroll < n; scroll = scroll+1)
	    {
	      if (fabs(sys_mod_local[scroll][pivot]) > fabs(sys_mod_local[pivot_finder][pivot]))
		{
		  pivot_finder = scroll;
		}
	    }
	  //pivot_finder holds the row of the largest value in the pivot column.
	  //Perform row switch for entire system:
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      //store value that is about to be replaced
	      holder = sys_mod_local[pivot][col];
	      //move appropriate row to pivot position
	      sys_mod_local[pivot][col] = sys_mod_local[pivot_finder][col];
	      //store previous [pivot][col] value in void from the row that was moved up
	      sys_mod_local[pivot_finder][col] = holder;
	    }
	  val_pivot = sys_mod_local[pivot][pivot]; //value at pivot position
	  //Divide pivot row by pivot value to get 1 in pivot position
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      sys_mod_local[pivot][col] = sys_mod_local[pivot][col] / val_pivot;
	    }
	  //Run through all rows below pivot, get zeros in pivot column
	  if (pivot < (n-1))
	    {
	      for (int row = pivot+1; row < n; row = row + 1)
		{
		  factor = sys_mod_local[row][pivot];
		  //Reduce pivot column to zeros outside of pivot row
		  for (int col = 0; col < n+1; col = col + 1)
		    {
		      sys_mod_local[row][col] = sys_mod_local[row][col] - sys_mod_local[pivot][col] * factor;
		    }
		}
	    }
	}
      //sys_mod now holds the sloution for PMx = PIc.
      //Done with vec_solve routine----------------------------------
      //Extract solution to Mx=Ic from sys_ori
      //vec_sys_extract(1); //populated sys_ori with solution to Mx=Ic
      //Myltiply bothy sides of PMx=PIc system by P^-1
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = 0;
	      for (int k = 0; k < n; k = k + 1)
		{
		  sum = sum + P_inverse_local[row][k] * sys_mod_local[k][col];
		}
	      sys_ori_local[row][col] = sum;
	    }
	}
      for (int k = 0; k < n; k = k + 1)
	{
	  sum = 0;
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = sum + P_inverse_local[k][col] * sys_mod_local[col][n];
	    }
	  sys_ori_local[k][n] = sum;
	}
      //Done-----------------------------------------------------------
      //Use sys_ori values to populate MAT_inverse_col :----------------
      //populate_Minverse(1);
      
      
      //Get b from right column of sys
      for (int row = 0; row < n; row = row + 1)
	{
	  b_ori[row] = sys_ori_local[row][n];
	}
      //sys_ori and b_ori now form the Mx=b system of interest
      //M is in reduced echelon form, but with disorganized diagonal
      //NOw, back substitute sys values to get x vector
      for (int row = n-1; row >= 0; row = row - 1)
	{
	  sum = 0;
	  for (int col = row+1; col < n; col = col + 1)
	    {
	      sum = sum + sys_ori_local[row][col]*MAT_inverse_col_local[col];
	    }
	  MAT_inverse_col_local[row] = (b_ori[row] - sum) / sys_ori_local[row][row];
	}
      //Done Populating relevant row of inverse----------------
      for (int row = 0; row < n; row = row + 1)
	{
	  output[row][oseer] = MAT_inverse_col_local[row];
	}
    }
  return output;
} 
D2_double invert_matrix_v2(D2_double input, int n)
{
  D2_double output = Allocate_2D_double(n,n);
  //Populate Appropriately sized Identity
  D2_double Iden_local = Allocate_2D_double(n,n);
  D2_double shift_Iden = Allocate_2D_double(n,n); //records preemptive row swithces
  D2_double MAT_local = Allocate_2D_double(n,n);
  int grab;
  double holder;
  double L1_cond[n];
  D2_double P_local = Allocate_2D_double(n,n);
  D2_double P_inverse_local = Allocate_2D_double(n,n);
  double sum;
  D2_double PMAT_local = Allocate_2D_double(n,n);
  double PIc_local[n];
  D2_double sys_mod_local = Allocate_2D_double(n,n+1);
  D2_double sys_ori_local = Allocate_2D_double(n,n+1);
  double MAT_inverse_col_local[n];

  double factor;
  int pivot_finder;
  double val_pivot;
  double b_ori[n];
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  Iden_local[row][col] = 0.0;
	}
      Iden_local[row][row] = 1.0;
    }

  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      //We work in one column of the full inversion at a time
      for (int row = 0; row < n; row = row + 1)
	{
	  shift_Iden[row][oseer] = Iden_local[row][oseer];
	}
      //Reset MAT to input matrix for peace of mind
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      MAT_local[row][col] = input[row][col];
	    }
	}
      //Perform preemptive row switch to prevent zeros in diagonal
      for (int row = 0; row < n; row = row + 1)
	{
	  grab = row;
	  while (fabs(MAT_local[grab][row]) <= mat_tol && grab < (n-1))
	    {
	      //printf("Found need for preemptive row switch for Cmat row %d\n",row);
	      grab = grab + 1;
	      //printf("grab = %d\n", grab);
	    }
	  
	  if (grab == n) //we've reached the end of the of the column without finding nonzero
	    {
	      //printf("Failed to find a nonzero entry going down the column\n");
	      grab = row-1;
	      while (fabs(MAT_local[grab][row]) <= mat_tol)
		{
		  grab = grab-1;
		}
	      //when the while loop terminates, we have found a nonzero in the column
	    }
	  
	  //Make appropriate row switch in MAT and Iden_col
	  if (grab >= row) //Normal case. We did forward track to make row switch
	    {
	      //printf("Performed forward/no row switch\n");
	      for (int col = 0; col < n; col = col + 1)
		{
		  holder = MAT_local[row][col];
		  MAT_local[row][col] = MAT_local[grab][col];
		  MAT_local[grab][col] = holder;
		}

	      holder = shift_Iden[row][oseer];
	      shift_Iden[row][oseer] = shift_Iden[grab][oseer];
	      shift_Iden[grab][oseer] = holder;
	    }
	  if (grab < row) //we had to backtrack. Perform row arithmetic to populate diagonal
	    {
	      //printf("Performed rearward row switch\n");
	      for (int col = 0; col < n; col = col + 1)
		{
		  MAT_local[row][col] = MAT_local[row][col] + MAT_local[grab][col];
		}
	      
	      shift_Iden[row][oseer]  = shift_Iden[row][oseer] + shift_Iden[grab][oseer];
	    }
	}
    }
      
  /*
    FILE*file99;
    file99 = fopen("switched_f_LHS.csv","w");
    for (int row = 0; row < n; row = row + 1)
    {
    for (int col = 0; col < n; col = col + 1)
    {
    fprintf(file99, "%30.16f,",MAT_local[row][col]);
    }
    fprintf(file99, "\n");
    }
    fclose(file99);
  */
  
  //Now, set preconditioner P based on MAT------------------------------------
  //my_Precon(1); //P and P_inverse are now properly populated
  
  //Preconditioning: move PMAT closer to diagonal dominance than MAT
  for (int col = 0; col < n; col = col + 1)
    {
      L1_cond[col] = 0;
      for (int row = 0; row < n; row = row + 1)
	{
	  if (row != col)
	    {
	      L1_cond[col] = L1_cond[col] + fabs(MAT_local[row][col]);
	    }
	}
      L1_cond[col] = L1_cond[col] / fabs(MAT_local[col][col]);
    }
  //L1_cond communicates the ratio between a column's entries and its diagonal
  for (int row = 0; row < n; row = row + 1)
    {
      //Set diagonal value to be exaggerated in PA
      P_local[row][row] = cond_mult * L1_cond[row] / MAT_local[row][row];
      P_inverse_local[row][row] = MAT_local[row][row] / (cond_mult*L1_cond[row]);
      // P_local[row][row] = 1;
      //P_inverse_local[row][row] = 1;
      for (int col = 0; col < n; col = col + 1)
	{
	  if (row != col)
	    {
	      P_local[row][col] = 0;
	      P_inverse_local[row][col] = 0;
	    }
	}
    }
  //variables P and P_inverse have been set for the current problem-------------
  //Set PMAT and P(Iden_col) as system inputs, populate sys_mod--------
  //populate_sys(1);
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * MAT_local[k][col];
	    }
	  PMAT_local[row][col] = sum;
	}
    }
  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      for (int row = 0; row < n; row = row + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * shift_Iden[k][oseer];
	    }
	  PIc_local[row] = sum;
	}
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sys_mod_local[row][col] = PMAT_local[row][col];
	    }
	  sys_mod_local[row][n] = PIc_local[row];
	}
      //done-----------------------------------------------------------
      //function takes sys_mod, returns solution x for PMx = PIc------
      //vec_sys_solve(1);
      
      
      //The term ''pivot'' identifies the column number where the pivot occurs
      for (int pivot = 0; pivot < n; pivot = pivot + 1)
	{
	  //First, find row with highest magnitude in pivot column, switch it to pivot row
	  pivot_finder = pivot;
	  val_pivot = sys_mod_local[pivot][pivot];
	  for (int scroll = pivot+1; scroll < n; scroll = scroll+1)
	    {
	      if (fabs(sys_mod_local[scroll][pivot]) > fabs(sys_mod_local[pivot_finder][pivot]))
		{
		  pivot_finder = scroll;
		}
	    }
	  //pivot_finder holds the row of the largest value in the pivot column.
	  //Perform row switch for entire system:
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      //store value that is about to be replaced
	      holder = sys_mod_local[pivot][col];
	      //move appropriate row to pivot position
	      sys_mod_local[pivot][col] = sys_mod_local[pivot_finder][col];
	      //store previous [pivot][col] value in void from the row that was moved up
	      sys_mod_local[pivot_finder][col] = holder;
	    }
	  val_pivot = sys_mod_local[pivot][pivot]; //value at pivot position
	  //Divide pivot row by pivot value to get 1 in pivot position
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      sys_mod_local[pivot][col] = sys_mod_local[pivot][col] / val_pivot;
	    }
	  //Run through all rows below pivot, get zeros in pivot column
	  if (pivot < (n-1))
	    {
	      for (int row = pivot+1; row < n; row = row + 1)
		{
		  factor = sys_mod_local[row][pivot];
		  //Reduce pivot column to zeros outside of pivot row
		  for (int col = 0; col < n+1; col = col + 1)
		    {
		      sys_mod_local[row][col] = sys_mod_local[row][col] - sys_mod_local[pivot][col] * factor;
		    }
		}
	    }
	}
      //sys_mod now holds the sloution for PMx = PIc.
      //Done with vec_solve routine----------------------------------
      //Extract solution to Mx=Ic from sys_ori
      //vec_sys_extract(1); //populated sys_ori with solution to Mx=Ic
      //Myltiply bothy sides of PMx=PIc system by P^-1
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = 0;
	      for (int k = 0; k < n; k = k + 1)
		{
		  sum = sum + P_inverse_local[row][k] * sys_mod_local[k][col];
		}
	      sys_ori_local[row][col] = sum;
	    }
	}
      for (int k = 0; k < n; k = k + 1)
	{
	  sum = 0;
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = sum + P_inverse_local[k][col] * sys_mod_local[col][n];
	    }
	  sys_ori_local[k][n] = sum;
	}
      //Done-----------------------------------------------------------
      //Use sys_ori values to populate MAT_inverse_col :----------------
      //populate_Minverse(1);
      
      
      //Get b from right column of sys
      for (int row = 0; row < n; row = row + 1)
	{
	  b_ori[row] = sys_ori_local[row][n];
	}
      //sys_ori and b_ori now form the Mx=b system of interest
      //M is in reduced echelon form, but with disorganized diagonal
      //NOw, back substitute sys values to get x vector
      for (int row = n-1; row >= 0; row = row - 1)
	{
	  sum = 0;
	  for (int col = row+1; col < n; col = col + 1)
	    {
	      sum = sum + sys_ori_local[row][col]*MAT_inverse_col_local[col];
	    }
	  MAT_inverse_col_local[row] = (b_ori[row] - sum) / sys_ori_local[row][row];
	}
      //Done Populating relevant row of inverse----------------
      for (int row = 0; row < n; row = row + 1)
	{
	  output[row][oseer] = MAT_inverse_col_local[row];
	}
    }
  return output;
} 
D2_double invert_matrix_v3(D2_double input, int n)
{
  //This is the "invert_matrix_better" algorithm from MAT.cpp in RDG_2D_NewBC
  D2_double output = Allocate_2D_double(n,n);
  //Populate Appropriately sized Identity
  D2_double Iden_local = Allocate_2D_double(n,n);
  D2_double shift_Iden = Allocate_2D_double(n,n); //records preemptive row swithces
  D2_double MAT_local = Allocate_2D_double(n,n);
  int grab;
  double holder;
  double L1_cond[n];
  D2_double P_local = Allocate_2D_double(n,n);
  D2_double P_inverse_local = Allocate_2D_double(n,n);
  double sum;
  D2_double PMAT_local = Allocate_2D_double(n,n);
  double PIc_local[n];
  D2_double sys_mod_local = Allocate_2D_double(n,n+1);
  D2_double sys_ori_local = Allocate_2D_double(n,n+1);
  double MAT_inverse_col_local[n];

  double factor;
  int pivot_finder;
  double val_pivot;
  double b_ori[n];
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  Iden_local[row][col] = 0.0;
	}
      Iden_local[row][row] = 1.0;
    }
  
  double test_value;
  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      //We work in one column of the full inversion at a time
      for (int row = 0; row < n; row = row + 1)
	{
	  shift_Iden[row][oseer] = Iden_local[row][oseer];
	}
      //Reset MAT to input matrix for peace of mind
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      MAT_local[row][col] = input[row][col];
	    }
	}
      //Perform preemptive row switch to prevent zeros in diagonal
      //Can't allow the subloop to attempt to access MAT_local[n][row]
      for (int row = 0; row < n; row = row + 1)
	{
	  grab = row;
	  test_value = MAT_local[grab][row];
	  while (fabs(test_value) <= mat_tol && grab < n)
	    {
	      grab = grab + 1;
	      if (grab == n)
		{
		  test_value = 0.0;
		}
	      else
		{
		  test_value = MAT_local[grab][row];
		}
	    }
	  
	  if (grab == n) //we've reached the end of the of the column without finding nonzero
	    {
	      grab = row-1;
	      while (fabs(MAT_local[grab][row]) <= mat_tol)
		{
		  grab = grab-1;
		}
	      //when the while loop terminates, we have found a nonzero in the column
	    }
	  
	  //Make appropriate row switch in MAT and Iden_col
	  if (grab >= row) //Normal case. We did forward track to make row switch
	    {
	      //printf("tracked forward to make row switch\n");
	      for (int col = 0; col < n; col = col + 1)
		{
		  holder = MAT_local[row][col];
		  MAT_local[row][col] = MAT_local[grab][col];
		  MAT_local[grab][col] = holder;
		}
	      
	      holder = shift_Iden[row][oseer];
	      shift_Iden[row][oseer] = shift_Iden[grab][oseer];
	      shift_Iden[grab][oseer] = holder;
	    }
	  if (grab < row) //we had to backtrack. Perform row arithmetic to populate diagonal
	    {						
	      //printf("tracked backwards to make row switch\n");
	      for (int col = 0; col < n; col = col + 1)
		{
		  MAT_local[row][col] = MAT_local[row][col] + MAT_local[grab][col];
		}
	      shift_Iden[row][oseer]  = shift_Iden[row][oseer] + shift_Iden[grab][oseer];
	    }
	}
    }
  
  //Now, set preconditioner P based on MAT------------------------------------
  //my_Precon(1); //P and P_inverse are now properly populated
  
  //Preconditioning: move PMAT closer to diagonal dominance than MAT
  for (int col = 0; col < n; col = col + 1)
    {
      L1_cond[col] = 0;
      for (int row = 0; row < n; row = row + 1)
	{
	  if (row != col)
	    {
	      L1_cond[col] = L1_cond[col] + fabs(MAT_local[row][col]);
	    }
	}
      L1_cond[col] = L1_cond[col] / fabs(MAT_local[col][col]);
    }
  //L1_cond communicates the ratio between a column's entries and its diagonal
  for (int row = 0; row < n; row = row + 1)
    {
      //Set diagonal value to be exaggerated in PA
      P_local[row][row] = cond_mult * L1_cond[row] / MAT_local[row][row];
      P_inverse_local[row][row] = MAT_local[row][row] / (cond_mult*L1_cond[row]);
      // P_local[row][row] = 1;
      //P_inverse_local[row][row] = 1;
      for (int col = 0; col < n; col = col + 1)
	{
	  if (row != col)
	    {
	      P_local[row][col] = 0;
	      P_inverse_local[row][col] = 0;
	    }
	}
    }

  //To cambat poor performance: Setting P=P^-1=I
  P_local = initialize_Iden(n,n);
  P_inverse_local = initialize_Iden(n,n);


  //variables P and P_inverse have been set for the current problem-------------
  //Set PMAT and P(Iden_col) as system inputs, populate sys_mod--------
  //populate_sys(1);
  
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * MAT_local[k][col];
	    }
	  PMAT_local[row][col] = sum;
	}
    }
  for (int oseer = 0; oseer < n; oseer = oseer + 1)
    {
      for (int row = 0; row < n; row = row + 1)
	{
	  sum = 0;
	  for (int k = 0; k < n; k = k + 1)
	    {
	      sum = sum + P_local[row][k] * shift_Iden[k][oseer];
	    }
	  PIc_local[row] = sum;
	}
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sys_mod_local[row][col] = PMAT_local[row][col];
	    }
	  sys_mod_local[row][n] = PIc_local[row];
	}
      //done-----------------------------------------------------------
      //function takes sys_mod, returns solution x for PMx = PIc------
      //vec_sys_solve(1);
      
      
      //The term ''pivot'' identifies the column number where the pivot occurs
      for (int pivot = 0; pivot < n; pivot = pivot + 1)
	{
	  //First, find row with highest magnitude in pivot column, switch it to pivot row
	  pivot_finder = pivot;
	  val_pivot = sys_mod_local[pivot][pivot];
	  for (int scroll = pivot+1; scroll < n; scroll = scroll+1)
	    {
	      if (fabs(sys_mod_local[scroll][pivot]) > fabs(sys_mod_local[pivot_finder][pivot]))
		{
		  pivot_finder = scroll;
		}
	    }
	  //pivot_finder holds the row of the largest value in the pivot column.
	  //Perform row switch for entire system:
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      //store value that is about to be replaced
	      holder = sys_mod_local[pivot][col];
	      //move appropriate row to pivot position
	      sys_mod_local[pivot][col] = sys_mod_local[pivot_finder][col];
	      //store previous [pivot][col] value in void from the row that was moved up
	      sys_mod_local[pivot_finder][col] = holder;
	    }
	  val_pivot = sys_mod_local[pivot][pivot]; //value at pivot position
	  //Divide pivot row by pivot value to get 1 in pivot position
	  for (int col = 0; col < n+1; col = col + 1)
	    {
	      sys_mod_local[pivot][col] = sys_mod_local[pivot][col] / val_pivot;
	    }
	  //Run through all rows below pivot, get zeros in pivot column
	  if (pivot < (n-1))
	    {
	      for (int row = pivot+1; row < n; row = row + 1)
		{
		  factor = sys_mod_local[row][pivot];
		  //Reduce pivot column to zeros outside of pivot row
		  for (int col = 0; col < n+1; col = col + 1)
		    {
		      sys_mod_local[row][col] = sys_mod_local[row][col] - sys_mod_local[pivot][col] * factor;
		    }
		}
	    }
	}
      //sys_mod now holds the sloution for PMx = PIc.
      //Done with vec_solve routine----------------------------------
      //Extract solution to Mx=Ic from sys_ori
      //vec_sys_extract(1); //populated sys_ori with solution to Mx=Ic
      //Myltiply bothy sides of PMx=PIc system by P^-1
      for (int row = 0; row < n; row = row + 1)
	{
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = 0;
	      for (int k = 0; k < n; k = k + 1)
		{
		  sum = sum + P_inverse_local[row][k] * sys_mod_local[k][col];
		}
	      sys_ori_local[row][col] = sum;
	    }
	}
      for (int k = 0; k < n; k = k + 1)
	{
	  sum = 0;
	  for (int col = 0; col < n; col = col + 1)
	    {
	      sum = sum + P_inverse_local[k][col] * sys_mod_local[col][n];
	    }
	  sys_ori_local[k][n] = sum;
	}
      //Done-----------------------------------------------------------
      //Use sys_ori values to populate MAT_inverse_col :----------------
      //populate_Minverse(1);
      
      
      //Get b from right column of sys
      for (int row = 0; row < n; row = row + 1)
	{
	  b_ori[row] = sys_ori_local[row][n];
	}
      //sys_ori and b_ori now form the Mx=b system of interest
      //M is in reduced echelon form, but with disorganized diagonal
      //NOw, back substitute sys values to get x vector
      for (int row = n-1; row >= 0; row = row - 1)
	{
	  sum = 0;
	  for (int col = row+1; col < n; col = col + 1)
	    {
	      sum = sum + sys_ori_local[row][col]*MAT_inverse_col_local[col];
	    }
	  MAT_inverse_col_local[row] = (b_ori[row] - sum) / sys_ori_local[row][row];
	}
      //Done Populating relevant row of inverse----------------
      for (int row = 0; row < n; row = row + 1)
	{
	  output[row][oseer] = MAT_inverse_col_local[row];
	}
    }
  
  return output;
}
D2_double invert_matrix(D2_double input, int n)
{
  return invert_matrix_v3(input, n);
}
D2_double clean_inverse(D2_double input, int n)
{
  D2_double input_inverse = Allocate_2D_double(n,n);
  input_inverse = invert_matrix(input,n); //get dirty inverse of input
  //I'm not sure if this cleaning process does any good
  int master = 0;
  double MAT_inverse_col_local[n];
  double Iden_local[n][n];
  double Iden_col_local[n];
  double sum;
  double vec2[n];
  double vec3[n];
  double vec_diff[n];
  for (int row = 0; row < n; row = row + 1)
    {
      for (int col = 0; col < n; col = col + 1)
	{
	  Iden_local[row][col] = 0.0;
	}
      Iden_local[row][row] = 1.0;
    }
  while (master < master_max)
    {
      master = master + 1;
      for (int col = 0; col < n; col = col + 1)
	{
	  for (int row = 0; row < n; row = row + 1)
	    {
	      MAT_inverse_col_local[row] = input_inverse[row][col];
	      //MAT_inverse_col[row] = Bmat_inverse[row][col];
	      Iden_col_local[row] = Iden_local[row][col];
	    }
	  for (int row = 0; row < n; row = row + 1)
	    {
	      sum = 0;
	      for (int col2 = 0; col2 < n; col2 = col2 + 1)
		{
		  sum = sum + input[row][col2] * MAT_inverse_col_local[col2];
		  //sum = sum + Bmat[row][col2] * MAT_inverse_col[col2];
		}
	      //B*x_old product
	      vec2[row] = sum;
	    }
	  for (int row = 0; row < n; row = row + 1)
	    {
	      vec_diff[row] = Iden_col_local[row] - vec2[row];
	    }
	  for (int row = 0; row < n; row = row + 1)
	    {
	      sum = 0;
	      for (int col2 = 0; col2 < n; col2 = col2 + 1)
		{
		  sum = sum + input_inverse[row][col2] * vec_diff[col2];
		  //sum = sum + Bmat_inverse[row][col2] * vec_diff[col2];
		}
	      vec3[row] = sum;
	    }
	  for (int row = 0; row < n; row = row + 1)
	    {
	      MAT_inverse_col_local[row] = MAT_inverse_col_local[row] + vec3[row];
	    }
	  for (int row = 0; row < n; row = row + 1)
	    {
	      input_inverse[row][col] = MAT_inverse_col_local[row];
	    }
	}
    }
  return input_inverse;
}
double mat_condition(D2_double A, int d1, int d2) //return condition number
{
  //Condition number of d1*d2 matrix
  double output = 0;
  double sums[d1];
  for (int row = 0; row < d1; row = row + 1)
    {
      double sum = 0;
      for (int col = 0; col < d2; col = col + 1)
	{
	  sum = sum + fabs(A[row][col]);
	}
      sums[row] = sum;
    }
  double row_min = sums[0];
  double row_max = sums[0];
  for (int row = 1; row < d1; row = row + 1)
    {
      row_max = fmax(row_max, sums[row]);
      row_min = fmin(row_min, sums[row]);
    }
  
  return row_max / row_min;
}
D1_double mat_by_col(D2_double A, D1_double b, int d_row, int d_col)
{
  D1_double output = Allocate_1D_double(d_col);
  double sum;
  for (int row = 0; row < d_row; row = row + 1)
    {
      sum = 0;
      for (int col = 0; col < d_col; col = col + 1)
	{
	  sum = sum + A[row][col]*b[col];
	}
      output[row] = sum;
    }
  return output;
}
D2_double mat_by_mat(D2_double A, D2_double B, int A_row, int A_col, int B_col)
{
  double sum;
  D2_double output = Allocate_2D_double(A_row,B_col);
  for (int row = 0; row < A_row; row = row + 1)
    {
      for (int col = 0; col < B_col; col = col + 1)
	{
	  sum = 0;
	  for (int k = 0; k < A_col; k = k + 1)
	    {
	      sum = sum + A[row][k] * B[k][col];
	    }
	  output[row][col] = sum;
	}
    }
  return output;
}
//SOme meshing stuff


//Some meshing algorithms
D1_cord get_mesh_tri(int dummy) //picks up the vertex locations from mesh file
{
  
  /*
  FILE*msh;
  msh = fopen("2d_tri_grid0.msh",'r');
  fclose(msh);
  */
  /*
  D1_cord MESH_node; //these will be mesh's node locations
  return D1_cord;
  */
  D1_cord _nodes;
  //the filename passed to this function must match the desired mesh
  //make sure to match it with element, node count in consts section (top of code)
  //std::ifstream input (fileName);
  //std::ifstream input ("2d_tri_grid0.msh");
  std::ifstream input (mesh_file);
  std::string line;
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Nodes")
    printf("invalid file format, line != $Nodes\n");
  int nnodes;
  input>>nnodes;
  printf("From get_mesh_tri: nnodes = %d\n",nnodes);
  RES_nodes = nnodes;
  //_nodes.resize (3, nnodes);
  _nodes = Allocate_1D_cord(nnodes);
  for (int i = 0; i < nnodes; i++) {
    int nodeId;
    if (!(input >> nodeId >> _nodes[i].rx/*_nodes (0, i)*/ >> _nodes[i].ry/*_nodes (1, i)*/ >> _nodes[i].rz/*_nodes (2, i)*/))
      printf("invalid file format\n");
  }
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndNodes")
    printf("invalid file format, line != $EndNodes\n");
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Elements")
    printf("invalid file format, line != $Elements\n");

  return _nodes;
  /*
  int nelements;
  input >> nelements;
  std::vector<int> enodes;
  _elements.resize(MSH_NUM_TYPE);      // holds elements in my partition
  _otherElements.resize(MSH_NUM_TYPE); // holds elements in other partitions
  getline (input, line);
  for (int i = 0; i < nelements; i++) {
    enodes.resize(0);
    int elementId, elementType, ntags, ptag, num, partition=1;
    getline (input, line);
    std::istringstream sline (line);
    sline >> elementId >> elementType >> ntags; 
    for (int j = 0; j < ntags; j++) {
      int tag;
      if      (j==0) sline >> ptag;      // physical tag of element
      else if (j==3) sline >> partition; // main partition of element
      else           sline >> tag;
      //printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);
    }
    int idNode;
    while (sline >> idNode) {
      enodes.push_back(idNode-1);
    }
    // Exit if the partition number is larger than the total number of processors
    if(_numprocs < partition){printf("Not enough processors for the mesh partitions. Exiting\n"); exit(1);}
    // Only store the element if it's in my partition
    if(_myid == partition-1) _elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
    // Otherwise store the element in the otherElements
    else _otherElements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
  }
  getline (input, line);
  if (line!="$EndElements")
    printf("invalid file format\n");
  */
}
D2_integer get_c2n_tri_SAM(int dummy) //pairs elements to the vertex locations
{
  //this actually repeats a good deal of the get_mesh procedure,
  //for future modularity i guess
  //std::ifstream input (fileName);
  //std::ifstream input ("2d_tri_grid0.msh");
  std::ifstream input (mesh_file);
  std::string line;
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Nodes")
    printf("invalid file format, line != $Nodes\n");
  int nnodes;
  input>>nnodes;
  //_nodes.resize (3, nnodes);
  D1_cord _nodes = Allocate_1D_cord(nnodes);
  for (int i = 0; i < nnodes; i++) {
    int nodeId;
    if (!(input >> nodeId >> _nodes[i].rx/*_nodes (0, i)*/ >> _nodes[i].ry/*_nodes (1, i)*/ >> _nodes[i].rz/*_nodes (2, i)*/))
      printf("invalid file format\n");
  }
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndNodes")
    printf("invalid file format, line != $EndNodes\n");
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Elements")
    printf("invalid file format, line != $Elements\n");

  
  int nelements;
  input >> nelements;
  printf("From get_c2n_tri: nelements = %d\n",nelements);
  printf("That mat be an error, depending on organization of gmsh file\n");
  //RES_cells = nelements;
  std::vector<int> enodes;
  //_elements.resize(MSH_NUM_TYPE);      // holds elements in my partition
  D2_integer elements = Allocate_2D_integer(nelements, N_N); //N_N nodes per element
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  int real_element = -1;
  for (int i = 0; i < nelements; i++) 
    {
      //   printf("\ni=%d, real_element = %d\n",i,real_element);
      enodes.resize(0); //erases everything in enodes
      int elementId, elementType, ntags, ptag, num, partition=1;
      getline (input, line);
      //    printf("present line is: \t%s\n", line.c_str());
      std::istringstream sline (line);
      sline >> elementId >> elementType >> ntags; 
      for (int j = 0; j < ntags; j++) {
	int tag;
	if      (j==0) sline >> ptag;      // physical tag of element
	else if (j==3) sline >> partition; // main partition of element
	else           sline >> tag;
	//printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);
      }
      int idNode;
      while (sline >> idNode) {
	enodes.push_back(idNode-1);
      }
      
      //_elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
      /*
	printf("enodes: ");
      for (int j = 0; j < enodes.size(); j++)
	{
	  printf("%d, ",enodes[j]);
	}
      printf("\n");
      */
      
      if (enodes.size() == Nsolut)
	{
	  //	  printf("%d nodes, so this is an actual element. Indexing real_element to %d\n", N_N,real_element+1);
	  real_element = real_element+1;
	  elements[real_element] = enodes;
	}
    }
  nelements = real_element+1;
  RES_cells = nelements;
  printf("Just set element count at RES_cells = %d\n", RES_cells);
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndElements")
    printf("invalid file format, line != $EndElements\n");
  return elements;
  
} 
D2_integer get_c2n_tri(int dummy) //pairs elements to the vertex locations
{
  //this actually repeats a good deal of the get_mesh procedure,
  //for future modularity i guess
  //std::ifstream input (fileName);
  //std::ifstream input ("2d_tri_grid0.msh");
  std::ifstream input (mesh_file);
  std::string line;
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Nodes")
    printf("invalid file format, line != $Nodes\n");
  int nnodes;
  input>>nnodes;
  //_nodes.resize (3, nnodes);
  D1_cord _nodes = Allocate_1D_cord(nnodes);
  for (int i = 0; i < nnodes; i++) {
    int nodeId;
    if (!(input >> nodeId >> _nodes[i].rx/*_nodes (0, i)*/ >> _nodes[i].ry/*_nodes (1, i)*/ >> _nodes[i].rz/*_nodes (2, i)*/))
      printf("invalid file format\n");
  }
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndNodes")
    printf("invalid file format, line != $EndNodes\n");
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Elements")
    printf("invalid file format, line != $Elements\n");

  
  int nelements;
  input >> nelements;
  printf("From get_c2n_tri: nelements = %d\n",nelements);
  RES_cells = nelements;
  std::vector<int> enodes;
  //_elements.resize(MSH_NUM_TYPE);      // holds elements in my partition
  D2_integer elements = Allocate_2D_integer(nelements, 3); //3 nodes per element
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  for (int i = 0; i < nelements; i++) 
    {
      enodes.resize(0); //erases everything in enodes
      int elementId, elementType, ntags, ptag, num, partition=1;
      getline (input, line);
      std::istringstream sline (line);
      sline >> elementId >> elementType >> ntags; 
      for (int j = 0; j < ntags; j++) {
	int tag;
	if      (j==0) sline >> ptag;      // physical tag of element
	else if (j==3) sline >> partition; // main partition of element
	else           sline >> tag;
	//printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);
      }
      int idNode;
      while (sline >> idNode) {
	enodes.push_back(idNode-1);
	int transfer = enodes.size();
	printf("i=%d: enodes size = %d\n", i,transfer);
      }
      
      //_elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
      elements[i] = enodes;
      
    }
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndElements")
    printf("invalid file format, line != $EndElements\n");
  return elements;
  
} 
D2_double grab_field_from_GMSH(int dummy) //pairs elements to the vertex locations
{
  //this actually repeats a good deal of the get_mesh procedure,
  //for future modularity i guess
  //std::ifstream input (fileName);
  //std::ifstream input ("2d_tri_grid0.msh");
  std::ifstream input (ux_file);
  std::string line;
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());


  D2_double _output = Allocate_2D_double(RES_cells, Nsolut); //Nsolut Lagrange nodes per element
  D1_double enodes;
  for (int om = 0; om < RES_cells; om++)
    {
      double messenger;
      //   printf("\ni=%d, real_element = %d\n",i,real_element);
      enodes.resize(0); //erases everything in enodes
      int elementId, elementType, ntags, ptag, num, partition=1;
      getline (input, line);
      printf("present line is: \t%s\n",line.c_str());
      printf("present line is: \t%f\n",atof(line.c_str()));
      std::istringstream sline (line);
      for (int j = 0; j < Nsolut+2; j++)
	{
	  int tag_discard;
	  if (j == 0)
	    {
	      sline >> tag_discard;
	    }
	  if (j == 1)
	    {
	      sline >> tag_discard;
	    }

	  else
	    {
	      sline >> enodes[j-2];
	    }
	  	  printf("j=%d, tag_discard=%d\n",j,tag_discard);
	}
	      
      
   
      sline >> elementId >> elementType >> ntags; 
      for (int j = 0; j < ntags; j++) {
	int tag;
	if      (j==0) sline >> ptag;      // physical tag of element
	else if (j==3) sline >> partition; // main partition of element
	else           sline >> tag;
	//printf("Id=%i: j=%i, tag=%i, ptag=%i, partition=%i \n",elementId, j, tag,ptag, partition);
      }
      int idNode;
      while (sline >> idNode) {
	enodes.push_back(idNode-1);
      }
      
      //_elements[elementType].push_back (simpleElement(elementId, ptag, partition-1, enodes));
      /*
	printf("enodes: ");
      for (int j = 0; j < enodes.size(); j++)
	{
	  printf("%d, ",enodes[j]);
	}
      printf("\n");
      */
      
      _output[om] = enodes;
	
    }
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$EndNodes")
    printf("invalid file format, line != $EndNodes\n");
  getline (input, line);
  printf("present line is %s\n", line.c_str());
  if (line!="$Elements")
    printf("invalid file format, line != $Elements\n");

  
 
  return _output;
  
} 

int main()
{

  D1_cord MESH_node;
  D2_integer CELL2NODE;
  MESH_node = get_mesh_tri(0);
  //CELL2NODE = get_c2n_tri(0);
  CELL2NODE = get_c2n_tri_SAM(0);

  printf("RES_cells = %d, RES_nodes = %d\n",RES_cells, RES_nodes);
  //return 0;
  
  //The supposed mesh information
  printf("---Global node locations, first 100:---\n");
  for (int n = 0; n < fmin(100,RES_nodes); n++)
    {
      printf("node %d: (%f, %f, %f)\n", n, MESH_node[n].rx, MESH_node[n].ry, MESH_node[n].rz);
    }
  printf("\n---Cell2Node pairing, first 100:---\n");
  for (int om = 0; om < fmin(100,RES_cells); om++)
    {
      //printf("element %d: {%d, %d, %d}\n", om, CELL2NODE[om][0], CELL2NODE[om][1], CELL2NODE[om][2]);
      //printf("element %d: {%d, %d, %d, %d, %d, %d, %d, %d}\n", om, CELL2NODE[om][0], CELL2NODE[om][1], CELL2NODE[om][2], CELL2NODE[om][3], CELL2NODE[om][4], CELL2NODE[om][5], CELL2NODE[om][6], CELL2NODE[om][7]);
      printf("element %d: {%d, %d, %d, %d, %d, %d, %d, %d, %d}\n", om, CELL2NODE[om][0], CELL2NODE[om][1], CELL2NODE[om][2], CELL2NODE[om][3], CELL2NODE[om][4], CELL2NODE[om][5], CELL2NODE[om][6], CELL2NODE[om][7], CELL2NODE[om][8]);
    }
  //return 0;
  //Let each cell know the coordinates of its nodes directly.
  //Makes initialization easier
  //D2_cord CELLVERT = Allocate_2D_cord(RES_cells, 8);
  D2_cord CELLVERT = Allocate_2D_cord(RES_cells, Nsolut);
  for (int om = 0; om < RES_cells; om++)
    {
      for (int s = 0; s < Nsolut; s++)
	{
	  CELLVERT[om][s] = MESH_node[CELL2NODE[om][s]];
	}
    }
  
  
  /*
  printf("\n---Cell Vertex Locations:---\n");
  for (int om = 0; om < RES_cells; om++)
    {
      if ((om % 10000) == 0)
	{
	  printf("element %d:\n",om);
	  for (int s = 0; s < N_N; s++)
	    {
	      printf("\tvertex %d: (%f, %f, %f)\n", om, CELLVERT[om][s].rx, CELLVERT[om][s].ry, CELLVERT[om][s].rz);
	    }
	  printf("\n");
	}
    }
  printf("\n---Full Cell node set locations, Lagrangean organization:---\n");
  for (int om = 0; om < RES_cells; om++)
    {
      if ((om % 10000) == 0)
	{
	  printf("element %d:\n",om);
	  for (int s = 0; s < Nsolut; s++)
	    {
	      printf("\tvertex %d: (%f, %f, %f)\n", om, MESH_node[CELL2NODE[om][s]].rx, MESH_node[CELL2NODE[om][s]].ry, MESH_node[CELL2NODE[om][s]].rz);
	    }
	  printf("\n");
	}
    }
  */
  //Populate the primitive degrees of freedom
  //  D3_double DOF = Allocate_3D_double(RES_cells, Nsolut, 5); //DOF primitives (5 field variables)
  D3_double DOF = Allocate_3D_double(RES_cells, Nsolut, Nfield); //DOF primitives (5 field variables)
  FILE*file9001 = fopen(rho_file,"r");
  FILE*file9002 = fopen(ux_file,"r");
  FILE*file9003 = fopen(uy_file,"r");
  //FILE*file9004 = fopen(uz_file,"r");
  FILE*file9005 = fopen(p_file,"r");
  //printf("opened the ux, uy, uz, rho, p files\n");
  printf("opened the ux, uy, rho, p files\n");
  //return 0;
  /*
  int discard_int;
  char discard_str;
  double discard_float;
  fscanf(file9001,"%c,\n", &discard_str);
  printf("discard = %c\n", discard_str);
  fclose(file9001);
  */
  int discard_int;
  char discard_str[50];
  char discard_space;
  double store_float;
  // D2_double b_poly_init = Allocate_2D_double(RES_cells, Nsolut); //primitive degrees of freedom
  for (int L = 0; L < 15; L++) //get past the headers of the files
    {
      fscanf(file9001,"%s,\n", discard_str);
      printf("discard = %s, ", discard_str);
      fscanf(file9002,"%s,\n", discard_str);
      printf("discard = %s, ", discard_str);
      fscanf(file9003,"%s,\n", discard_str);
      printf("discard = %s, ", discard_str);
      //fscanf(file9004,"%s,\n", discard_str);
      //printf("discard = %s ", discard_str);
      fscanf(file9005,"%s,\n", discard_str);
      printf("discard = %s\n ", discard_str);
    }
  // return 0;
  for (int om = 0; om < RES_cells; om++) //Now, read in the primitive DOF
    {
      //Each element om is represented by a row. Throw out the first two entries (cell id, Nsolut)
      for (int col = 0; col < 2; col++)
	{
	  fscanf(file9001, "%d  ", &discard_int);
	  fscanf(file9002, "%d  ", &discard_int);
	  fscanf(file9003, "%d  ", &discard_int);
	  //  fscanf(file9004, "%d  ", &discard_int);
	  fscanf(file9005, "%d  ", &discard_int);
	  //printf("discard = %d,", discard_int);
	}
      for (int col = 2; col < 2+Nsolut; col++) //Now, grab the Nsolut DOF from the element
	{
	  /*
	  fscanf(file9001, "%c ", &discard_space);
	  fscanf(file9001, "%lf ", &store_float);
	  */
	  fscanf(file9001, "%s,", discard_str); //rho
	  store_float = atof(discard_str);
	  DOF[om][col-2][0] = store_float;

	  fscanf(file9002, "%s,", discard_str); //ux
	  store_float = atof(discard_str);
	  DOF[om][col-2][1] = store_float;

	  fscanf(file9003, "%s,", discard_str); //uy
	  store_float = atof(discard_str);
	  DOF[om][col-2][2] = store_float;
	  /*
	  fscanf(file9004, "%s,", discard_str); //uz
	  store_float = atof(discard_str);
	  DOF[om][col-2][3] = store_float;
	  */
	  fscanf(file9005, "%s,", discard_str); //pressure
	  store_float = atof(discard_str);
	  DOF[om][col-2][3] = store_float;
	  //printf("(discard = %s, ", discard_str);
	  //printf("store_float = %f),",store_float);
	}
      fscanf(file9001,"\n");
      fscanf(file9002,"\n");
      fscanf(file9003,"\n");
      //fscanf(file9004,"\n");
      fscanf(file9005,"\n");
      //printf("\n");
    }
  fclose(file9001);
  fclose(file9002);
  fclose(file9003);
  //fclose(file9004);
  fclose(file9005);
  printf("Finished reading in the primitives\n");

  //Report the variables in friendly csv format
  FILE*RHO_OUT = fopen("friendly_rho.csv","w");
  FILE*UX_OUT = fopen("friendly_ux.csv","w");
  FILE*UY_OUT = fopen("friendly_uy.csv","w");
  //FILE*UZ_OUT = fopen("friendly_uz.csv","w");
  FILE*P_OUT = fopen("friendly_p.csv","w");
  for (int om = 0; om < RES_cells; om++)
    {
      for (int k = 0; k < Nsolut; k++)
	{
	  //only the rho file will have geometry information, to save space
	  double xLocal = MESH_node[CELL2NODE[om][k]].rx;
	  double yLocal = MESH_node[CELL2NODE[om][k]].ry;
	  double zLocal = MESH_node[CELL2NODE[om][k]].rz;
	  fprintf(RHO_OUT,"%f,%f,%f,%f\n",xLocal,yLocal,zLocal,DOF[om][k][0]);
	  fprintf(UX_OUT,"%f\n",DOF[om][k][1]);
	  fprintf(UY_OUT,"%f\n",DOF[om][k][2]);
	  //fprintf(UZ_OUT,"%f\n",DOF[om][k][3]);
	  fprintf(P_OUT,"%f\n",DOF[om][k][3]);
	}
    }
  fclose(RHO_OUT);
  fclose(UX_OUT);
  fclose(UY_OUT);
  //fclose(UZ_OUT);
  fclose(P_OUT);

  //Finished reporting primitive variables. Now, calculate vorticity
  // D3_double vorticity = Allocate_3D_double(RES_cells, Nsolut, 3);
  //Use basis functions to calculate the vorticity at each 
  //solution point.
  /*
  for (int om = 0; om < RES_cells; om++)
    {
      //Get the x,y,z coordinates of the nodes
      D1_double XYZNodes = Allocate_2D_double(Nsolut, 3);
      for (int j = 0; j < Nsolut; j++)
	{
	  XYZNodes[j][0] = MESH_node[CELL2NODE[om][j]].rx;
	  XYZNodes[j][1] = MESH_node[CELL2NODE[om][j]].ry;
	  XYZNodes[j][2] = MESH_node[CELL2NODE[om][j]].rz;
	}
      //Nodes known. Use knowledge of dof layout to calculate
      //shape functions
    }
  */


  /*

  //Now, report the vorticity
  FILE*filev;
  filev = fopen("vorticity_post.csv","w");
  for (int om = 0; om < RES_cells; om++)
    {
      for (int i = 0; i < GQout; i++)
	{
	  fprintf(filev,"%f,%f,%f\n",phys_node_GQvol[om][i].rx,phys_node_GQvol[om][i].ry, vorticity[om][i]);
	}
    }
  fclose(filev);
  
  */
 
  return 0;
}


