#include "deck.h"
#include <fstream>
#include <sstream>
#include <iostream>


void deck::readDeck(const char *fileName)
{
  std::ifstream input;
  input.open(fileName,std::ifstream::in);
  if(input.is_open()==0){
    printf("No file named %s. Defaulting to deck.inp\n",fileName);
    input.open("deck.inp",std::ifstream::in);
  }
  std::string line;
  getline(input,line);
  if (line!="#Dimensions")
    printf("Invalid file format (at dimensions)\n");
  input>>_dims;
  getline(input,line); getline(input,line);
  if (line!="#Options")
    printf("Invalid file format (at options)\n");
  getline(input,line);
  if (line=="blas")
    _blas = true;
  else _blas = false;
  getline(input,line);
  if (line!="#time integration")
    printf("Invalid file format (at time integration)\n");
  getline(input,_timeMeth);
  getline(input,line);
  if (line!="#time step size")
    printf("Invalid file format (at time step size)\n");
  input>>_Dt;
  getline(input,line); getline(input,line);
  if (line!="#number of time steps")
    printf("Invalid file format (at number of time steps)\n");
  input>>_nt;
  getline(input,line); getline(input,line);
  if (line!="#output factor")
    printf("Invalid file format (at output factor)\n");
  input>>_ofactor;
  getline(input,line); getline(input,line);
  if (line!="#spatial method")
    printf("Invalid file format (at spatial method)\n");
  getline(input,_spatialMeth);
  getline(input,line);
  if (line!="#order")
    printf("Invalid file format (at order)\n");
  input>>_order;
  getline(input,line); getline(input,line);
  if (line!="#flux")
    printf("Invalid file format (at flux)\n");
  getline(input,_flux);
  getline(input,line);
  if (line!="#mesh file")
    printf("Invalid file format (at mesh file)\n");
  getline(input,_meshfile);
  getline(input,_elemType);
  getline(input,line);
  if (line!="#problem")
    printf("Invalid file format (at problem)\n");
  getline(input,_problem);
  getline(input,line);
  if (line!="#model")
    printf("Invalid file format (at model)\n");
  getline(input,_model);
  getline(input,line);
  if (line!="#limiter")
    printf("Invalid file format (at limiter)\n");
  getline(input,_limiter);
  getline(input,line);
  if (line!="#number of fields")
    printf("Invalid file format (at number of fields)\n");
  input>>_nf;  
  getline(input,line); getline(input,line);
  if (line!="#initial condition")
    printf("Invalid file format (at initial condition)\n");
  getline(input,_ic);
  getline(input,line);
  if (line!="#boundary condition")
    printf("Invalid file format (at boundary condition)\n");
  getline(input,_bc);

  input.close();
}
