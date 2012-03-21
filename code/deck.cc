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
  if (line!="#Options")
    printf("Invalid file format\n");
  getline(input,_proc);
  getline(input,line);
  if (line=="debug")
    _debug = true;
  else _debug = false;
  getline(input,line);
  if (line=="blas")
    _blas = true;
  else _blas = false;
  getline(input,line);
  if (line!="#time integration")
    printf("Invalid file format\n");
  getline(input,_timeMeth);
  getline(input,line);
  if (line!="#time step size")
    printf("Invalid file format\n");
  input>>_Dt;
  getline(input,line); getline(input,line);
  if (line!="#number of time steps")
    printf("Invalid file format\n");
  input>>_nt;
  getline(input,line); getline(input,line);
  if (line!="#output factor")
    printf("Invalid file format\n");
  input>>_ofactor;
  getline(input,line); getline(input,line);
  if (line!="#spatial method")
    printf("Invalid file format\n");
  getline(input,_spatialMeth);
  getline(input,line);
  if (line!="#order")
    printf("Invalid file format\n");
  input>>_order;
  getline(input,line); getline(input,line);
  if (line!="#flux")
    printf("Invalid file format\n");
  getline(input,_flux);
  getline(input,line);
  if (line!="#mesh file")
    printf("Invalid file format\n");
  getline(input,_meshfile);
  getline(input,_elemType);
  getline(input,line);
  if (line!="#problem")
    printf("Invalid file format\n");
  getline(input,_problem);
  getline(input,line);
  if (line!="#model")
    printf("Invalid file format\n");
  getline(input,_model);
  getline(input,line);
  if (line!="#limiter")
    printf("Invalid file format\n");
  getline(input,_limiter);
  getline(input,line);
  if (line!="#number of fields")
    printf("Invalid file format\n");
  input>>_nf;  
  getline(input,line); getline(input,line);
  if (line!="#initial condition")
    printf("Invalid file format\n");
  getline(input,_ic);
  getline(input,line);
  if (line!="#boundary condition")
    printf("Invalid file format\n");
  getline(input,_bc);

  input.close();
}
