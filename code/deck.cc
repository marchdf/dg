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
  while (line!="#Code options"){
    getline(input,line);
  }
  getline(input,line);
  if (line!="#time integration")
    printf("Invalid file format (at time integration)\n");
  getline(input,_timeMeth);
  getline(input,line);
  if (line!="#output time step size")
    printf("Invalid file format (at output time step size)\n");
  input>>_Dt;
  getline(input,line); getline(input,line);
  if (line!="#final time")
    printf("Invalid file format (at final time)\n");
  input>>_tf;
  getline(input,line); getline(input,line);
  if (line!="#Courant-Friedrichs-Lewy condition")
    printf("Invalid file format (at Courant-Friedrichs-Lewy condition)\n");
  input>>_cfl;
  getline(input,line); getline(input,line);
  if (line!="#order")
    printf("Invalid file format (at order)\n");
  input>>_order;
  getline(input,line); getline(input,line);
  if (line!="#mesh file")
    printf("Invalid file format (at mesh file)\n");
  getline(input,_meshfile);
  getline(input,_elemType);
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
