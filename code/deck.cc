/*!
  \file deck.cc
  \brief Functions for the deck class
  \copyright Copyright (C) 2014, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "deck.h"


void deck::readDeck(const char *fileName)
{
  /*!
    \brief Read the input deck
    \param[in] fileName deck filename
  */    
  std::ifstream input;
  input.open(fileName,std::ifstream::in);
  if(input.is_open()==0){
    printf("No file named %s. Defaulting to deck.inp\n",fileName);
    input.open("deck.inp",std::ifstream::in);
    if(input.is_open()==0){
      printf("No file named deck.inp. Exiting.\n");
      exit(1);
    }
  }

  // Ignore the lines before the code options
  std::string line;
  getline(input,line);
  while (line!="#Code options"){
    getline(input,line);
  }

  // Parse the deck
  int cnt = 0;
  while(getline(input, line)){
    //std::cout<< "Parsing line: "<<line<<std::endl;
    if      (line=="#time integration")                 {getline(input,_timeMeth); cnt++;}
    else if (line=="#output time step size")            {input>>_Dt; getline(input,line); cnt++;}
    else if (line=="#final time")                       {input>>_tf;getline(input,line); cnt++;}
    else if (line=="#Courant-Friedrichs-Lewy condition"){input>>_cfl; getline(input,line); cnt++;}
    else if (line=="#order")                            {input>>_order; getline(input,line); cnt++;}
    else if (line=="#mesh file")                        {getline(input,_meshfile); cnt++; getline(input,_elemType); cnt++;}
    else if (line=="#limiter")                          {getline(input,_limiter); cnt++;}
    else if (line=="#initial condition")                {input>>_ic; cnt++; double t; while(input>> t){_ic_inputs.push_back(t);} input.clear();}
    else if (line=="#sensor thresholds")                {double t; while(input>> t){_thresholds.push_back(t);} input.clear();}
    else if (line=="#lagrange particles")               {double t; while(input>> t){_lagrange_particles.push_back(t);} input.clear();}
    else if (line=="#restart step")                     {input>>_restart_step; getline(input,line);}
    else{
      std::cout<<"Unrecognized option in deck: "<< line << std::endl;
      std::cout<<"Ignoring for now."<< std::endl;
    } // if on line
  }

  // Check to make sure we read all the mandatory options (9 right now)
  if(cnt < 9){
    std::cout<<"Input deck is incomplete. Exiting."<< std::endl;
    exit(1);
  }

  input.close();
}
