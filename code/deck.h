/*!
  \file deck.h
  \brief Class deals with the input deck.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef DECK_H
#define DECK_H
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>

class deck {

 private:
  std::string _timeMeth;
  double _Dt;
  double _tf;
  double _cfl;
  int _order;
  std::string _meshfile;
  std::string _elemType;
  std::string _limiter;
  std::string _ic;
  std::vector<double> _ic_inputs;  
  std::vector<double> _thresholds;
  std::vector<double> _lagrange_particles;
  int _restart_step;
  
 public:
  /*!\brief Constructor defaults*/ 
  deck(){_restart_step = 0;};
  inline std::string getTimeMeth() {/*!Return time integration method*/ return _timeMeth;}
  inline double getOutputTimeStep() {/*!Return output time step*/return _Dt;}
  inline double getFinalTime() {/*!Return final time*/return _tf;}
  inline double getCFL() {/*!Return CFL number*/return _cfl;}
  inline int getOrder() {/*!Return DG order*/return _order;}
  inline std::string getMeshfile() {/*!Return the name of the mesh file*/return _meshfile;}
  inline std::string getElemType() {/*!Return type of element in the mesh*/return _elemType;}
  inline std::string getLimiter() {/*!Return limiting method*/return _limiter;}
  inline std::string getInitialCondition() {/*!Return initial condition*/return _ic;}
  inline const std::vector<double> & getInitialConditionInputs() const {/*!Return initial condition inputs*/return _ic_inputs;}
  inline const std::vector<double> & getThresholds() const {/*!Return thresholds*/return _thresholds;}
  inline const std::vector<double> & getLagrangeParticles() const {/*!Return Lagrange particles*/return _lagrange_particles;}
  inline int getRestartStep()  {/*!Return restart step*/ return _restart_step;}
  void readDeck(const char *fileName);
};
#endif
