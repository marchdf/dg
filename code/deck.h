/*!
  \file deck.h
  \brief Class deals with the input deck.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
*/
#ifndef DECK_H
#define DECK_H
#include <string>

class deck {

 private:
  std::string _timeMeth;
  double _Dt;
  double _tf;
  double _cfl;
  int _order;
  std::string _flux;
  std::string _meshfile;
  std::string _elemType;
  std::string _limiter;
  std::string _ic;
  std::string _bc;
  
 public:
  inline std::string getTimeMeth() {/*!Return time integration method*/ return _timeMeth;}
  inline double getOutputTimeStep() {/*!Return output time step*/return _Dt;}
  inline double getFinalTime() {/*!Return final time*/return _tf;}
  inline double getCFL() {/*!Return CFL number*/return _cfl;}
  inline int getOrder() {/*!Return DG order*/return _order;}
  inline std::string getMeshfile() {/*!Return the name of the mesh file*/return _meshfile;}
  inline std::string getElemType() {/*!Return type of element in the mesh*/return _elemType;}
  inline std::string getLimiter() {/*!Return limiting method*/return _limiter;}
  inline std::string getInitialCondition() {/*!Return initial condition*/return _ic;}
  inline std::string getBoundaryCondition() {/*!Return BC (not used anymore)*/return _bc;}
  void readDeck(const char *fileName);
};
#endif
