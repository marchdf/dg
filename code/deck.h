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
  int _nf;
  std::string _ic;
  std::string _bc;
  
 public:
  inline std::string getTimeMeth() { return _timeMeth;}
  inline double getOutputTimeStep() {return _Dt;}
  inline double getFinalTime() {return _tf;}
  inline double getCFL() {return _cfl;}
  inline int getOrder() {return _order;}
  inline std::string getMeshfile() { return _meshfile;}
  inline std::string getElemType() { return _elemType;}
  inline std::string getLimiter() { return _limiter;}
  inline int getNumFields() {return _nf;}
  inline std::string getInitialCondition() { return _ic;}
  inline std::string getBoundaryCondition() { return _bc;}
  void readDeck(const char *fileName);
};
#endif
