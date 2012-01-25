#ifndef DECK_H
#define DECK_H
#include <string>

class deck {

 private:
  std::string _proc;
  bool _debug;
  bool _blas;
  std::string _timeMeth;
  double _Dt;
  int _nt;
  int _ofactor;
  std::string _spatialMeth;
  int _order;
  std::string _flux;
  std::string _meshfile;
  std::string _elemType;
  std::string _problem;
  std::string _model;
  int _nf;
  std::string _ic;
  std::string _bc;
  
 public:
  inline std::string getProc() { return _proc;}
  inline bool getDebug() {return _debug;}
  inline bool getBlas() {return _blas;}
  inline std::string getTimeMeth() { return _timeMeth;}
  inline double getTimeStep() {return _Dt;}
  inline int getNumberStep() {return _nt;}
  inline int getOutputFactor() {return _ofactor;}
  inline std::string getSpatialMeth() { return _spatialMeth;}
  inline int getOrder() {return _order;}
  inline std::string getFlux() { return _flux;}
  inline std::string getMeshfile() { return _meshfile;}
  inline std::string getElemType() { return _elemType;}
  inline std::string getProblem() { return _problem;}
  inline std::string getModel() { return _model;}
  inline int getNumFields() {return _nf;}
  inline std::string getInitialCondition() { return _ic;}
  inline std::string getBoundaryCondition() { return _bc;}
  void readDeck(const char *fileName);
};
#endif
