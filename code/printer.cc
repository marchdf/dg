/*!
  \file printer.cc
  \brief Function definitions for the PRINTER class.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup printer
*/
#include "printer.h"

void PRINTER::set_names(){
  _sname.push_back("sensor"); // sensor name
#ifdef ONED
#ifdef PASSIVE
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux"); 
  _names.push_back("P");     _fnames.push_back("p");
  _names.push_back("PhiC");  _fnames.push_back("phic");
  _names.push_back("PhiNC"); _fnames.push_back("phinc");
#elif MULTIFLUID
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux"); 
  _names.push_back("G");     _fnames.push_back("g");  
  _names.push_back("P");     _fnames.push_back("p");  
#elif STIFFENED
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux"); 
  _names.push_back("G");     _fnames.push_back("g");
  _names.push_back("Pinf");  _fnames.push_back("pinf");  
  _names.push_back("P");     _fnames.push_back("p");  
#endif // problem type
#elif TWOD
#ifdef PASSIVE
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux");
  _names.push_back("Uy");    _fnames.push_back("uy"); 
  _names.push_back("P");     _fnames.push_back("p");
  _names.push_back("PhiC");  _fnames.push_back("phic");
  _names.push_back("PhiNC"); _fnames.push_back("phinc");
#elif MULTIFLUID
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux"); 
  _names.push_back("Uy");    _fnames.push_back("uy"); 
  _names.push_back("G");     _fnames.push_back("g");  
  _names.push_back("P");     _fnames.push_back("p");  
#elif STIFFENED
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux");
  _names.push_back("Uy");    _fnames.push_back("uy"); 
  _names.push_back("G");     _fnames.push_back("g");
  _names.push_back("Pinf");  _fnames.push_back("pinf");  
  _names.push_back("P");     _fnames.push_back("p");  
#endif // problem type
#endif // dimensions

  // Mass fraction names
  char buffer1 [5]; char buffer2 [2];
  std::string fname;
  std::string name;
#include "loopstart.h"
#define LOOP_END N_Y
#define MACRO(x) sprintf(buffer2, "Y%i", x); name = buffer2; sprintf(buffer1, "y%i", x); fname = buffer1;  _names.push_back(name); _fnames.push_back(fname);
#include "loop.h"

} // end set names


void PRINTER::print(scalar* U, const int step, const double time){
  /*!
    \brief Output solution for the PRINTER class
  */

  // format output variables
#ifdef USE_CPU
  Lformater(_N_s,_N_E,U,_output);
#elif USE_GPU
  Lformater(_N_s,_N_E,U,_d_output);
  cudaMemcpy(_output, _d_output, _N_s*_N_E*N_F*sizeof(scalar), cudaMemcpyDeviceToHost);
#endif

  
  // print to the output file
  _m.writeSolution(_output, _N_s, _N_E, _elem_type, _fnames, _names, step, time);
  
} // end print_dg


void PRINTER::print_sensor(SENSOR &sensor, const int step, const double time){
  /*!
    \brief Output sensor if needed
  */

  // Check to make sure a sensor is on
  if(sensor.isSensor()){
  
  // format output variables
#ifdef USE_CPU
    Lformat_sensor(_N_s,_N_E,sensor.getSensors(),_output);
#elif USE_GPU
    Lformat_sensor(_N_s,_N_E,sensor.getSensors(),_d_output);
    cudaMemcpy(_output, _d_output, _N_s*_N_E*sizeof(scalar), cudaMemcpyDeviceToHost);
#endif

  
    // print to the output file
    _m.writeSolution(_output, _N_s, _N_E, _elem_type, _sname, _sname, step, time);
  }
} 

