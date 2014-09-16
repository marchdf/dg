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
#elif SINGLEFLUID
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux"); 
  _names.push_back("P");     _fnames.push_back("p");
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
#elif SINGLEFLUID
  _names.push_back("Rho");   _fnames.push_back("rho");
  _names.push_back("Ux");    _fnames.push_back("ux");
  _names.push_back("Uy");    _fnames.push_back("uy"); 
  _names.push_back("P");     _fnames.push_back("p");
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


void PRINTER::print(scalar* U, LAGRANGE_PARTICLES &particles, const int step, const double time){
  /*!
    \brief Output solution for the PRINTER class
    \param[in] U solution to output
    \param[in] particles lagrange particles to output
    \param[in] step time step number
    \param[in] time time value
  */

  // format output variables
  _timers.start_timer(23);
#ifdef USE_CPU
  Lformater(_N_s,_N_E,U,_output);
#elif USE_GPU
  Lformater(_N_s,_N_E,U,_d_output);
  cudaMemcpy(_output, _d_output, _N_s*_N_E*N_F*sizeof(scalar), cudaMemcpyDeviceToHost);
#endif
  _timers.stop_timer(23);
    
  // print to the output file
  _timers.start_timer(24);
  _m.writeSolution(_output, _N_s, _N_E, _elem_type, _fnames, _names, step, time);
  _timers.stop_timer(24);

  // Output the particles if you have any
  _timers.start_timer(27);
  if(particles.haveParticles()){particles.printParticles(time,_output);}
  _timers.stop_timer(27);
  
} // end print

void PRINTER::read(const int step, double &time, scalar* U){
  /*!
    \brief Read solution from output files for the PRINTER class
    \param[in] step time step number
    \param[out] time time of the output solution
    \param[out] U solution to populate from files
  */

  // Read the output files
  _m.readSolution(_N_s, _N_E, _elem_type, _fnames, _names, step, time, _output);

  // format the output to the main solution
#ifdef USE_CPU
  Lformater(_N_s,_N_E,U,_output,true);
#elif USE_GPU
  cudaMemcpy(_d_output, _output, _N_s*_N_E*N_F*sizeof(scalar), cudaMemcpyHostToDevice);
  Lformater(_N_s,_N_E,U,_d_output,true);
#endif
  
} // end read


void PRINTER::print_sensor(SENSOR &sensor, const int step, const double time){
  /*!
    \brief Output sensor if needed
    \param[in] sensor sensor to output
    \param[in] step time step number
    \param[in] time time value
  */

  // Check to make sure a sensor is on
  if(sensor.isSensor()){
  
    // format output variables
    _timers.start_timer(25);
#ifdef USE_CPU
    Lformat_sensor(_N_s,_N_E,sensor.getSensors(),_output);
#elif USE_GPU
    Lformat_sensor(_N_s,_N_E,sensor.getSensors(),_d_output);
    cudaMemcpy(_output, _d_output, _N_s*_N_E*sizeof(scalar), cudaMemcpyDeviceToHost);
#endif
    _timers.stop_timer(25);
    
    // print to the output file
    _timers.start_timer(26);
    _m.writeSolution(_output, _N_s, _N_E, _elem_type, _sname, _sname, step, time);
    _timers.stop_timer(26);
  }
} 

