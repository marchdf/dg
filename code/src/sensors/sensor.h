/*!
  \file sensor.h
  \brief Class deals with the sensor for local limiting.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup sensor
*/
#ifndef SENSOR_H
#define SENSOR_H
#include <vector>
#include "macros.h"
#include "misc.h"
#include "scalar_def.h"
#include "sensor_kernels.h"
#include "timers.h"
#include "mem_counter.h"
#ifdef USE_GPU
#include <cublas.h>
#endif

class SENSOR {


 private:
  int _N_s;
  int _N_E;
  int _N_N;
  scalar one_div_N_s;
  bool _sensor1; scalar _thresh1;
  bool _sensor2; scalar _thresh2;
  bool _sensor3; scalar _thresh3;
  bool _isSensor; // true if a sensor is in place
  int* _sensors;
  scalar* _Uavg;
  scalar* _ones;
  TIMERS &_timers;
  
 public:
  
  /*!\brief Constructor */
 SENSOR(int N_s, int N_E, int N_N, TIMERS &timers, MEM_COUNTER &mem_counter, const std::vector<double> &thresholds = std::vector<double>()) : _N_s(N_s), _N_E(N_E), _N_N(N_N), _timers(timers) {

    _sensor1=false;
    _sensor2=false;
    _sensor3=false;

    // Thresholds from the input deck
    for(int k=0; k<thresholds.size(); k++){
      if      (k==0){_sensor1 = true;_thresh1 = thresholds[k];}
      else if (k==1){_sensor2 = true;_thresh2 = thresholds[k];}
      else if (k==2){_sensor3 = true;_thresh3 = thresholds[k];}
    }

    // Are there any sensors on?
    if (_sensor1 || _sensor2 || _sensor3){ _isSensor = true;}
    else                                 { _isSensor = false;}
    
    one_div_N_s = 1.0/(scalar)_N_s;
#ifdef USE_CPU
    _sensors = new int[_N_E];                                 mem_counter.addToCPUCounter(_N_E*sizeof(int));
    _Uavg    = new scalar[_N_E*N_F];                          mem_counter.addToCPUCounter(_N_E*N_F*sizeof(scalar));
    _ones    = new scalar[_N_s];                              mem_counter.addToCPUCounter(_N_s*sizeof(scalar));
    for(int e=0; e<_N_E; e++){ _sensors[e] = 0;}
    for(int i=0; i<_N_s; i++){ _ones[i] = 1;}
#elif USE_GPU
    cudaMalloc((void**) &_sensors,_N_E*sizeof(int));         mem_counter.addToGPUCounter(_N_E*sizeof(int));	   
    cudaMalloc((void**) &_Uavg,_N_E*N_F*sizeof(scalar));     mem_counter.addToGPUCounter(_N_E*N_F*sizeof(scalar));
    cudaMalloc((void**) &_ones,_N_s*sizeof(scalar));	     mem_counter.addToGPUCounter(_N_s*sizeof(scalar));    
    cudaMemset(_sensors, 0, _N_E*sizeof(int));

    // Roundabout way of setting _ones to 1 (bc memset is a bitch)
    scalar* tmp = new scalar[_N_s]; for(int i=0; i<_N_s; i++){ tmp[i] = 1;}
    cudaMemcpy(_ones, tmp, _N_s*sizeof(scalar), cudaMemcpyHostToDevice);
    delete[] tmp;
#endif
  }

  /*!\brief Destructor */
  ~SENSOR(){
    if(_sensors) del(_sensors);
    if(_Uavg)    del(_Uavg);
    if(_ones)    del(_ones);
  }

  int* getSensors()const {/*!\brief Return sensor array*/return _sensors;};
  bool isSensor()const {/*!\brief Return true if using sensors*/return _isSensor;};
  void sensing(int* neighbors, scalar* U);
  void copy_detected_elements(scalar* Uold, scalar* U);
};
#endif 
