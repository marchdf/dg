/*!
  \file sensor.h
  \brief Class deals with the sensor for local limiting.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup sensor Sensor
  \ingroup sensor
*/
#ifndef SENSOR_H
#define SENSOR_H
#include <macros.h>
#include <misc.h>
#include <scalar_def.h>
#include <sensor_kernels.h>
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
  bool _isSensor; // true if a sensor is in place
  int* _sensors;
  scalar* _Uavg;
  scalar* _ones;
  
 public:
  
  /*!\brief Constructor */
 SENSOR(int N_s, int N_E, int N_N, bool sensor1=false, bool sensor2=false) : _N_s(N_s), _N_E(N_E), _N_N(N_N), _sensor1(sensor1), _sensor2(sensor2) {

    // Are there any sensors on?
    if (_sensor1 || _sensor2){ _isSensor = true;}
    else                    { _isSensor = false;}

    // Thresholds
    _thresh1 = 0.05;
    _thresh2 = 0.1;
    
    one_div_N_s = 1.0/(scalar)_N_s;
#ifdef USE_CPU
    _sensors = new int[_N_E];
    _Uavg    = new scalar[_N_E*N_F];
    _ones    = new scalar[_N_s];
    for(int e=0; e<_N_E; e++){ _sensors[e] = 0;}
    for(int i=0; i<_N_s; i++){ _ones[i] = 1;}
#elif USE_GPU
    cudaMalloc((void**) &_sensors,_N_E*sizeof(int));
    cudaMalloc((void**) &_Uavg,_N_E*N_F*sizeof(scalar));
    cudaMalloc((void**) &_ones,_N_s*sizeof(scalar));
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

  int* getSensors()const {/*! Return sensor array*/return _sensors;};
  bool isSensor()const {return _isSensor;};
  void sensing(int* neighbors, scalar* U);
  void copy_detected_elements(scalar* Uold, scalar* U);
};
#endif 
