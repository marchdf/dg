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
#include <scalar_def.h>
#include <sensor_kernels.h>
#ifdef USE_GPU
#include <cuda_runtime_api.h>
#endif

class SENSOR {


 private:
  int _N_E;
  int _N_s;
  bool _sensor1;
  bool _sensor2;
  bool _isSensor; // true if a sensor is in place
  int* _sensors;
  
 public:
  
  /*!\brief Constructor */
 SENSOR(int N_s, int N_E, bool sensor1=false, bool sensor2=false) : _N_s(N_s), _N_E(N_E), _sensor1(sensor1), _sensor2(sensor2) {

    // Are there any sensor on?
    if (_sensor1 || _sensor2){ _isSensor = true;}
    else                    { _isSensor = false;}
    
#ifdef USE_CPU
    _sensors = new int[_N_E]; for(int e=0; e<_N_E; e++){_sensors[e] = 0;}    
#elif USE_GPU
    cudaMalloc((void**) &_sensors,_N_E*sizeof(int));
    cudaMemset(_sensors, 0, _N_E*sizeof(int));
#endif
  }

  /*!\brief Destructor */
  ~SENSOR(){
    if(_sensors) del(_sensors);
  }

  int* getSensors()const {/*! Return sensor array*/return _sensors;}
  bool isSensor()const {return _isSensor;}
  void sensing(scalar* U);
  void copy_detected_elements(scalar* Uold, scalar* U);
};
#endif 
