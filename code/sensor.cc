/*!
  \file sensor.cc
  \brief Function definitions for SENSOR class.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup sensor
*/
#include <sensor.h>

void SENSOR::sensing(int* neighbors, scalar* U){
  /*!
    \brief Main sensor function to drive all the sensor calculations
    \param[in] neighbors array of neighbor element indices
    \param[in] U main solution
  */

  _timers.start_timer(21);
  
  // Set the sensors to 0
#ifdef USE_CPU
  for(int e=0; e<_N_E; e++){ _sensors[e] = 0;}
#elif USE_GPU
  cudaMemset(_sensors, 0, _N_E*sizeof(int));
#endif
  
  // Calculate the average solution in each cell. This is NOT the same
  // thing as a cell average.
  blasGemm('N','N', 1, _N_E*N_F, _N_s, one_div_N_s, _ones, 1, U, _N_s, 0.0, _Uavg, 1);
  
  // Call sensor calculation function
  Lcalc_sensors(_N_E,_N_N,_sensor1,_thresh1,_sensor2,_thresh2,neighbors,_Uavg,_sensors);

  _timers.stop_timer(21);
}


void SENSOR::copy_detected_elements(scalar* Uold, scalar* U){
  /*!
    \brief Copy only the elements which were detected with a sensor and limited.
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
  */
  Lcopy_detected(_N_s,_N_E,_sensors,Uold,U);
}
