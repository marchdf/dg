/*!
  \file sensor.cc
  \brief Function definitions for SENSOR class.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup sensor
*/
#include <sensor.h>

void SENSOR::sensing(scalar* U){
  /*!
    \brief Main sensor function to drive all the sensor calculations
    \param[in] U main solution
  */

  // Use the first sensor
  if (_sensor1){
    Lsensor1(_N_s,_N_E,U,_sensors);
  }

  // Use the second sensor
  if (_sensor2){
    
  }

}


void SENSOR::copy_detected_elements(scalar* Uold, scalar* U){
  /*!
    \brief Copy only the elements which were detected with a sensor and limited.
    \param[in] Uold solution containing the limited solutions
    \param[out] U solution to receive the limited elements
  */
  Lcopy_detected(_N_s,_N_E,_sensors,Uold,U);
}
