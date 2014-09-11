/*!
  \file lagrange_particle.cc
  \brief Function definitions for LAGRANGE_PARTICLE class.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup lagrange_particles
*/
#include "lagrange_particles.h"

void LAGRANGE_PARTICLES::advectParticles(scalar Dt, scalar* U){
  /*!
    \brief Function which advects a set of particles
    \param[in] Dt time step over which to advect the particle
    \param[in] U main solution
  */
  _timers.start_timer(28);

  // Loop on the particles
  for(int k=0; k<_NP; k++){

    // Get the element it's in
    int e = Lget_element_belong(&_positions[k],_prev_el[k],_m.getNeighbors(),_XYZNodes,_nvert,_N_N,_N_E);
    
    // Get the velocity at that position
    
    
    // Update the position of the particle
    
    // Update the element it belongs to
    _prev_el[k] = e;
  }
  
  _timers.stop_timer(28);
}
