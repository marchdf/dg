/*!
  \file lagrange_particle.cc
  \brief Function definitions for LAGRANGE_PARTICLE class.
  \copyright Copyright (C) 2014, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
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
    _prev_el[k] = Lget_element_belong(&_positions[k*D],_prev_el[k],_m.getNeighbors(),_XYZNodes,_nvert,_N_N,_N_E);
    
    // Get the velocity at that position
    Lget_velocity_at_position(&_positions[k*D],_prev_el[k],_N_s,U,_solution,_avg_velocity);
    
    // Update the position of the particle
    for(int alpha=0; alpha<D; alpha++){
      _positions[k*D+alpha] = Dt*_avg_velocity[alpha] + _positions[k*D+alpha];
      _avg_velocity[alpha] = 0.0;
    }
  }
  
  _timers.stop_timer(28);
}

void LAGRANGE_PARTICLES::printParticles(const double time, scalar* output){
  /*!
    \brief Function to output the particle positions and field values
    \param[in] time output time
    \param[in] output formatted solution from the printer class
    \section Description
    This outputs the time step, the particle position, and the average
    solution in the element containing the particle. A more accurate
    output would interpolate the solution at the particle position
    (instead of averaging) but that's complicated (see Lget_velocity_at_position
    function).
  */
  for(int k = 0; k < _NP; k++){

    // open the file   
    _ofile = fopen(_pnames[k].c_str(),"a");

    // write the time
    fprintf(_ofile,"%20.16E",time);

    // write the particle position
    for(int alpha=0; alpha<D; alpha++){fprintf(_ofile,",%20.16E",_positions[k*D+alpha]);}

    // Average the solution in the element by getting the formatted
    // solution in the printing step
    scalar avg = 0;
    if(_prev_el[k]!=-1){

      // Loop on the fields
      for(int fc=0; fc<N_F; fc++){
	// calculate the average
	for(int i=0; i<_N_s; i++){avg += output[(_prev_el[k]*N_F+fc)*_N_s+i];}

	// output it and set it back to zero
	fprintf(_ofile,",%20.16E",avg/_N_s); avg = 0;
      }
    }
    else{
      for(int fc=0; fc<N_F; fc++){ fprintf(_ofile,",%20.16E",0.0);}
    }

    // end the line and close the file
    fprintf(_ofile,"\n");
    fclose(_ofile);
  } 
}
