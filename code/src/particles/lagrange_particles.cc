/*!
  \file lagrange_particle.cc
  \brief Function definitions for LAGRANGE_PARTICLE class.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
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

#ifdef USE_MPI
 
  // Share the particle information to all procs
  MPI_Bcast(_prev_part,_NP,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(_prev_el,  _NP,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(_positions,_NP*D,MPI_SCALAR,0,MPI_COMM_WORLD);
  
  // Loop on all the particles
  for(int k=0; k<_NP; k++){

    // if the particle left the domain, skip all this
    if (_prev_el[k] != -1){

      // reset the found flags on all processors
      if (_myid==0){for(int k=0; k<_numprocs; k++){_finding_proc_array[k] = -1;}}
      int finding_proc = -1;

      //
      // First, find where the particle is
      //
    
      // first check the previous partition (most likely the particle is still there)
      if (_myid == _prev_part[k]){
      	_prev_el[k] = Lget_element_belong(&_positions[k*D],_prev_el[k],_m.getNeighbors(),_XYZNodes,_nvert,_N_N,_N_E);
      	if (_prev_el[k] != -1){finding_proc = _myid;}
      }
      
      // let everyone else know if the particle was found
      MPI_Bcast(&finding_proc, 1, MPI_INT, _prev_part[k], MPI_COMM_WORLD);
      
      // if wasn't found by the previous partition, check all the
      // other partitions (starting at default element 0)
      if ((finding_proc==-1) && (_myid != _prev_part[k])){
      	_prev_el[k] = Lget_element_belong(&_positions[k*D],0,_m.getNeighbors(),_XYZNodes,_nvert,_N_N,_N_E);
      	if (_prev_el[k] != -1){finding_proc = _myid;}
      }
	    
      // Gather all the processors' finding flag so that proc 0 can
      // figure out who found the particle
      MPI_Gather(&finding_proc, 1, MPI_INT, _finding_proc_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(_myid==0){

	// get which processor found the particle. Assume only one proc found the particle.
      	for(int i=0; i<_numprocs; i++){ if (_finding_proc_array[i] != -1) finding_proc = _finding_proc_array[i];}

	// assign the particle to that partition
	_prev_part[k] = finding_proc;

      }

      // If you didn't find the particle, just move on to the next
      // particle. If none of the processors found the particle, this
      // condition makes sure that the element is set to -1 and
      // that the particle will no longer be accounted for.
      if (finding_proc == -1){_prev_el[k]=-1; continue;}

      // if you are the processor who found the particle, get the velocity field
      if (_myid == finding_proc){
	Lget_velocity_at_position(&_positions[k*D],_prev_el[k],_N_s,U,_solution,_avg_velocity);
      }
      
      // If the finding proc is processor 0, don't do any MPI
      // communications because the data is already there. Otherwise,
      // you need to send data to proc 0
      if (finding_proc !=0){
	if (_myid == finding_proc){
      	  MPI_Isend(&_prev_el[k], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &_request[0]);
      	  MPI_Isend(_avg_velocity, D, MPI_SCALAR, 0, 1, MPI_COMM_WORLD, &_request[1]);
      	}
      	else if (_myid == 0){
      	  MPI_Irecv(&_prev_el[k], 1, MPI_INT, finding_proc, 0, MPI_COMM_WORLD, &_request[0]);
      	  MPI_Irecv(_avg_velocity, D, MPI_SCALAR, finding_proc, 1, MPI_COMM_WORLD, &_request[1]);
      	}
      	// Wait for communication to end
      	MPI_Waitall(2, _request, _status);
      }
      
      //
      // Next, advect the particle. Only proc 0 does this
      //
      if (_myid == 0){
	for(int alpha=0; alpha<D; alpha++){
	  _positions[k*D+alpha] = Dt*_avg_velocity[alpha] + _positions[k*D+alpha];
	  _avg_velocity[alpha] = 0.0;
	}
      }
    }
  }
  
#else
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
#endif
  
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

  // First, get output solution from other processors (if necessary) on to processor 0.
  // If no MPI, just copy the solution locally

#ifdef USE_MPI

  // make sure everyone is on the same page and has the updated partition information
  MPI_Bcast(_prev_part,_NP,MPI_INT,0,MPI_COMM_WORLD);
  int cnt = 0;

  // Loop on particles
  for(int k = 0; k < _NP; k++){

    // if the particle left the domain, skip all this
    if (_prev_el[k] != -1){
    
      // If the data is already on proc 0 and I am proc 0, no need to communicate
      if ((_myid == _prev_part[k]) && (_myid==0)){
	for(int fc=0; fc<N_F; fc++){
	  for(int i=0; i<_N_s; i++){
	    _local_output[(k*N_F+fc)*_N_s+i] = output[(_prev_el[k]*N_F+fc)*_N_s+i];
	  }
	}
      }
    
      // Otherwise, you have to get the data to proc 0
      else{
	if (_myid == _prev_part[k]){
	  MPI_Isend(&output[_prev_el[k]*N_F*_N_s], N_F*_N_s, MPI_SCALAR, 0, k, MPI_COMM_WORLD, &_request[cnt]); cnt++;
	}
	else if (_myid == 0){
	  MPI_Irecv(&_local_output[k*N_F*_N_s], N_F*_N_s, MPI_SCALAR, _prev_part[k], k, MPI_COMM_WORLD, &_request[cnt]); cnt++;
	}
      }
    }
  }
  // Wait for communication to end
  MPI_Waitall(cnt, _request, _status);

    
#else
  for(int k = 0; k < _NP; k++){
    for(int fc=0; fc<N_F; fc++){
      for(int i=0; i<_N_s; i++){
	_local_output[(k*N_F+fc)*_N_s+i] = output[(_prev_el[k]*N_F+fc)*_N_s+i];	
      }
    }
  }
#endif

  // Now, output the solution (only for proc 0)
  if (_myid == 0){
    for(int k = 0; k < _NP; k++){

      // Only do all this if the particle is still in the domain
      if(_prev_el[k]!=-1){
	
	// open the file   
	_ofile = fopen(_pnames[k].c_str(),"a");

	// write the time
	fprintf(_ofile,"%20.16E",time);

	// write the particle position
	for(int alpha=0; alpha<D; alpha++){fprintf(_ofile,",%20.16E",_positions[k*D+alpha]);}
	
	// Average the solution in the element by getting the formatted
	// solution in the printing step
	scalar avg = 0;
	// Loop on the fields
	for(int fc=0; fc<N_F; fc++){
	  // calculate the average
	  for(int i=0; i<_N_s; i++){
	    avg += _local_output[(k*N_F+fc)*_N_s+i];
	  }	

	  // output it and set it back to zero
	  fprintf(_ofile,",%20.16E",avg/_N_s); avg = 0;
	}

	// end the line and close the file
	fprintf(_ofile,"\n");
	fclose(_ofile);
      }
    }
  }
}
