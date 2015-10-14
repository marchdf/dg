/*!
  \file lagrange_particles.h
  \brief Class deals with the lagrange particles
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup lagrange_particles

  This is a pure CPU implementation of Lagrange particles. The reasons
  for this are detailed in my notes 10/09/14 but it boils down to the
  following:
  
  - the GPU would need to know XYZ coordinates of the elements (lots
    of memory use!!) to locate the particle in an element

  - it would be necessary to implement a reduction/search in parallel
    on the GPU to (1) find which elements the particles belong to and
    (2) do this efficiently. This can be complicated to implement well
    on a GPU.

  - The main downside: at every time step, to advect the particle on
    the CPU, I will need to communicate the field quantities from the
    GPU to the CPU. This array will be small (N_F * N_s) but
    still... I think I can live with this because I will most likely
    not be running many GPU simulations which require Lagrange
    particles.  
*/
#ifndef LAGRANGE_PARTICLES_H
#define LAGRANGE_PARTICLES_H
#include <vector>
#include "macros.h"
#include "misc.h"
#include "scalar_def.h"
#include "timers.h"
#include "mem_counter.h"
#include "simpleMesh.h"
#include "fullMatrix.h"
#include "lagrange_particles_kernels.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

class LAGRANGE_PARTICLES {

 private:
  bool _have_particles;
  int _NP;
  int _nvert;
  int _N_N;
  int _N_E;
  int _N_s;
  int _myid;
  int _numprocs;
  scalar* _positions;
  int* _prev_el;
  int* _prev_part;
  int* _neighbors;
  scalar* _solution;
  scalar* _avg_velocity;
  scalar* _local_output;
  int* _finding_proc_array;
  TIMERS &_timers;
  simpleMesh &_m;
  const fullMatrix<scalar> &_XYZNodes;
  std::vector<std::string> _pnames;
  FILE* _ofile;

#ifdef USE_MPI
  int _sendtag;
  int _recvtag;
  MPI_Status *_status;
  MPI_Request *_request;
#endif

  
 public:
  /*!\brief Constructor defaults*/
 LAGRANGE_PARTICLES(TIMERS &timers, MEM_COUNTER &mem_counter, simpleMesh &m, const fullMatrix<scalar> &XYZNodes, int nvert, int N_N, int N_E, int N_s, int myid, int numprocs, const std::vector<double> &input_particles = std::vector<double>()) : _timers(timers), _m(m), _XYZNodes(XYZNodes), _nvert(nvert), _N_N(N_N), _N_E(N_E), _N_s(N_s), _myid(myid), _numprocs(numprocs){

    // Initialization
    _have_particles = false; // no particles by default
    _sendtag = 0;
    _recvtag = 0;

    
    if (input_particles.size() != 0){
      _have_particles = true;
      
      // Get the number of particles we are tracking
      _NP = input_particles[0];
      
      // Get the initial positions of these particles
      _positions = new scalar[_NP*D]; mem_counter.addToCPUCounter(_NP*D*sizeof(scalar));
      _prev_el = new int[_NP];        mem_counter.addToCPUCounter(_NP*sizeof(int));
      _prev_part = new int[_NP];      mem_counter.addToCPUCounter(_NP*sizeof(int));
      for(int k=0; k<_NP; k++){
	_prev_el[k] = 0;
	_prev_part[k] = 0;
	for(int alpha = 0; alpha < D; alpha ++){
	  _positions[k*D+alpha] = input_particles[1+k*D+alpha]; // 1+ because the first element is the number of particles
	}
      }

      // Initialize a vector to store the solution in an element and the velocities
      _solution     = new scalar[(1+D)*_N_s]; // rho, ux, uy
      _avg_velocity = new scalar[D];  for(int alpha=0; alpha<D; alpha++){_avg_velocity[alpha] = 0.0;}

      // init some other pointers to null (to avoid weirdness)
      _local_output = NULL;
      _finding_proc_array  = NULL;
#ifdef USE_MPI
      _status = new MPI_Status [_NP];
      _request = new MPI_Request [_NP];
#endif
     
      // Only the zeroth processor prints stuff
      if (_myid == 0){

	// Output particle location (just to see them)
	printf("There are %i Lagrange particles to track. They are located at:\n",_NP);
	for(int k=0; k<_NP; k++){ 
	  printf("\tparticle %i:",k);
	  for(int alpha = 0; alpha < D; alpha ++){
	    printf(" %8.6f, ",_positions[k*D+alpha]);
	  }
	  printf("\n");
	}
	
	// initialize a local output vector to store things from other processors
	_local_output = new scalar[_NP*N_F*_N_s];  mem_counter.addToCPUCounter(_NP*N_F*_N_s*sizeof(scalar));
	_finding_proc_array  = new int[_numprocs]; mem_counter.addToCPUCounter(_numprocs*sizeof(int)); 
	
	// Initialize the particle files names
	char buffer [256];
	std::string name;
	for(int k=0; k<_NP; k++){
	  sprintf(buffer, "particle%03i.dat", k);
	  name = buffer;
	  _pnames.push_back(name);
	}
      
	// Output a comment line for each particle
	for(int k = 0; k < _NP; k++){
	  // open the file   
	  _ofile = fopen(_pnames[k].c_str(),"w");
	  // write the comment line
	  fprintf(_ofile,"# time, position, and average solution for particle %i.\n",k);
	  // end the line and close the file
	  fclose(_ofile);
	}
      }
    }
  }
  
  /*!\brief Destructor */
  ~LAGRANGE_PARTICLES(){
    if (_have_particles){
      if(_positions)             del(_positions);
      if(_prev_el)               del(_prev_el);
      if(_prev_part)             del(_prev_part);
      if(_solution)              del(_solution);
      if(_avg_velocity)          del(_avg_velocity);
      if(_local_output)          del(_local_output);
      if(_finding_proc_array)    del(_finding_proc_array);
    }
  }

  bool haveParticles()const {/*!\brief Return true if you have Lagrange particles*/return _have_particles;};
  void advectParticles(scalar Dt, scalar* U);
  void printParticles(scalar time, scalar* output);
  
};
#endif
