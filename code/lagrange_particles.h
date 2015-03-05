/*!
  \file lagrange_particles.h
  \brief Class deals with the lagrange particles
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \defgroup lagrange_particles Lagrange particles
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
#include <macros.h>
#include <misc.h>
#include <scalar_def.h>
#include <timers.h>
#include <mem_counter.h>
#include "simpleMesh.h"
#include "fullMatrix.h"
#include <lagrange_particles_kernels.h>

class LAGRANGE_PARTICLES {

 private:
  bool _have_particles;
  int _NP;
  int _nvert;
  int _N_N;
  int _N_E;
  int _N_s;
  scalar* _positions;
  int* _prev_el;
  int* _neighbors;
  scalar* _solution;
  scalar* _avg_velocity;
  TIMERS &_timers;
  simpleMesh &_m;
  const fullMatrix<scalar> &_XYZNodes;
  std::vector<std::string> _pnames;
  FILE* _ofile;
    
 public:
  /*!\brief Constructor defaults*/
 LAGRANGE_PARTICLES(TIMERS &timers, MEM_COUNTER &mem_counter, simpleMesh &m, const fullMatrix<scalar> &XYZNodes, int nvert, int N_N, int N_E, int N_s, const std::vector<double> &input_particles = std::vector<double>()) : _timers(timers), _m(m), _XYZNodes(XYZNodes), _nvert(nvert), _N_N(N_N), _N_E(N_E), _N_s(N_s){

    if (input_particles.size() == 0){
      _have_particles = false;
    }
    else{
      _have_particles = true;

      // Get the number of particles we are tracking
      _NP = input_particles[0];
      printf("There are %i Lagrange particles to track. They are located at:\n",_NP);
      
      // Get the initial positions of these particles
      _positions = new scalar[_NP*D]; mem_counter.addToCPUCounter(_NP*D*sizeof(scalar));
      _prev_el = new int[_NP];        mem_counter.addToCPUCounter(_NP*sizeof(int));
      for(int k=0; k<_NP; k++){
	_prev_el[k] = 0;
	printf("\tparticle %i:",k);
	for(int alpha = 0; alpha < D; alpha ++){
	  _positions[k*D+alpha] = input_particles[1+k*D+alpha]; // 1+ because the first element is the number of particles
	  printf(" %8.6f, ",_positions[k*D+alpha]);
	}
	printf("\n");
      }

      // Initialize a vector to store the solution in an element and the velocities
      _solution     = new scalar[(1+D)*_N_s]; // rho, ux, uy
      _avg_velocity = new scalar[D];  for(int alpha=0; alpha<D; alpha++){_avg_velocity[alpha] = 0.0;}

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
  
  /*!\brief Destructor */
  ~LAGRANGE_PARTICLES(){
    if (_have_particles){
      if(_positions) del(_positions);
      if(_prev_el)   del(_prev_el);
      if(_solution)  del(_solution);
      if(_avg_velocity)  del(_avg_velocity);
    }
  }

  bool haveParticles()const {/*!\brief Return true if you have Lagrange particles*/return _have_particles;};
  void advectParticles(scalar Dt, scalar* U);
  void printParticles(scalar time, scalar* output);
  
};
#endif
