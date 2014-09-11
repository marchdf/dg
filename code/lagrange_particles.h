/*!
  \file lagrange_particles.h
  \brief Class deals with the lagrange particles
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
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
  scalar* _positions;
  int* _prev_el;
  int* _neighbors;
  int _nvert;
  int _N_N;
  int _N_E;
  TIMERS &_timers;
  simpleMesh &_m;
  const fullMatrix<scalar> &_XYZNodes;
  
 public:
  /*!\brief Constructor defaults*/
 LAGRANGE_PARTICLES(TIMERS &timers, MEM_COUNTER &mem_counter, simpleMesh &m, const fullMatrix<scalar> &XYZNodes, int nvert, int N_N, int N_E, const std::vector<double> &input_particles = std::vector<double>()) : _timers(timers), _m(m), _XYZNodes(XYZNodes), _nvert(nvert), _N_N(N_N), _N_E(N_E){

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
    }
  }
  
  /*!\brief Destructor */
  ~LAGRANGE_PARTICLES(){
    if(_positions) del(_positions);
    if(_prev_el)   del(_prev_el);
  }

  bool haveParticles()const {/*!\brief Return true if you have Lagrange particles*/return _have_particles;};
  void advectParticles(scalar Dt, scalar* U);
  
};
#endif
