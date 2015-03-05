/*!
  \file lagrange_particles_kernels.cu
  \brief Kernels used by the LAGRANGE_PARTICLES class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup lagrange_particles
*/
#include <lagrange_particles_kernels.h>
#include <stdio.h>

//==========================================================================
//
// Internal prototype function definitions
//
//==========================================================================
void get_verts_of_element(const int e, const int nvert, const fullMatrix<scalar> XYZNodes, scalar* verts);
int pnpoly(int nvert, scalar *verts, scalar* position);

//==========================================================================
//
// Kernel definitions
//
//==========================================================================


//==========================================================================
//
//  Host C functions
//
//==========================================================================
int Lget_element_belong(scalar* position, int prev_el, const int* neighbors, const fullMatrix<scalar> XYZNodes, const int nvert, const int N_N, const int N_E){
  /*!
    \brief Host C function to launch figure out which element the point is in
    \param[in] position coordinates of particle
    \param[in] prev_el element the point belonged to previously
    \param[in] neighbors index of neighboring elements
    \param[in] XYZNodes matrix of the nodal coordinates
    \param[in] nvert number of vertices per element
    \param[in] N_N number of neighbors
    \param[in] N_E number of elements
    \return elnum index to the element the point belongs to
  */

  //
  // If the element the point previously belonged to is -1 (point left
  // the domain previously), don't do anything
  //
  if (prev_el == -1){return -1;}

  //
  // Otherwise find the element the point belongs to. We do these
  // tests in order of most probable location
  //
  
  // Test if the point is in the current element
  scalar* verts = new scalar[nvert*D];
  get_verts_of_element(prev_el,nvert,XYZNodes,verts);
  if (pnpoly(nvert, verts, position)){return prev_el;}

  // Test if the point is in the neighbor elements
  int neighbor = 0;
  for(int nn=0; nn<N_N; nn++){
    neighbor = neighbors[prev_el*N_N+nn]; // neighbor index
    if (neighbor >= 0){
      get_verts_of_element(neighbor,nvert,XYZNodes,verts);
      if (pnpoly(nvert, verts, position)){ return neighbor;}
    }
  }

  // Test all elements
  for(int e=0; e<N_E; e++){
    get_verts_of_element(e,nvert,XYZNodes,verts);
    if (pnpoly(nvert, verts, position)){ return e;}
  } 
  
  // Default to -1 (not found in the domain)
  return -1;
};

void Lget_velocity_at_position(scalar* position, int el, int N_s, scalar* U, scalar* solution, scalar* avg_velocity){
  /*!
    \brief Host C function to get the velocity at a given position
    \param[in] position coordinates of particle
    \param[in] el element index the particle belongs to
    \param[in] N_s number of nodes
    \param[in] U solution in all the elements
    \param[in] solution the nodal solution in that element
    \param[out] avg_velocity the velocity at the position (average in element)
    \section Description
    We will return the average velocity in the element to which the
    particle belongs to. To be more accurate we should really be
    interpolation and evaluating the basis functions at that point
    precisely. Doing that is more expensive and harder to code, so we
    won't be doing velocity interpolation. In theory we could do all
    this directly on the GPU (to avoid some memory transfer).
  */

  // Skip if the particle does not belong to anything
  if (el==-1){
    for(int alpha=0; alpha<D; alpha++){ avg_velocity[alpha] = 0;}
  }

  else{
    // Copy the velocity locally
#ifdef USE_CPU
    for(int k=0; k<(1+D)*N_s; k++){ solution[k] = U[el*N_F*N_s+k];}
#elif USE_GPU
    cudaMemcpy(solution, &U[el*N_F*N_s], (1+D)*N_s*sizeof(scalar), cudaMemcpyDeviceToHost);
#endif

    // Separate the fields
    scalar rho;
    for(int i=0; i<N_s; i++){
      rho               = solution[0*N_s+i];
      for(int alpha=0; alpha<D; alpha++){ avg_velocity[alpha]  += solution[(1+alpha)*N_s+i]/rho;}
    }

    // average the velocities
    for(int alpha=0; alpha<D; alpha++){ avg_velocity[alpha] =  avg_velocity[alpha]/N_s;}
  }  
}

//==========================================================================
//
//  Internal functions
//
//==========================================================================

//==========================================================================
void get_verts_of_element(const int e, const int nvert, const fullMatrix<scalar> XYZNodes, scalar* verts){
  /*!
    \brief Get the vertex coordinates of an element
    \param[in] e element index
    \param[in] nvert number of vertices in the polygon
    \param[in] XYZNodes matrix of all the nodal coordinates
    \param[out] verts array containing the coordinates of the polygon's vertices.
  */
  for(int alpha=0; alpha<D; alpha++){
    for(int k=0; k<nvert; k++){
      verts[alpha*nvert+k] = XYZNodes(k, e*D+alpha);
    }
  }
}

//==========================================================================
int pnpoly(int nvert, scalar *verts, scalar* position){
  /*!
    \brief Find out if a point is in a polygon
    \param[in] nvert number of vertices in the polygon
    \param[in] verts array containing the coordinates of the polygon's vertices.
    \param[in] position coordinates of the test point.
    \return c 1 if the point is in the polygon, 0 if not
    \section Description
    The 1D portion is pretty straightforward.
    For the 2D part of this function: 
    It is stolen from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    Uses ray tracing combined with the Jordan curve theorem.
    Works for 2D polygon (convex, concave, holes, etc)
    Copyright (c) 1970-2003, Wm. Randolph Franklin

    Works with:
    scalar* vertx = new scalar[N];
    vertx[0] = 0;  vertx[1] = 1;  vertx[2] = 1;  vertx[3] = 0;
    scalar* verty = new scalar[N];
    verty[0] = 0;  verty[1] = 0;  verty[2] = 1;  verty[3] = 1;
    printf("1: %d (1 expected)\n", pnpoly(4, vertx, verty, 0.1, 0.5));
    printf("1: %d (0 expected)\n", pnpoly(4, vertx, verty, 0.1,-0.5));
    OR
    scalar x1[] = { 0, 1, 1, 0 };
    scalar y1[] = { 0, 0, 1, 1 };
    printf("1: %d (1 expected)\n", pnpoly(4, x1, y1, 0.8, 0.8));
    printf("1: %d (0 expected)\n", pnpoly(4, x1, y1, -0.8, -0.8));
  */

#ifdef ONED
  scalar min = MIN(verts[0],verts[1]);
  scalar max = MAX(verts[0],verts[1]);
  if ((min < position[0]) && (position[0] < max)){ return 1;}
  return 0;
  
#elif TWOD
  scalar* vertx = &verts[0];
  scalar* verty = &verts[nvert];
  scalar testx = position[0];
  scalar testy = position[1];
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  vertx = NULL; verty = NULL;
  return c;
#endif
}
