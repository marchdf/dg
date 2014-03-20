/*!
  \file communicator.h
  \class COMMUNICATOR communicator.h
  \brief Class to communicate the faces between cpus and gpus. 
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup cpugroup Group CPU
  \defgroup gpugroup Group GPU
  \section Description
  Class to communicate the faces between cpus and gpus.
*/

#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include "simpleMesh.h"
#ifdef USE_MPI
#include <scalar_def.h>
#include "mpi.h"
#include "misc.h"
#ifdef USE_GPU
#include <cuda_runtime_api.h>
#include "misc_cuda.h"
#endif

//================================================================================
//
// Communicator class for MPI
//
//================================================================================
class COMMUNICATOR
{
 private:
  int _N_ghosts;
  int _N_s;

  MPI_Status *_status;
  MPI_Request *_request;

  int* _ghostElementSend;
  int* _ghostElementRecv;

#ifdef USE_GPU
  scalar* _h_bufferSend;
  scalar* _h_bufferRecv;
  int* _d_ghostElementSend; // these will only hold the element number, no need to tags or dest
  int* _d_ghostElementRecv;
  scalar* _d_buffer;
#endif
  
 public:
  /*!
    \brief Constructor sets default number of ghosts and nodes in an element
    \param N_ghosts number of ghosts
    \param N_s number of nodes in an element
    \param m mesh operated on
  */     
  COMMUNICATOR(int N_ghosts, int N_s, simpleMesh &m): _N_ghosts(N_ghosts), _N_s(N_s){
    // Initialize MPI things
    _status = new MPI_Status [2*N_ghosts];
    _request = new MPI_Request [2*N_ghosts];

    // Allocate/initialize CPU arrays
    _ghostElementSend = new int[3*N_ghosts]; 
    _ghostElementRecv = new int[3*N_ghosts];     
    memcpy(_ghostElementSend, m.getGhostElementSend() , 3*N_ghosts*sizeof(int));
    memcpy(_ghostElementRecv, m.getGhostElementRecv() , 3*N_ghosts*sizeof(int));

#ifdef USE_GPU
    // Allocate/initialize CPU/GPU arrays
    checkCuda(cudaMallocHost((void**)&_h_bufferSend, _N_s*_N_ghosts*N_F*sizeof(scalar)));
    checkCuda(cudaMallocHost((void**)&_h_bufferRecv, _N_s*_N_ghosts*N_F*sizeof(scalar)));
    cudaMalloc((void**) &_d_ghostElementSend, _N_ghosts*sizeof(int));
    cudaMalloc((void**) &_d_ghostElementRecv, _N_ghosts*sizeof(int));
    cudaMalloc((void**) &_d_buffer, _N_s*_N_ghosts*N_F*sizeof(scalar));

    // Copy necessary data to GPU memory (using a tmp array)
    int* tmp = new int[_N_ghosts];
    for(int k=0; k<_N_ghosts; k++){tmp[k] = _ghostElementSend[k*3+0];}
    cudaMemcpy(_d_ghostElementSend, tmp , _N_ghosts*sizeof(int) , cudaMemcpyHostToDevice);
    for(int k=0; k<_N_ghosts; k++){tmp[k] = _ghostElementRecv[k*3+0];}
    cudaMemcpy(_d_ghostElementRecv, tmp , _N_ghosts*sizeof(int) , cudaMemcpyHostToDevice);
    delete[] tmp;
#endif
  };

  /*! Destructor */
  ~COMMUNICATOR(){
    delete[] _ghostElementSend;
    delete[] _ghostElementRecv;

#ifdef USE_GPU
    cudaFreeHost(_h_bufferSend);
    cudaFreeHost(_h_bufferRecv);
    cudaFree(_d_ghostElementSend);
    cudaFree(_d_ghostElementRecv);
    cudaFree(_d_buffer);
#endif
  };

  void CommunicateGhosts(int Nfields, scalar* U);
};

//================================================================================
//
// If there is no MPI, make an empty class
//
//================================================================================
#else
class COMMUNICATOR{
 public:
  /*! \brief Empty constructor
  */     
  COMMUNICATOR(int a, int b, simpleMesh &m){};
  void CommunicateGhosts(int Nfields, scalar* U){};
};
#endif

#endif // COMMUNICATOR_H

