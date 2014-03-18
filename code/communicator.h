//
// Class to communicate the faces between cpus and gpus
//

#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#ifdef USE_MPI
#include <scalar_def.h>
#include "mpi.h"
#include "simpleMesh.h"
#include "misc.h"

//================================================================================
//
// Communicator class with CPUs only
//
//================================================================================
#ifdef USE_CPU
class COMMUNICATOR
{
 private:
  int _N_ghosts;
  int _N_s;
  int* _ghostElementSend;
  int* _ghostElementRecv;
  MPI_Status *_status;
  MPI_Request *_request;

 public:
  // constructors
  COMMUNICATOR(int N_ghosts, int N_s, simpleMesh &m): _N_ghosts(N_ghosts), _N_s(N_s){
    // Initialize MPI things
    _status = new MPI_Status [2*N_ghosts];
    _request = new MPI_Request [2*N_ghosts];

    _ghostElementSend = new int[3*N_ghosts];
    _ghostElementRecv = new int[3*N_ghosts];
    // I thought this would help the double free... it didn't
    /* hostDeepCopyArray(m.getGhostElementSend(), _ghostElementSend, 3*N_ghosts); */
    /* hostDeepCopyArray(m.getGhostElementRecv(), _ghostElementRecv, 3*N_ghosts); */
    memcpy(_ghostElementSend, m.getGhostElementSend() , 3*N_ghosts*sizeof(int));
    memcpy(_ghostElementRecv, m.getGhostElementRecv() , 3*N_ghosts*sizeof(int));
  };

  // destructor
  ~COMMUNICATOR(){
    delete[] _ghostElementSend;
    delete[] _ghostElementRecv;
  };

  // Communicate the ghost elements
  void CommunicateGhosts(int Nfields, scalar* U);
};


//================================================================================
//
// Communicator class with GPUs only
//
//================================================================================
#elif USE_GPU
#include <cuda_runtime_api.h>
#include "misc_cuda.h"
class COMMUNICATOR
{
 private:
  int _N_ghosts;
  int _N_s;

  MPI_Status *_status;
  MPI_Request *_request;

  int* _ghostElementSend;
  int* _ghostElementRecv;
  scalar* _h_bufferSend;
  scalar* _h_bufferRecv;

  int* _d_ghostElementSend; // these will only hold the element number, no need to tags or dest
  int* _d_ghostElementRecv;
  scalar* _d_buffer;

 public:
  // constructors
  COMMUNICATOR(int N_ghosts, int N_s, simpleMesh &m): _N_ghosts(N_ghosts), _N_s(N_s){
    // Initialize MPI things
    _status = new MPI_Status [2*N_ghosts];
    _request = new MPI_Request [2*N_ghosts];

    // Allocate/initialize CPU arrays
    _ghostElementSend = new int[3*N_ghosts]; 
    _ghostElementRecv = new int[3*N_ghosts];     
    memcpy(_ghostElementSend, m.getGhostElementSend() , 3*N_ghosts*sizeof(int));
    memcpy(_ghostElementRecv, m.getGhostElementRecv() , 3*N_ghosts*sizeof(int));
    checkCuda(cudaMallocHost((void**)&_h_bufferSend, _N_s*_N_ghosts*N_F*sizeof(scalar)));
    checkCuda(cudaMallocHost((void**)&_h_bufferRecv, _N_s*_N_ghosts*N_F*sizeof(scalar)));
    
    // Allocate/initialize GPU arrays
    cudaMalloc((void**) &_d_ghostElementSend, _N_ghosts*sizeof(int));
    cudaMalloc((void**) &_d_ghostElementRecv, _N_ghosts*sizeof(int));
    cudaMalloc((void**) &_d_buffer, _N_s*_N_ghosts*N_F*sizeof(int));

    // Copy necessary data to GPU memory (using a tmp array)
    int* tmp = new int[_N_ghosts];
    for(int k=0; k<_N_ghosts; k++){tmp[k] = _ghostElementSend[k*3+0];}
    cudaMemcpy(_d_ghostElementSend, tmp , _N_ghosts*sizeof(int) , cudaMemcpyHostToDevice);
    for(int k=0; k<_N_ghosts; k++){tmp[k] = _ghostElementRecv[k*3+0];}
    cudaMemcpy(_d_ghostElementRecv, tmp , _N_ghosts*sizeof(int) , cudaMemcpyHostToDevice);
    delete[] tmp;
   
  };

  // destructor
  ~COMMUNICATOR(){
    delete[] _ghostElementSend;
    delete[] _ghostElementRecv;
    cudaFreeHost(_h_bufferSend);
    cudaFreeHost(_h_bufferRecv);

    cudaFree(_d_ghostElementSend);
    cudaFree(_d_ghostElementRecv);
    cudaFree(_d_buffer);
  };

  // Communicate the ghost elements
  void CommunicateGhosts(int Nfields, scalar* U);
};
#endif // on MPI processor type

//================================================================================
//
// If there is no MPI, make an empty class
//
//================================================================================
#else
class COMMUNICATOR{
 public:
  COMMUNICATOR(int a, int b, simpleMesh &m){};
  void CommunicateGhosts(int Nfields, scalar* U){};
};
#endif

#endif // COMMUNICATOR_H

