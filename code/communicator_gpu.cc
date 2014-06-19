/*!
  \file communicator_gpu.cc
  \brief GPU implementation of CommunicateGhosts
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \ingroup gpugroup
*/

#ifdef USE_MPI
#ifdef USE_GPU
#include <communicator.h>
#include <communicator_gpu_kernels.h>

void COMMUNICATOR::CommunicateGhosts(int Nfields, scalar* U){
  /*!
    \brief Communicate the ghosts (GPU implementation)
    \param[in] Nfields number of fields to act on
    \param[in] U solution to communicate
    \ingroup gpugroup
    \section Description
    This function communicates the elements which are not in my
    partition. Basically you send the elements of U that other
    partitions need and you receive the elements from other
    partitions that you will need. You store these in the back
    columns of U.
    
    Since this is the GPU version, we do it using host buffers, which
    involves more steps:     
    1) copy data from device U to device buffer
    2) copy device buffer to pinned host send buffer
    3) loop on k for immediate send/receives
    4) wait end of communications
    5) copy pinned host recv buffer to device buffer
    6) kernel to copy data from device buffer to device U
    
    Had to declare Nfields when we operate on just one field instead of N_F

    See also \link communicator_cpu.cc CommunicateGhosts() \endlink the CPU implementation
  */

  _timers.start_timer(22);
  
  // wait until every process gets here
  MPI_Barrier(MPI_COMM_WORLD);

  // 1) copy data from device U to device buffer
  Lpackager(_N_s, _N_ghosts, Nfields, _d_ghostElementSend, _d_buffer, U);
  
  // 2) copy device buffer to pinned host send buffer
  cudaMemcpy(_h_bufferSend, _d_buffer, _N_s*_N_ghosts*Nfields*sizeof(scalar), cudaMemcpyDeviceToHost);
  
  // 3) loop on k for immediate send/receives
  int e, dest, source, tag;  
  for(int k=0; k<_N_ghosts; k++){

    // Send info 
    e    = _ghostElementSend[k*3+0]; // local index of U to send
    dest = _ghostElementSend[k*3+1]; // destination processor
    tag  = _ghostElementSend[k*3+2]; // global idx of element
    MPI_Isend(&_h_bufferSend[k*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    
    // Recv info
    e      = _ghostElementRecv[k*3+0]; // local index of U to recv (back columns)
    source = _ghostElementRecv[k*3+1]; // destination processor
    tag    = _ghostElementRecv[k*3+2]; // global idx of element
    MPI_Irecv(&_h_bufferRecv[k*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
  }

  // 4) wait end of communications
  MPI_Waitall(2*_N_ghosts, _request, _status);

  // 5) copy pinned host recv buffer to device buffer
  cudaMemcpy(_d_buffer, _h_bufferRecv, _N_s*_N_ghosts*Nfields*sizeof(scalar), cudaMemcpyHostToDevice);
  
  // 6) kernel to copy data from device buffer to device U
  Lunpackager(_N_s, _N_ghosts, Nfields, _d_ghostElementRecv, _d_buffer, U);

  _timers.stop_timer(22);
}
#endif
#endif
