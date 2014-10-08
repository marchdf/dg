/*!
  \file communicator_cpu.cc
  \brief CPU implementation of CommunicateGhosts
  \copyright Copyright (C) 2014, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup cpugroup
*/

#ifdef USE_MPI
#ifdef USE_CPU
#include <communicator.h>

void COMMUNICATOR::CommunicateGhosts(int Nfields, scalar* U){
  /*!
    \brief Communicate the ghosts (CPU implementation)
    \param[in] Nfields number of fields to act on
    \param[in] U solution to communicate
    \ingroup cpugroup
    \section Description
    This function communicates the elements which are not in my
    partition. Basically you send the elements of U that other
    partitions need and you receive the elements from other
    partitions that you will need. You store these in the back
    columns of U.
    
    Had to declare Nfields when we operate on just one field instead of N_F

    See also \link communicator_gpu.cc CommunicateGhosts() \endlink the GPU implementation
  */

  _timers.start_timer(22);
  
  // wait until every process gets here
  MPI_Barrier(MPI_COMM_WORLD);

  _timers.start_timer(31);
  int e, dest, source, tag;
  for(int k=0; k<_N_ghosts; k++){

    // Send info 
    e    = _ghostElementSend[k*3+0]; // local index of U to send
    dest = _ghostElementSend[k*3+1]; // destination processor
    tag  = _ghostElementSend[k*3+2]; // global idx of element
    MPI_Isend(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    
    // Recv info
    e      = _ghostElementRecv[k*3+0]; // local index of U to recv (back columns)
    source = _ghostElementRecv[k*3+1]; // destination processor
    tag    = _ghostElementRecv[k*3+2]; // global idx of element
    MPI_Irecv(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
  }
  _timers.stop_timer(31);
  
  // Wait for communication to end
  MPI_Waitall(2*_N_ghosts, _request, _status);

  _timers.stop_timer(22);
}
#endif
#endif
