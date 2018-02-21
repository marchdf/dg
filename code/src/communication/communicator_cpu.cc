/*!
  \file communicator_cpu.cc
  \brief CPU implementation of CommunicateGhosts
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup cpugroup
*/

#ifdef USE_MPI
#ifdef USE_CPU
#include "communicator.h"

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
  _timers.start_timer(52);
  MPI_Barrier(MPI_COMM_WORLD);
  _timers.stop_timer(52);

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

void COMMUNICATOR::CommunicateCords(scalar* XYZNodes)
{
  /*!
    \brief communicate solution node coordinates across partitions
   */
  // wait until every process gets here
  MPI_Barrier(MPI_COMM_WORLD);

  int e, dest, source, tag;
  for(int k=0; k<_N_ghosts; k++){
    
    // Send info 
    e    = _ghostElementSend[k*3+0]; // local index of U to send
    dest = _ghostElementSend[k*3+1]; // destination processor
    tag  = _ghostElementSend[k*3+2]; // global idx of element
    //printf("Send command: k=%d/%d, e=%d, dest=%d, tag=%d\n", k, _N_ghosts, e, dest, tag);
    MPI_Isend(&XYZNodes[e*D*_N_s], D*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
   //MPI_Isend(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    
    // Recv info
    e      = _ghostElementRecv[k*3+0]; // local index of U to recv (back columns)
    source = _ghostElementRecv[k*3+1]; // destination processor
    tag    = _ghostElementRecv[k*3+2]; // global idx of element
    //printf("Receive command: k=%d/%d, e=%d, source=%d, tag=%d\n", k, _N_ghosts, e, source, tag);
    MPI_Irecv(&XYZNodes[e*D*_N_s], D*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
    //MPI_Irecv(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
  }  
    // Wait for communication to end
    MPI_Waitall(2*_N_ghosts, _request, _status);
    
}

void COMMUNICATOR::CommunicateSensor(int* SensorTag_Aug)
{
  /*!
    \brief communicate sensor output across partition psuedo-boundaries
   */
  // wait until every process gets here
  _timers.start_timer(22);
  MPI_Barrier(MPI_COMM_WORLD);

  int e, dest, source, tag;
  for(int k=0; k<_N_ghosts; k++){
    
    // Send info 
    e    = _ghostElementSend[k*3+0]; // local index of U to send
    dest = _ghostElementSend[k*3+1]; // destination processor
    tag  = _ghostElementSend[k*3+2]; // global idx of element
    //printf("Send command: k=%d/%d, e=%d, dest=%d, tag=%d\n", k, _N_ghosts, e, dest, tag);
    MPI_Isend(&SensorTag_Aug[e], 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    //MPI_Isend(&XYZNodes[e*D*_N_s], D*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
   //MPI_Isend(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    
    // Recv info
    e      = _ghostElementRecv[k*3+0]; // local index of U to recv (back columns)
    source = _ghostElementRecv[k*3+1]; // destination processor
    tag    = _ghostElementRecv[k*3+2]; // global idx of element
    //printf("Receive command: k=%d/%d, e=%d, source=%d, tag=%d\n", k, _N_ghosts, e, source, tag);
    MPI_Irecv(&SensorTag_Aug[e], 1, MPI_INT, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
    //    MPI_Irecv(&XYZNodes[e*D*_N_s], D*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
    //MPI_Irecv(&U[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
  }  
    // Wait for communication to end
    MPI_Waitall(2*_N_ghosts, _request, _status);
    _timers.stop_timer(22);
}

#endif
#endif
