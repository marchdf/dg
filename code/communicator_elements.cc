#include <communicator_elements.h>

void COMMUNICATOR_ELEMENTS::CommunicateGhosts(int Nfields, int* ghostElementSend, int* ghostElementRecv, scalar* A){

  /* This function, used by the limiting procedure, communicates the
     elements which are not in my partition. Basically you send the
     elements of A that other partitions need and you receive the
     elements from other partitions that you will need. You store
     these in the back columns of A.
     
     Had to declare Nfields when we operate on just one field instead of N_F*/

  // wait until every process gets here
  MPI_Barrier(MPI_COMM_WORLD);
      
  int e, dest, source, tag;
  
  for(int k=0; k<_N_ghosts; k++){

    // Send info 
    e   = ghostElementSend[k*3+0]; // local index of A to send
    dest = ghostElementSend[k*3+1]; // destination processor
    tag = ghostElementSend[k*3+2]; // global idx of element
    MPI_Isend(&A[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, dest, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    
    // Recv info
    e      = ghostElementRecv[k*3+0]; // local index of A to recv (back columns)
    source = ghostElementRecv[k*3+1]; // destination processor
    tag    = ghostElementRecv[k*3+2]; // global idx of element
    MPI_Irecv(&A[e*Nfields*_N_s], Nfields*_N_s, MPI_SCALAR, source, tag, MPI_COMM_WORLD, &_request[2*k+1]);
  }

  // Wait for communication to end
  MPI_Waitall(2*_N_ghosts, _request, _status);

}
