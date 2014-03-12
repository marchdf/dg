#ifdef USE_MPI
#include <communicator_faces.h>

void COMMUNICATOR_FACES::mapGhostFace(int* ghostInterfaces, scalar* UF){
  /* Communicate the ghost edges in the mesh to and from the other
     partitions

     Just CPU version for now
  */

  // wait until every process gets here
  MPI_Barrier(MPI_COMM_WORLD);

  // attach the buffer
  MPI_Buffer_attach( _buf, _bufsize );
  
  int t, dest_source, tag,id;
  //MPI_Comm_rank(MPI_COMM_WORLD,&id);
  //printf("_M_ghosts=%i\n",_M_ghosts);
  for(int k = 0; k < _M_ghosts; k++){
    t           = ghostInterfaces[k*3+0];
    dest_source = ghostInterfaces[k*3+1];
    tag         = ghostInterfaces[k*3+2];

    //printf("CPU id %i is sending/receiving interface %i to proc %i with tag %i\n",id,t,dest_source,tag);

    // VERSION 1: non-blocking send and receive
    // Non-blocking send of my interface
    //MPI_Isend(&UF[t*N_F*2*M_s], 1, _strided, dest_source, tag, MPI_COMM_WORLD, &_request[2*k+0]);
    // Non-blocking receive of the other side of that interface
    //MPI_Irecv(&UF[t*N_F*2*M_s+M_s], 1, _strided, dest_source, tag, MPI_COMM_WORLD, &_request[2*k+1]);

    // VERSION 2: buffered send and non-blocking receive
    MPI_Bsend(&UF[t*N_F*2*_M_s+_M_s], 1, _strided, dest_source, tag, MPI_COMM_WORLD);
    MPI_Irecv(&UF[t*N_F*2*_M_s+_M_s], 1, _strided, dest_source, tag, MPI_COMM_WORLD, &_request[k]);
  }

  // Wait for all communications to end
  MPI_Waitall(_M_ghosts, _request, _status);
  //printf("Sending/Receiving done by cpu %i\n",id);

  // Detach the buffer
  MPI_Buffer_detach( &_buf, &_bufsize);

  // Not necessary I think? wait until every process gets here
  //MPI_Barrier(MPI_COMM_WORLD); 
}
#endif
