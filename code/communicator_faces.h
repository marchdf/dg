//
// Class to communicate the faces between cpus and gpus
//
#ifndef COMMUNICATOR_FACES_H
#define COMMUNICATOR_FACES_H

#include <scalar_def.h>
#include "mpi.h"

class COMMUNICATOR_FACES
{
 private:
  int _M_ghosts;
  int _M_s;
  MPI_Status *_status;//[M_ghosts];
  MPI_Request *_request;//[M_ghosts];
  
  MPI_Datatype _strided; // for strided access

  int _bufsize;
  scalar *_buf;
  
 public:
  // constructors
 COMMUNICATOR_FACES(int M_ghosts,int M_s): _M_ghosts(M_ghosts), _M_s(M_s){

    // Initialize MPI things
    _status = new MPI_Status [M_ghosts];
    _request = new MPI_Request [M_ghosts];
    
    // make a strided datatype to access UF (only access one side of the interface)
    // from http://stackoverflow.com/questions/15483360/mpi-sending-segments-of-an-array
    // and http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Type_vector.html
    MPI_Type_vector (N_F, M_s, 2*M_s, MPI_SCALAR, &_strided);
    MPI_Type_commit(&_strided);

    // Initialize a buffer
    _bufsize =  M_ghosts * (MPI_BSEND_OVERHEAD + N_F*M_s) + MPI_BSEND_OVERHEAD;
    _buf = new scalar[_bufsize];
  };
  
  // destructor
  ~COMMUNICATOR_FACES(){
    //MPI_Type_free(&_strided); // Is this necessary? doesn't look like it
    if(_buf)       delete[] _buf;
  };

  // Communicate the faces
  void mapGhostFace(int* ghostInterfaces, scalar* UF);
  
};


#endif // COMMUNICATOR_H
