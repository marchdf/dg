//
// Class to communicate the faces between cpus and gpus
//

#ifndef COMMUNICATOR_ELEMENTS_H
#define COMMUNICATOR_ELEMENTS_H

#ifdef USE_MPI
#include <scalar_def.h>
#include "mpi.h"

class COMMUNICATOR_ELEMENTS
{
 private:
  int _N_ghosts;
  int _N_s;
  MPI_Status *_status;
  MPI_Request *_request;

 public:
  // constructors
 COMMUNICATOR_ELEMENTS(){};
 COMMUNICATOR_ELEMENTS(int N_ghosts, int N_s): _N_ghosts(N_ghosts), _N_s(N_s){
    // Initialize MPI things
    _status = new MPI_Status [2*N_ghosts];
    _request = new MPI_Request [2*N_ghosts];
  };

  // destructor
  ~COMMUNICATOR_ELEMENTS(){};

  // Communicate the ghost elements
  void CommunicateGhosts(int Nfields, int* ghostElementSend, int* ghostElementRecv, scalar* A);

};
#else
class COMMUNICATOR_ELEMENTS{
 public:
  COMMUNICATOR_ELEMENTS(){};
  COMMUNICATOR_ELEMENTS(int a, int b){};
};
#endif

#endif // COMMUNICATOR_ELEMENTS_H

