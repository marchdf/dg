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
//#include <cstring> // for some reason I need to include this or I get a memcpy not defined

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
#else
class COMMUNICATOR{
 public:
  COMMUNICATOR(int a, int b, simpleMesh &m){};
  void CommunicateGhosts(int Nfields, scalar* U){};
};
#endif

#endif // COMMUNICATOR_H

