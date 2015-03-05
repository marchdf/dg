/*!
  \file timers.cc
  \brief Function definitions for the TIMERS class.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \ingroup timers
*/
#include "timers.h"

void TIMERS::set_timer_names(){
  /*!
    \brief Set the names to all the timers
  */

  _names.push_back("main time");               // timer 0

  // RK routines
  _names.push_back("RK time");                 // timer 1
  _names.push_back("output time");             // timer 2
  _names.push_back("solver time");             // timer 3
  _names.push_back("new_dt time");             // timer 4

  // DG routines
  _names.push_back("DG time");                 // timer 5
  _names.push_back("mapToFace time");          // timer 6
  _names.push_back("rflctiveBoundary time");   // timer 7
  _names.push_back("collocationU time");       // timer 8
  _names.push_back("collocationdU time");      // timer 9
  _names.push_back("collocationUF time");      // timer 10
  _names.push_back("evaluate_sf time");        // timer 11
  _names.push_back("evaluate_q time");         // timer 12
  _names.push_back("redistribute_sf time");    // timer 13
  _names.push_back("gemm_s time");             // timer 14
  _names.push_back("gemm_f time");             // timer 15
  _names.push_back("redistribute_q time");     // timer 16
  _names.push_back("gemm_q time");             // timer 17
  _names.push_back("mapToElement time");       // timer 18
  _names.push_back("addSFQ time");             // timer 19

  // Limiting routines
  _names.push_back("limiting time");           // timer 20
  _names.push_back("sensors time");            // timer 21

  // Communication routines
  _names.push_back("communication time");      // timer 22
  
  // Output routines
  _names.push_back("format_output time");      // timer 23
  _names.push_back("write_output time");       // timer 24
  _names.push_back("format_sensor time");      // timer 25
  _names.push_back("write_sensor time");       // timer 26
  _names.push_back("write_particles time");    // timer 27
  
  // Lagrange particle routines
  _names.push_back("advect_particles time");   // timer 28

  // More timers
  _names.push_back("comm_packaging time");     // timer 29
  _names.push_back("comm_memcpy time");        // timer 30
  _names.push_back("comm_sendrecv time");      // timer 31
  _names.push_back("comm_unpackaging time");   // timer 32
}

void TIMERS::start_timer(int k){
  /*!
    \brief Start a timer
    \param[in] k the timer id to start
  */

  _starters[k] = std::clock();
  _counters[k]++;
}


void TIMERS::stop_timer(int k){
  /*!
    \brief Stop a timer
    \param[in] k the timer id to stop
  */

  // If you use a GPU, you need to wait for the kernel to end to time
  // it. This is because kernels are launched asynchronously by the
  // CPU (though executed sequentially on the GPU if they are placed
  // in the same stream).
#ifdef USE_GPU
  cudaDeviceSynchronize();
#endif
  
  _times[k] += ( std::clock() - _starters[k] )/ (double) CLOCKS_PER_SEC;  
}

void TIMERS::print_timers(){
  /*!
    \brief Print all the timers to a file
  */

  fprintf(_ofile, "# timers for this run: time in seconds, number of calls\n");
  for(int k=0;k<_N;k++){
    fprintf(_ofile,"# timer %3i: %s\n%20.16E\t%20d\n",k,_names[k].c_str(),_times[k],_counters[k]);
  }
}
