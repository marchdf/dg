/*!
  \file timers.cc
  \brief Function definitions for the TIMERS class.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
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

  //Phil's viscous DG routines
  _names.push_back("VIS_CUGUF time");           // timer 33
  _names.push_back("Unused");       // timer 34
  _names.push_back("Unused");         // timer 35
  _names.push_back("Unused");   // timer 36
  _names.push_back("VIS_Uhat_to_Sigma time");                  // timer 37
  _names.push_back("VIS_UhCommon_to_Sigma time");              // timer 38
  _names.push_back("VIS_adjust_sf time");                      // timer 39
  _names.push_back("VIS_adjust_q time");                       // timer 40
  
  //Phil's icb reconstruction stuff for advective flux components
  _names.push_back("ICB_BuildSolution time");        //timer 41
  _names.push_back("Unused");        //timer 42

  //Phil's additions for artificial dissipation (like regular viscous flow plus Cparam routine)
  _names.push_back("AD_Cparams time");                        // timer 43
  _names.push_back("AD_Uhat_to_GradCommonHalf time");         // timer 44
  _names.push_back("AD_GradCommonHalf_to_GradCommon time");   // timer 45
  _names.push_back("AD_adjust_sf time");                      // timer 46
  _names.push_back("AD_adjust_q time");                       // timer 47
  _names.push_back("AD_Uhat_to_UhCommonHalf time");           // timer 48
  _names.push_back("AD_UhCommonHalf_to_UhCommon time");       // timer 49
  _names.push_back("AD_Uhat_to_Sigma time");                  // timer 50
  _names.push_back("AD_UhCommon_to_Sigma time");              // timer 51

  //A parallel efficieny check:
  _names.push_back("Ghost Communication barrier time");        // timer 52
  _names.push_back("C-parameter AllReduce");                   //timer 53

  //Initialization timers:
  _names.push_back("Initialization: PSIxR for diffusion");     //timer 54
  _names.push_back("Initialization: PSIxR for advection");    // timer 55
  _names.push_back("Initialization: SigFaceMatrices");        //timer 56
  _names.push_back("Initialization: BuildSigMatrices and BuildAuxHatMatrices"); //timer 57
  _names.push_back("Initialization: Copy matrix to pointer"); //timer 58
  _names.push_back("Inside PSIxR_oneFace: dets and precords"); //timer 59 
  _names.push_back("Inside PSIxR_oneFace: Reco Cords and Psi"); //timer 60
  _names.push_back("Inside PSIxR_oneFace: Matrix formation and solve"); //timer 61
  _names.push_back("Inside PSIxR_biased: dets and precords"); //timer 62
  _names.push_back("Inside PSIxR_biased: NodalEqCord"); //timer 63
  _names.push_back("Inside PSIxR_biased: Reco Cords and Psi"); //timer 64
  _names.push_back("Inside PSIxR_biased: Some basis stuff, rg"); //timer 65
  _names.push_back("Inside PSIxR_biased: Matrix formation and solve"); //timer 66
  
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

double TIMERS::print_single_timer(int k)
{
  return _times[k];
}
