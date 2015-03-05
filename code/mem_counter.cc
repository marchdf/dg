/*!
  \file mem_counter.cc
  \brief Functions for the MEM_COUNTER class
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#include "mem_counter.h"

void MEM_COUNTER::outputCounters(){
  /*!
    \brief Format/output bytes counter (from http://thomasfischer.biz/how-to-format-bytes-as-a-human-readable-string-without-branching/)
  */

  // Initializations
  char tmp[128] = "";
  const char *si_prefix[] = { "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB" };
  const int base = 1024;
  unsigned long bytes = 0;
  int c = 0;

  // CPU counter
  bytes = _counter_CPU;
  if (bytes==0){std::cout << "Estimated CPU allocations 0 bytes"<<std::endl;}
  else{
    c = std::min((int)(log((double)bytes)/log((double)base)), (int)sizeof(si_prefix) - 1);
    sprintf(tmp, "%1.2f %s", bytes / pow((double)base, c), si_prefix[c]);
    std::cout << "Estimated CPU allocations "<<std::string(tmp) << " (" << bytes << " bytes)"<<std::endl;
  }

  bytes = _counter_GPU;
  if (bytes==0){std::cout << "Estimated GPU allocations 0 bytes"<<std::endl;}
  else{
    c = std::min((int)(log((double)bytes)/log((double)base)), (int)sizeof(si_prefix) - 1);
    sprintf(tmp, "%1.2f %s", bytes / pow((double)base, c), si_prefix[c]);
    std::cout << "Estimated GPU allocations "<<std::string(tmp) << " (" << bytes << " bytes)"<<std::endl;
  }
}
