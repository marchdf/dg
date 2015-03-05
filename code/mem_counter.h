/*!
  \file mem_counter.h
  \brief Class deals with counting the number of bytes allocated.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef MEM_COUNTER_H
#define MEM_COUNTER_H
#include <string>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <scalar_def.h>

class MEM_COUNTER{

 private:
  size_t _counter_CPU; // counter in bytes (CPU)
  size_t _counter_GPU; // counter in bytes (GPU)

 public:
  /*!\brief Constructor defaults*/ 
  MEM_COUNTER(){_counter_CPU=0;_counter_GPU=0;};
  void addToCPUCounter(size_t new_bytes) {/*!\brief Add bytes to CPU counter*/ _counter_CPU+=new_bytes;}
  void addToGPUCounter(size_t new_bytes) {/*!\brief Add bytes to GPU counter*/ _counter_GPU+=new_bytes;}
  void outputCounters();
};
#endif
