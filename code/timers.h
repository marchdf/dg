/*!
  \file timers.h
  \brief Class which deals with the timers 
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup timers Timers
  \ingroup timers
*/
#ifndef TIMERS_H
#define TIMERS_H

#include <ctime>
#include <vector>
#include <string>
#include <stdio.h>

class TIMERS{

 private:
  int _N;                          // number of timers
  double* _times;                  // holds the times
  int* _counters;                  // holds the number of timer calls (i.e. function calls)
  std::clock_t* _starters;         // holds the timer start times
  std::vector<std::string> _names; // holds the timer names
  FILE* _ofile;
  std::string _fname;
  void set_timer_names();
  
 public:

  /*!
    \brief Constructor
    \param[in] myid processor id
   */     
  TIMERS(int myid){

    // Set the names and number of the different timers
    set_timer_names();

    // Initialize the data
    _N = _names.size();
    _times = new double[_N]; for(int k=0;k<_N;k++){_times[k]=0.0;}
    _counters = new int[_N]; for(int k=0;k<_N;k++){_counters[k]=0;}
    _starters = new std::clock_t[_N];

    // Set the filename, open the file
    _fname = "timers.dat";
#ifdef USE_MPI
    char myidstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(myidstr, "%d", myid);
    _fname += myidstr;
#endif
    _ofile = fopen(_fname.c_str(),"w");
  }

  /*! Destructor */
  ~TIMERS(){
    delete[] _times;
    delete[] _counters;
    delete[] _starters;
    fclose(_ofile);    
  }

  void start_timer(int k);
  void stop_timer(int k);
  void print_timers();
};

#endif
