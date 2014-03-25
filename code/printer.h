/*!
  \file printer.h
  \brief Class deals with the solution output.
  \author Marc T. Henry de Frahan <marchdf@gmail.com>
  \defgroup printer Printer
  \ingroup printer
*/
#ifndef PRINTER_H
#define PRINTER_H
#include "fullMatrix.h"
#include "simpleMesh.h"
#include <scalar_def.h>
#include <stdlib.h>
#include <constants.h>
#include <vector>
#include <string>
#include <stdio.h>
#include "printer_kernels.h"
#ifdef USE_GPU
#include "misc_cuda.h"
#endif
class PRINTER{

 private:
  int _N_s;
  int _N_E;
  int _elem_type;
  scalar* _output;
  scalar* _d_output;
  std::vector<std::string> _names;
  std::vector<std::string> _fnames;
  simpleMesh &_m; 
    
 public:

  /*!
    \brief Constructor
    \param[in] N_s number of nodes per element
    \param[in] N_E number of elements
    \param[in] elem_type type of element to output
    \param[in] m mesh we are operating on
  */     
 PRINTER(int N_s, int N_E, int elem_type, simpleMesh &m): _N_s(N_s), _N_E(N_E), _elem_type(elem_type), _m(m){
#ifdef USE_CPU
    _output = new scalar[_N_s*_N_E*N_F];
#elif USE_GPU
    checkCuda(cudaMallocHost((void**)&_output, _N_s*_N_E*N_F*sizeof(scalar)));
    cudaMalloc((void**) &_d_output, _N_s*_N_E*N_F*sizeof(scalar));
#endif
  };

  /*! Destructor */
  ~PRINTER(){
#ifdef USE_CPU
    delete[] _output;
#elif USE_GPU
    cudaFreeHost(_output);
    cudaFree(_d_output);
#endif
  }

  /*!\brief Set the filenames and field names for the PRINTER class*/
  void set_names();
  
  void print(scalar* U, const int step, const double time);
};
#endif