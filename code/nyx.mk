##
## \file nyx.mk
## \brief Makefile options for NYX/CAC UMICH specific paths
## \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
## \license
## \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
##

CUDA_INSTALL_PATH := $(CUDA_ROOT)

CXX               := icc

#BLAS_LINKER       := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
BLAS_LINKER       := -openmp -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#

#LAPACK_LINKER     := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
LAPACK_LINKER       := -openmp -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#

GSL_INC	          := $(GSL_ROOT)/include
GSL_LIB           := $(GSL_ROOT)/lib

