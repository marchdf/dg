##
## \file flux.mk
## \brief Makefile options for FLUX ARC UMICH specific paths
## \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
## \license This project is released under the GNU Public License. See LICENSE.
## \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
## \ingroup makefiles
## 
## Flux is University of Michigan's cluster (formely Nyx) and operated
## by Advanced Research Computing (formely Center for Advanced
## Computing). This makefile adjusts some paths for that cluster.
##
## @cond

CUDA_INSTALL_PATH := /gpfs/gpfs0/software/rhel72/packages/cuda/8.0

CXX               := g++

BLAS_LIB          := /gpfs/gpfs0/software/rhel72/packages/openblas/0.2.19/gcc-5.4.0/lib
BLAS_INC          := /gpfs/gpfs0/software/rhel72/packages/openblas/0.2.19/gcc-5.4.0/include
BLAS_LINKER       := /gpfs/gpfs0/software/rhel72/packages/openblas/0.2.19/gcc-5.4.0/lib/libopenblas.so.0 #/usr/lib64/libblas.so.3.4.2 #you can control the number of threads: export MKL_

LAPACK_LIB        := /gpfs/gpfs0/software/rhel72/packages/lapack/3.6.1/gcc-5.4.0/lib
LAPACK_LINKER       := /gpfs/gpfs0/software/rhel72/packages/lapack/3.6.1/gcc-5.4.0/lib/liblapack.so.3.6.1 #/usr/lib64/liblapack.so.3.4.2 #you can control the number of threads: export 

GSL_INC	          := $(GSL_ROOT)/include
GSL_LIB           := $(GSL_ROOT)/lib



## @endcond
