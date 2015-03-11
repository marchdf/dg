##
## \file tacc.mk
## \brief Makefile options for TACC specific paths
## \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
## \license This project is released under the GNU Public License. See LICENSE.
## \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
## \ingroup makefiles
##
## TACC is the Texas Advanced Computing Center at the U. of
## Texas-Austin. Through XSEDE we received computing time on the
## Stampede cluster. This makefile adjusts some paths for that
## cluster.
##
## @cond

CUDA_INSTALL_PATH := $(TACC_CUDA_DIR)
CXX               := icc

BLAS_INC	  := -I$(TACC_MKL_INC)
BLAS_LIB	  := -L$(TACC_MKL_LIB)
BLAS_LINKER       := -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

LAPACK_INC	  := -I$TACC_MKL_INC
LAPACK_LIB	  := -L$TACC_MKL_LIB
LAPACK_LINKER     := -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

GSL_LIB           := $(TACC_GSL_LIB)
GSL_INC	          := $(TACC_GSL_INC)

CUFLAGS           := -arch=compute_35 -code=sm_35

## @endcond