# Makefile options for TACC specific paths

CUDA_INSTALL_PATH := $(TACC_CUDA_DIR)
CXX               := icc

BLAS_INC	  := -I$(TACC_MKL_INC)
BLAS_LIB	  := -L$(TACC_MKL_LIB)
BLAS_LINKER       := -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

LAPACK_INC	  := -I$TACC_MKL_INC
LAPACK_LIB	  := -L$TACC_MKL_LIB
LAPACK_LINKER     := -Wl,-rpath,$TACC_MKL_LIB -L$TACC_MKL_LIB -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

GSL_LIB           := $(TACC_GSL_LIB)
GSL_INC	          := $(TACC_GSL_INC)

CUFLAGS           := -arch=compute_20 -code=sm_20