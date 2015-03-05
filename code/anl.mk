##
## \file anl.mk
## \brief Makefile options for ANL specific paths
## \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
## \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
##

CXX               := mpicxx
# If you want to use the xl compilers, uncomment these lines
#CXX               := mpixlcxx_r	
#CC                := mpixlc_r
#DEPFLAGS          := -qmakedep=gcc

# ALCF Netlib built with XLF - assuming no OpenMP
XLF_LIBS             := -L$(IBM_MAIN_DIR)/xlf/bg/11.1/bglib -lxlf90_r -lxlfmath -lxl -L$(IBM_MAIN_DIR)/xlsmp/bg/1.7/bglib -lxlomp_ser -lpthread
BLAS_LIB        := /soft/apps/current/BLAS
BLAS_LINKER     := -lblas_bgp $(XLF_LIBS)
# ESSL - assuming OpenMP but not using threaded BLAS (since it is slower for smaller operations)
#XLF_LIBS          := -L$(IBM_MAIN_DIR)/xlf/bg/11.1/bglib -lxlf90_r -lxlfmath -lxl -L$(IBM_MAIN_DIR)/xlsmp/bg/1.7/bglib -lxlsmp -lpthread
#BLAS_LIB          := /soft/apps/ESSL-4.4.1-1/lib
#BLAS_LINKER       := -lesslbg $(XLF_LIBS)

LAPACK_LIB        := /soft/apps/current/LAPACK
# ALCF Netlib built with XLF
LAPACK_LINKER     := -llapack_bgp
# ESSL - linker acts left-to-right so this will pick up ESSL LAPACK symbols first then grab Netlib ones if the app needs some that are missing in ESSL
#LAPACK_LINKER     := -lesslbg -llapack_bgp
# This is the wrong way to do things if one is calling the handful of functions in ESSL that have the same name (symbol) but different call syntax than LAPACK
# _gesvd is a function that seems to require LAPACK be linked first

# For ALCF/ANL HPCTW profiling
TRACE		  := -L/soft/apps/current/ibm-hpct/lib -lmpitrace -llicense -lgetarg -L/soft/apps/ibmcmp-jan2013/vac/bg/9.0/bglib -L/soft/apps/ibmcmp-jan2013/vacpp/bg/9.0/bglib -lxlopt -lxl -libmc++