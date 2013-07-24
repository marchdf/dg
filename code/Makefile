# defaults
ARCH=USE_CPU
PARA=NO_MPI
PREC=USE_DOUBLE
DIMS=ONED
PROB=MULTIFLUID
GAMM=GAMNCON
FLUX=ROE
LIEU=HOME
NFLD=4
MASS=0

# configuration
directories       := . quadratures
BINNAME           := dg1d

# paths to libs
CUDA_INSTALL_PATH := /usr/local/cuda
CUDA_SDK_PATH     := /home/marchdf/NVIDIA_GPU_Computing_SDK/C
ifeq ($(LIEU),NYX)
CUDA_INSTALL_PATH := $(CUDA_ROOT)
endif

# compilers
NVCC              := nvcc
CXX               := g++
DEPFLAGS          := -M
ifeq ($(LIEU),NYX)
CXX               := icc
endif
ifeq ($(PARA),USE_MPI)
CXX               := mpicxx
endif
ifeq ($(LIEU),ANL)
CXX               := mpicxx
#CXX               := mpixlcxx_r	
#CC                := mpixlc_r
DEPFLAGS          := -M #-qmakedep=gcc
endif

# defines
DEFINES           := -DHAVE_BLAS -DHAVE_LAPACK -D$(ARCH) -D$(LIEU) -D$(PREC) -D$(PROB) -D$(GAMM) -D$(FLUX) -D$(DIMS) -D$(PARA) -DN_Y=$(MASS) -DN_F=$(NFLD)
ifeq ($(DIMS),ONED)
DEFINES           := $(DEFINES) -DD=1 -DDF=1
endif
ifeq ($(DIMS),TWOD)
DEFINES           := $(DEFINES) -DD=2 -DDF=1
endif
ifeq ($(DIMS),THREED)
DEFINES           := $(DEFINES) -DD=3 -DDF=2
endif

# blas libs
BLAS_INC          := . # non-empty vars because XL doesn't like "-I -L "
BLAS_LIB          := .
BLAS_LINKER       := -lblas
ifeq ($(LIEU),NYX)
#BLAS_LINKER       := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
BLAS_LINKER       := -openmp -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#
endif
ifeq ($(LIEU),ANL)
# ALCF Netlib built with XLF - assuming no OpenMP
XLF_LIBS             := -L$(IBM_MAIN_DIR)/xlf/bg/11.1/bglib -lxlf90_r -lxlfmath -lxl -L$(IBM_MAIN_DIR)/xlsmp/bg/1.7/bglib -lxlomp_ser -lpthread
BLAS_LIB        := /soft/apps/current/BLAS
BLAS_LINKER     := -lblas_bgp $(XLF_LIBS)
# ESSL - assuming OpenMP but not using threaded BLAS (since it is slower for smaller operations)
#XLF_LIBS          := -L$(IBM_MAIN_DIR)/xlf/bg/11.1/bglib -lxlf90_r -lxlfmath -lxl -L$(IBM_MAIN_DIR)/xlsmp/bg/1.7/bglib -lxlsmp -lpthread
#BLAS_LIB          := /soft/apps/ESSL-4.4.1-1/lib
#BLAS_LINKER       := -lesslbg $(XLF_LIBS)
endif

# lapack libs
LAPACK_INC        := . # non-empty vars because XL doesn't like "-I -L "
LAPACK_LIB        := .
LAPACK_LINKER     := -llapack
ifeq ($(LIEU),NYX)
#LAPACK_LINKER     := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
LAPACK_LINKER       := -openmp -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#
endif
ifeq ($(LIEU),ANL)
LAPACK_LIB        := /soft/apps/current/LAPACK
# ALCF Netlib built with XLF
LAPACK_LINKER     := -llapack_bgp
# ESSL - linker acts left-to-right so this will pick up ESSL LAPACK symbols first then grab Netlib ones if the app needs some that are missing in ESSL
#LAPACK_LINKER     := -lesslbg -llapack_bgp
# This is the wrong way to do things if one is calling the handful of functions in ESSL that have the same name (symbol) but different call syntax than LAPACK
# _gesvd is a function that seems to require LAPACK be linked first
endif

# cuda libs
CUDA_INC          := $(CUDA_INSTALL_PATH)/include
CUDA_LIB          := $(CUDA_INSTALL_PATH)/lib64
CUDA_LINKER       :=-lcudart

# cuda libs
CUDASDK_INC       := $(CUDA_SDK_PATH)/common/inc
CUDASDK_LIB       := $(CUDA_SDK_PATH)/lib

# cutil
CUTIL_LINKER      :=-lcutil_x86_64

# cublas
CUBLAS_LINKER     :=-lcublas

# All together
LIBS              := $(LIBS) -L$(BLAS_LIB) -L$(LAPACK_LIB) -L$(CUDA_LIB) -L$(CUDASDK_LIB)
INCLUDES          += -I. -I$(BLAS_INC) -I$(LAPACK_INC)

MATH_LIB          := -lm
LDFLAGS           := $(LAPACK_LINKER) $(BLAS_LINKER) $(MATH_LIB)	
CFLAGS 	          := -03
CXXFLAGS          := -I. -O3 $(DEFINES)
ifeq ($(ARCH),USE_GPU)
INCLUDES          += -I. -I$(CUDA_INC) -I$(CUDASDK_INC)
LDFLAGS           := $(LDFLAGS) $(CUDA_LINKER) $(CUBLAS_LINKER)
CUFLAGS	          := --ptxas-options=-v $(DEFINES)
endif
ifeq ($(PREC),USE_DOUBLE)
ifeq ($(LIEU),HOME)
CUFLAGS	          := $(CUFLAGS) -arch sm_13
endif
ifeq ($(LIEU),NYX)
CUFLAGS	          := $(CUFLAGS) -arch=compute_20
endif
endif


.DELETE_ON_ERROR: *.o
.IGNORE: *.h

OBJECTS_DIR?=objects
DEPENDENCIES_DIR?=dependencies

# Path where objects and dependencies will be put
ODIR:=objects
DDIR:=dependencies


srcc:=$(subst ./,,$(basename $(foreach dir,$(directories),$(wildcard $(dir)/*.c) $(wildcard $(dir)/*.cc))))
ifneq ($(LIEU),ANL)
srcu:=$(subst ./,,$(basename $(foreach dir,$(directories),$(wildcard $(dir)/*.cu))))
endif
srch:=$(foreach dir,$(directories),$(wildcard $(dir)/*.h))
deps:=$(foreach file,$(srcc) $(srcu),$(DDIR)/$(file).d)
objs:=$(foreach file,$(srcc) $(srcu),$(ODIR)/$(file).o)


$(DDIR)/%.d : %.cc
	@echo "Dependencies of $*.cc"
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$@

$(DDIR)/%.d : %.c
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
	@$(CC)  $(INCLUDES) $(CXXFLAGS) $(CFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o $@|'>$@	

ifneq ($(LIEU),ANL)
$(DDIR)/%.d : %.cu
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(INCLUDES) $(CUFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@
endif
ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$@
endif
endif

all: $(BINNAME)


sinclude $(deps)
%.h :
	@echo "Dummy call for suppressed header $@"

$(ODIR)/%.o : %.c
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) $(CFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.c=.d)
	@echo "Compilation of $<"
	$(CC) $(INCLUDES) $(CXXFLAGS) $(CFLAGS) -c $< -o $@

$(ODIR)/%.o : %.cc
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cc=.d)
	@echo "Compilation of $<"
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

ifneq ($(LIEU),ANL)
$(ODIR)/%.o : %.cu
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(INCLUDES) $(CUFLAGS) $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
	@echo "Compilation of $<"
ifeq ($(ARCH),USE_GPU)
	$(NVCC) $(INCLUDES) $(CUFLAGS) -c $< -o $@
endif
ifeq ($(ARCH),USE_CPU)
	$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ -c $< -o $@
endif
endif

$(BINNAME) : $(objs)
	@echo "Link of $@"
	@mkdir -p $(@D)
	$(CXX) -o$@ $(objs) $(INCLUDES) $(LIBS) $(LDFLAGS) 

.PHONY : clean cleanall all install

clean :
	@echo "Clean of $(DDIR) and $(ODIR)"
	@-rm -rf $(DDIR) $(ODIR)

cleanall :
	@-rm -rv $(OBJECTS_DIR) $(DEPENDENCIES_DIR)

cleandep :
	@echo "Cleaning dependencies"
	@rm -rf $(DDIR)
