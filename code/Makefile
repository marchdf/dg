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
CXX               := gcc
ifeq ($(LIEU),NYX)
CXX               := icc
endif
ifeq ($(PARA),USE_MPI)
CXX               := mpicc
endif
ifeq ($(LIEU),ANL)
CXX               := mpixlcxx
endif



# defines
DEFINES           := -DHAVE_BLAS -DHAVE_LAPACK -D$(ARCH) -D$(PREC) -D$(PROB) -D$(GAMM) -D$(FLUX) -D$(DIMS) -D$(PARA) -DN_Y=$(MASS) -DN_F=$(NFLD)
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
BLAS_INC          := 
BLAS_LIB          := 
BLAS_LINKER       := -lblas
ifeq ($(LIEU),NYX)
#BLAS_LINKER       := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
BLAS_LINKER       := -openmp -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#
endif
ifeq ($(LIEU),ANL)
BLAS_LIB        := /soft/apps/current/blas-lapack-lib -L/soft/apps/current/BLAS
BLAS_LINKER     := -lblas_bgp
endif

# lapack libs
LAPACK_INC        := 
LAPACK_LIB        := 
LAPACK_LINKER     := -llapack
ifeq ($(LIEU),NYX)
#LAPACK_LINKER     := -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
LAPACK_LINKER       := -openmp -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm #you can control the number of threads: export MKL_NUM_THREADS=#
endif
ifeq ($(LIEU),ANL)
LAPACK_LIB        := /soft/apps/current/LAPACK
LAPACK_LINKER     := -llapack_bgp
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

# flags
LDFLAGS           := -lrt -lm -lstdc++ $(LAPACK_LINKER) $(BLAS_LINKER)
CFLAGS 	          := -03
CXXFLAGS          := -I. -w -O3 -ansi $(DEFINES)
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
srcu:=$(subst ./,,$(basename $(foreach dir,$(directories),$(wildcard $(dir)/*.cu))))
srch:=$(foreach dir,$(directories),$(wildcard $(dir)/*.h))
deps:=$(foreach file,$(srcc) $(srcu),$(DDIR)/$(file).d)
objs:=$(foreach file,$(srcc) $(srcu),$(ODIR)/$(file).o)

$(DDIR)/%.d : %.cc
	@echo "Dependencies of $*.cc"
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$@

$(DDIR)/%.d : %.c
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
	@$(CC)  $(INCLUDES) $(CXXFLAGS) $(CFLAGS) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o $@|'>$@

$(DDIR)/%.d : %.cu
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(INCLUDES) $(CUFLAGS) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@
endif
ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$@
endif


all: $(BINNAME)


sinclude $(deps)
%.h :
	@echo "Dummy call for suppressed header $@"

$(ODIR)/%.o : %.c
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) $(CFLAGS) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.c=.d)
	@echo "Compilation of $<"
	$(CC) $(INCLUDES) $(CXXFLAGS) $(CFLAGS) -c $< -o $@

$(ODIR)/%.o : %.cc
	@mkdir -p $(@D)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cc=.d)
	@echo "Compilation of $<"
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

$(ODIR)/%.o : %.cu
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(INCLUDES) $(CUFLAGS) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
	@echo "Compilation of $<"
ifeq ($(ARCH),USE_GPU)
	$(NVCC) $(INCLUDES) $(CUFLAGS) -c $< -o $@
endif
ifeq ($(ARCH),USE_CPU)
	$(CXX) $(INCLUDES) $(CXXFLAGS) -x c++ -c $< -o $@
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
