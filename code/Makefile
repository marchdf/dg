# defaults
ARCH=USE_CPU
PREC=USE_DOUBLE
DIMS=ONED
PROB=MULTIFLUID
GAMM=GAMNCON
FLUX=ROE

# configuration
directories       := . quadratures
BINNAME           := dg1d

# paths to libs
CUDA_INSTALL_PATH := /usr/local/cuda
CUDA_SDK_PATH     := /home/marchdf/NVIDIA_GPU_Computing_SDK/C

# compilers
NVCC              := nvcc
CXX               := g++

# defines
DEFINES           := -DHAVE_BLAS -DHAVE_LAPACK -D$(ARCH) -D$(PREC) -D$(PROB) -D$(GAMM) -D$(FLUX) -D$(DIMS)

# libs
LIBS              := 

# blas libs
BLAS_INC          := 
BLAS_LIB          := 
BLAS_LINKER       := -lblas

# lapack libs
LAPACK_INC        := 
LAPACK_LIB        := 
LAPACK_LINKER     := -llapack

# cuda libs
CUDA_INC          := $(CUDA_INSTALL_PATH)/include
CUDA_LIB          := $(CUDA_INSTALL_PATH)/lib64
CUDA_LINKER       := -lcudart

# cuda libs
CUDASDK_INC       := $(CUDA_SDK_PATH)/common/inc
CUDASDK_LIB       := $(CUDA_SDK_PATH)/lib

# cutil
CUTIL_LINKER      := -lcutil_x86_64

# cublas
CUBLAS_LINKER     := -lcublas

# All together
LIBS              := $(LIBS) -L$(BLAS_LIB) -L$(LAPACK_LIB) -L$(CUDA_LIB) -L$(CUDASDK_LIB)

# flags
LDFLAGS           := -lrt -lm -lstdc++ $(LAPACK_LINKER) $(BLAS_LINKER)
CFLAGS 	          := -03
CXXFLAGS          := -I. -w -ansi $(DEFINES)
ifeq ($(ARCH),USE_GPU)
INCLUDES          += -I. -I$(CUDA_INC) -I$(CUDASDK_INC)
LDFLAGS           := $(LDFLAGS) $(CUDA_LINKER) $(CUBLAS_LINKER)
CUFLAGS	          := --ptxas-options=-v $(DEFINES)
endif
ifeq ($(PREC),USE_DOUBLE)
CUFLAGS	          := $(CUFLAGS) -arch sm_13
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
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$@

$(DDIR)/%.d : %.c
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
	@$(CC) $(CXXFLAGS) $(CFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o $@|'>$@

$(DDIR)/%.d : %.cu
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@
else ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$@
endif


all: $(BINNAME)


sinclude $(deps)
%.h :
	@echo "Dummy call for suppressed header $@"

$(ODIR)/%.o : %.c
	@mkdir -p $(@D)
	@$(CXX) $(CXXFLAGS) $(CFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.c=.d)
	@echo "Compilation of $<"
	$(CC) $(CXXFLAGS) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(ODIR)/%.o : %.cc
	@mkdir -p $(@D)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cc=.d)
	@echo "Compilation of $<"
	$(CXX) $(CXXFLAGS)  $(INCLUDES) -c $< -o $@

$(ODIR)/%.o : %.cu
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
else ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
	@echo "Compilation of $<"
ifeq ($(ARCH),USE_GPU)
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $< -o $@
else ifeq ($(ARCH),USE_CPU)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -x c++ -c $< -o $@
endif

$(BINNAME) : $(objs)
	@echo "Link of $@"
	@mkdir -p $(@D)
	$(CXX) -o$@ $(objs) $(LDFLAGS) $(INCLUDES) $(LIBS)

.PHONY : clean cleanall all install

clean :
	@echo "Clean of $(DDIR) and $(ODIR)"
	@-rm -rf $(DDIR) $(ODIR)

cleanall :
	@-rm -rv $(OBJECTS_DIR) $(DEPENDENCIES_DIR)

cleandep :
	@echo "Cleaning dependencies"
	@rm -rf $(DDIR)
