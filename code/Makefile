# defaults
ARCH=USE_CPU
PREC=USE_DOUBLE
PROB=MULTIFLUID
GAMM=GAMNCON

# configuration
directories       := . quadratures
BINNAME           := dg1d


# flags
CUDA_INSTALL_PATH := /usr/local/cuda
CUDA_SDK_PATH     := /home/marchdf/NVIDIA_GPU_Computing_SDK/C
#INCLUDES += -I. -I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/common/inc
LIBS              := 
LDFLAGS           := -lrt -lm -lblas -lstdc++ -llapack 
ifeq ($(ARCH),USE_GPU)
INCLUDES += -I. -I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/common/inc
LIBS              := $(LIBS) -L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_SDK_PATH)/lib -lcutil_x86_64
LDFLAGS           := $(LDFLAGS) -lcudart -lcublas 
endif
CFLAGS 	          := -03
CPPFLAGS          := -I. -DHAVE_BLAS -DHAVE_LAPACK -w -ansi -D$(ARCH) -D$(PREC) -D$(PROB) -D$(GAMM)
CUFLAGS	          := --ptxas-options=-v -D$(ARCH) -D$(PREC) -D$(PROB) -D$(GAMM)
ifeq ($(PREC),USE_DOUBLE)
CUFLAGS	          := $(CUFLAGS) -arch sm_13
endif


# compilers
NVCC              := nvcc

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
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$@

$(DDIR)/%.d : %.c
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
	@$(CC) $(CPPFLAGS) $(CFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o $@|'>$@

$(DDIR)/%.d : %.cu
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@
else ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$@
endif


all: $(BINNAME)


sinclude $(deps)
%.h :
	@echo "Dummy call for suppressed header $@"

$(ODIR)/%.o : %.c
	@mkdir -p $(@D)
	@$(CXX) $(CPPFLAGS) $(CFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.c=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.c=.d)
	@echo "Compilation of $<"
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(ODIR)/%.o : %.cc
	@mkdir -p $(@D)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cc=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cc=.d)
	@echo "Compilation of $<"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(ODIR)/%.o : %.cu
	@mkdir -p $(@D)
ifeq ($(ARCH),USE_GPU)
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
else ifeq ($(ARCH),USE_CPU)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -x c++ -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
endif
	@echo "Compilation of $<"
ifeq ($(ARCH),USE_GPU)
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $< -o $@
else ifeq ($(ARCH),USE_CPU)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -x c++ -c $< -o $@
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
