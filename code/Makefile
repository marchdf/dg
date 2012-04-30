# configuration
directories       := . quadratures
BINNAME           := dg1d



# flags
CUDA_INSTALL_PATH := /usr/local/cuda
CUDA_SDK_PATH     := /home/marchdf/NVIDIA_GPU_Computing_SDK/C
INCLUDES += -I. -I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/common/inc
LIBS              := -L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_SDK_PATH)/lib -lcutil_x86_64
CFLAGS 	          := -03
CPPFLAGS          := -I. -DHAVE_BLAS -DHAVE_LAPACK -w -ansi
CUFLAGS	          := --ptxas-options=-v #-arch sm_13 
LDFLAGS           := -lrt -lm -lblas -lstdc++ -llapack -lcudart -lcublas 

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
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@

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
	@$(NVCC) $(CUFLAGS) $(INCLUDES) -M $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
	@echo "Compilation of $<"
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $< -o $@

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
