##
## \file Makefile
## \brief Main makefile for the code
## \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
## \license This project is released under the GNU Public License. See LICENSE.
## \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
## \ingroup makefiles
##
## \param[in] ARCH processor to use
## \param[in] PARA to use or not to use parallel computing
## \param[in] PREC double or single precision (only double works right now but single could be made to work easily)
## \param[in] PROB equation to use
## \param[in] GAMM gamma model to use
## \param[in] FLUX Riemann solver to use
## \param[in] LIEU Compile location
## \param[in] NFLD Total number of fields
## \param[in] MASS Number of mass fraction fields
## \param[in] PROF to use or not to use profiling
## \param[in] ERROR to output or not to output final solution error (L1, L2, Linfinity)
## \param[in] CONS to output or not to ouptut conservation errors
##
## Parameter      | Default values (or other valid values)
## :------------- | :-------------
## ARCH           | USE_CPU (or USE_GPU)
## PARA           | NO_MPI (or USE_MPI)
## PREC           | USE_DOUBLE (or USE_FLOAT)
## PROB           | MULTIFLUID (or SINGLEFLUID or PASSIVE)
## GAMM           | GAMNCON (or GAMCONS)
## FLUX           | ROE (or HLL or LLF)
## LIEU           | HOME (or TACC or ANL or FLUX)
## NFLD           | -
## MASS           | 0-10
## PROF           | OFF (or ON)
## ERRR           | NOERROR (or ERROR)
## CONS           | NOCONS (or CONS)
## @cond           

# default options
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
PROF=OFF
ERRR=NOERROR
CONS=NOCONS

# configuration
SRCDIR            := src
SRCDIRS           := $(shell find $(SRCDIR) -type d -name '*')
MKDIR             := custom_makefiles
BINDIR            := bin
TARGET            := dgexec


#================================================================================
#
# DEFAULTS AND INITIALIZATION
#
#================================================================================

# default compilers
NVCC              := nvcc
CXX               := g++
DEPFLAGS          := -M

# default paths to cuda libs
CUDA_INSTALL_PATH := .
CUDA_SDK_PATH     := .

# default blas libs
BLAS_INC          := .
BLAS_LIB          := .
BLAS_LINKER       := -lblas

# default lapack libs
LAPACK_INC        := . 
LAPACK_LIB        := .
LAPACK_LINKER     := -llapack

# default GSL stuff
GSL_INC	          := .
GSL_LIB	          := .
GSL_LINKER        := -lgsl -lgslcblas

# default math lib
MATH_LINKER       := -lm

# default cuda flags
ifeq ($(PREC),USE_FLOAT)
CUFLAGS	          := -arch sm_13
endif
ifeq ($(PREC),USE_DOUBLE)
CUFLAGS	          := -arch sm_20
endif

# default cublas linker
CUBLAS_LINKER     :=-lcublas

# Other variable initializations
TRACE		  := 

#================================================================================
#
# INCLUDE MACHINE SPECIFIC FILES
#
#================================================================================
ifeq ($(LIEU),HOME)
include $(MKDIR)/home.mk
endif

ifeq ($(LIEU),FLUX)
include $(MKDIR)/flux.mk
endif

ifeq ($(LIEU),TACC)
include $(MKDIR)/tacc.mk
endif

ifeq ($(LIEU),ANL)
include $(MKDIR)/anl.mk
endif

#================================================================================
#
# COMMON MODIFIERS TO PREVIOUSLY DEFINED VARIABLES
#
#================================================================================

# if you are using MPI, modify the cxx compiler
ifeq ($(PARA),USE_MPI)
CXX               := mpicxx
endif

# Profiling with gprof
ifeq ($(PROF),ON)
CXX               := $(CXX) -pg
CC                := $(CC) -pg
endif

# defines
DEFINES           := -DHAVE_BLAS -DHAVE_LAPACK -D$(ARCH) -D$(LIEU) -D$(PREC) -D$(PROB) -D$(GAMM) -D$(FLUX) -D$(DIMS) -D$(PARA) -DN_Y=$(MASS) -DN_F=$(NFLD) -D$(ERRR) -D$(CONS)
ifeq ($(DIMS),ONED)
DEFINES           := $(DEFINES) -DD=1 -DDF=1
endif
ifeq ($(DIMS),TWOD)
DEFINES           := $(DEFINES) -DD=2 -DDF=1
endif
ifeq ($(DIMS),THREED)
DEFINES           := $(DEFINES) -DD=3 -DDF=2
endif

# cuda libs
CUDA_INC          := $(CUDA_INSTALL_PATH)/include
CUDA_LIB          := $(CUDA_INSTALL_PATH)/lib64
CUDA_LINKER       :=-lcudart

# cuda sdk libs
CUDASDK_INC       := $(CUDA_SDK_PATH)/common/inc
CUDASDK_LIB       := $(CUDA_SDK_PATH)/lib

#================================================================================
#
# PUT IT ALL TOGETHER
#
#================================================================================
LIBS              := $(LIBS) -L$(BLAS_LIB) -L$(LAPACK_LIB) -L$(GSL_LIB)
INCLUDES          := -I. -I$(BLAS_INC) -I$(LAPACK_INC) -I$(GSL_INC) $(foreach dir,$(SRCDIRS),-I$(dir))
LDFLAGS           := $(BLAS_LINKER) $(LAPACK_LINKER) $(MATH_LINKER) $(GSL_LINKER)
CFLAGS 	          := -O3
CXXFLAGS          := -I. -O3 $(DEFINES)
CUFLAGS	          := $(CUFLAGS) $(DEFINES) --ptxas-options=-v


# .cu compiler and libraries
ifeq ($(ARCH),USE_GPU)
LIBS              := $(LIBS) -L$(CUDA_LIB) -L$(CUDASDK_LIB)
INCLUDES          := $(INCLUDES) -I$(CUDA_INC) -I$(CUDASDK_INC)
LDFLAGS           := $(LDFLAGS) $(CUDA_LINKER) $(CUBLAS_LINKER) 
CUCOMPILER        := $(NVCC) $(CUFLAGS) 
endif
ifeq ($(ARCH),USE_CPU)
CUCOMPILER := $(CXX) $(CXXFLAGS) -x c++ 
endif

#================================================================================
#
# ACTUAL MAKE COMMANDS. TRY NOT TO MODIFY BELOW THIS SECTION
#
#================================================================================

.DELETE_ON_ERROR: *.o
.IGNORE: *.h

# Path where objects and dependencies will be put
ODIR?=objects
DDIR?=dependencies

srcc:=$(subst ./,,$(basename $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.c) $(wildcard $(dir)/*.cc))))
srcu:=$(subst ./,,$(basename $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cu))))
srch:=$(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.h))
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

$(DDIR)/%.d : %.cu
	@echo "Dependencies of $<"
	@mkdir -p $(@D)
	@$(CUCOMPILER) $(INCLUDES) $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o $@|'>$@

all: $(BINDIR)/$(TARGET)

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

$(ODIR)/%.o : %.cu
	@mkdir -p $(@D)
	@$(CUCOMPILER) $(INCLUDES) $(DEPFLAGS) $< |sed 's|$(<F:.cu=\.o)|$(ODIR)/$*.o |'>$(DDIR)/$(<:.cu=.d)
	@echo "Compilation of $<"
	@$(CUCOMPILER) $(INCLUDES) -c $< -o $@

$(BINDIR)/$(TARGET) : $(objs)
	@echo "Link of $@"
	@mkdir -p $(@D)
	$(CXX) -o$@ $(objs) $(INCLUDES) $(LIBS) $(LDFLAGS) $(TRACE)
	@echo "SUCCESS"

doc : $(objs)
	@echo "Generating documentation with Doxygen."
	doxygen Doxyfile

test : $(objs)
	@echo "Running the test suite"
	@./tests/runtests.py -f;

clean :
	@echo "Clean of $(DDIR), $(ODIR), and $(BINDIR)"
	@-rm -rf $(DDIR) $(ODIR) $(BINDIR)

cleanall :
	@-rm -rv $(ODIR) $(DDIR) $(BINDIR)

cleandep :
	@echo "Cleaning dependencies"
	@rm -rf $(DDIR)

.PHONY : clean cleanall all install doc test

## @endcond
