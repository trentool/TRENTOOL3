
#MEX=/usr/local/matlab2007a/bin/mex
MEX=/usr/local/MATLAB/R2013b/bin/mex



# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# Location of the CUDA Toolkit binaries and libraries
CUDA_PATH       ?= /usr/local/cuda
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_LIB_PATH   ?= $(CUDA_PATH)/lib
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin

# Common binaries
NVCC            ?= $(CUDA_BIN_PATH)/nvcc
GCC             ?= g++

# Extra user flags
EXTRA_NVCCFLAGS ?=
EXTRA_LDFLAGS   ?=
EXTRA_CCFLAGS   ?=
MARIO_CACHE_FLA	?=

# CUDA code generation flags
GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30
#Mario
GENCODE_SM35	:= -gencode arch=compute_35,code=sm_35

ifeq ($(nocache),1)
	GENCODE_FLAGS	:=  $(GENCODE_SM10)
else
	GENCODE_FLAGS   :=  $(GENCODE_SM20) $(GENCODE_SM30) $(GENCODE_SM35)
endif


# OS-specific build flags
ifneq ($(DARWIN),) 
      LDFLAGS   := -Xlinker -rpath $(CUDA_LIB_PATH) -L$(CUDA_LIB_PATH) -lcudart
      CCFLAGS   := -arch $(OS_ARCH) 
else
  ifeq ($(OS_SIZE),32)
      LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart
      CCFLAGS   := -m32
  else
      LDFLAGS   := -L$(CUDA_LIB_PATH)64 -lcudart
      CCFLAGS   := -m64
  endif
endif

# OS-architecture specific flags
ifeq ($(OS_SIZE),32)
      NVCCFLAGS := -m32
else
      NVCCFLAGS := -m64
endif

# Debug build flags
ifeq ($(dbg),1)
      CCFLAGS   += -g 
      NVCCFLAGS += -g -G
      TARGET    := debug
      MEX += -g #mario debug mex
else
      TARGET    := release
endif

#ENABLE OR DISABLE L1 CACHE FLAGS
ifeq ($(nocache),1)
	#-Xptxas -dlcm=cg
	MARIO_CACHE_FLA = -Xptxas -dlcm=cg
endif

# Common includes and paths for CUDA
INCLUDES      := -I. -I.. -I../../common/inc -I$(CUDA_INC_PATH) -I../../../shared/inc -I/usr/local/cuda-5.0/samples/common/inc/ -I/usr/local/cuda-5.0/samples/shared/inc
LDFLAGS       += $(LIBPATH_OPENGL)

# Target rules
all:  mex

gpuKnnLibrary.o: gpuKnnLibrary.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) -Xcompiler -fPIC -c gpuKnnLibrary.cu

libraryknn: gpuKnnLibrary.o
	${AR} -r libgpuKnnLibrary.a gpuKnnLibrary.o


mex: libraryknn	
	${MEX} -L. -lgpuKnnLibrary -v fnearneigh_gpu.cpp -L$(CUDA_PATH)/lib64 -lcudart -lcusparse -lcublas
	${MEX} -L. -lgpuKnnLibrary -v range_search_all_gpu.cpp -L$(CUDA_PATH)/lib64 -lcudart -lcusparse -lcublas
	

clean:
	rm -f *.o ../bin/knn_mario 

install: all
	rm -f *.cpp *.cu *.o *.a compileformatlab.m
