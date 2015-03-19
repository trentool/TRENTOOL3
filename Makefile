
#MEX=/usr/local/matlab2007a/bin/mex
MEX=/usr/local/matlab/bin/mex
CUDA_PATH       ?= /usr/local/cuda
all: mex
mex: libgpuKnnLibrary.a	
	${MEX} -L. -lgpuKnnLibrary -v fnearneigh_gpu.cpp -L$(CUDA_PATH)/lib64 -lcudart 
	${MEX} -L. -lgpuKnnLibrary -v range_search_all_gpu.cpp -L$(CUDA_PATH)/lib64 -lcudart 
	
clean:
	rm -f *.mexa64  
