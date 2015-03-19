
CUDA_LIB_PATH='/usr/local/cuda/lib64';

mex('-v',['-L' '.'], '-lgpuKnnLibrary', ['-L' CUDA_LIB_PATH],'-lcudart', 'fnearneigh_gpu.cpp');
mex('-v',['-L' '.'], '-lgpuKnnLibrary', ['-L' CUDA_LIB_PATH],'-lcudart', 'range_search_all_gpu.cpp');
