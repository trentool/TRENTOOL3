/*
 * gpuKnnBF_SoAmex.cpp
 *
 *  Created on: 28/11/2012
 *      Author: marmar
 */


#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int cudaFindKnn(int* indexes, float* distances, float* pointset, float* queryset, int kth, int thelier, int nchunks, int pointdim, int signallength);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    
	//Arguments: single(pointset),single(queryset),kth,thelier,nchunks
	if(nrhs<5){
		mexErrMsgTxt("Incorrect number of input arguments.");
	}
	if(nlhs>2){
		mexErrMsgTxt("Too many output arguments.");
	}
	for(unsigned int i=0; i<2; i++){
		if(!mxIsSingle(prhs[i])||mxIsComplex(prhs[i])){
			mexErrMsgTxt("First two input arguments must be real single vectors.");
		}
	}
	for(unsigned int i=2;i<5;i++){
		if(!mxIsNumeric(prhs[i])||mxIsComplex(prhs[i])){
			mexErrMsgTxt("Last three input arguments must be numeric real values.");
		}
	}

	double kth,thelier,nchunks;
	
	//input
		
	size_t mscalar = mxGetM(prhs[2]);
	size_t nscalar = mxGetN(prhs[2]);
	if( !(mscalar==1 && nscalar==1) ) {
	   	mexErrMsgTxt("Input argument kth must be a scalar.");
	}else{
		kth = mxGetScalar(prhs[2]);
	}
	mscalar = mxGetM(prhs[3]);
	nscalar = mxGetN(prhs[3]);
	if( !(mscalar==1 && nscalar==1) ) {
	   	mexErrMsgTxt("Input argument thelier must be a scalar.");
	}else{
		thelier = mxGetScalar(prhs[3]);
	}
	mscalar = mxGetM(prhs[4]);
	nscalar = mxGetN(prhs[4]);
	if( !(mscalar==1 && nscalar==1) ) {
	   	mexErrMsgTxt("Input argument nchunks must be a scalar.");
	}else{
		nchunks=mxGetScalar(prhs[4]);
	}
           
	//fprintf(stderr,"%d %d %d",(int)kth,(int)thelier,(int)nchunks);
	//double thelier = mxGetScalar(prhs[4]);
	//double nchunks = mxGetScalar(prhs[5]);

    	float* pointset = (float*)mxGetData(prhs[0]);
    	int dimMp = mxGetM(prhs[0]);
    	int dimNp = mxGetN(prhs[0]);
    	
	float* query = (float*)mxGetData(prhs[1]);
    	int dimMq = mxGetM(prhs[1]);  //signallength
    	int dimNq = mxGetN(prhs[1]);  //pointdim

	if(dimNp!=dimNq) mexErrMsgTxt("Pointset and queryset must have points of same dimension (#columns)");

	//output
	mwSize outdim[2];
	outdim[0]=dimMq;
	outdim[1]=(int)kth;
	plhs[0]=mxCreateNumericArray(2,outdim,mxINT32_CLASS,mxREAL);  //indexes. TODO: int32 or int64
	plhs[1]=mxCreateNumericArray(2,outdim,mxSINGLE_CLASS,mxREAL); //distances

	int* indexes = (int*)mxGetData(plhs[0]);
	float* distances = (float*)mxGetData(plhs[1]);

	int success=cudaFindKnn(indexes,distances,pointset,query,(int)kth,(int)thelier,(int)nchunks,dimNq,dimMq);
	if(!success) mexErrMsgTxt("Error detected in the GPU (check console)");

}
