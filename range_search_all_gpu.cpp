/*
 *
 *  Created on: 28/11/2012
 *      Author: marmar
 */


#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int cudaFindRSAll(int* h_bf_npointsrange, float* h_pointset, float* h_query, float* vecradius, int thelier, int nchunks, int pointdim, int signallength);

/*
 * Range search with varying radius for every point in queryset
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    
	//Arguments: single(pointset),single(queryset),single(vecradius),thelier,nchunks
	//vecradius is a vector of length queryset
	if(nrhs<5){
		mexErrMsgTxt("Incorrect number of input arguments.");
	}
	if(nlhs>1){
		mexErrMsgTxt("Too many output arguments.");
	}
	for(unsigned int i=0; i<3; i++){
		if(!mxIsSingle(prhs[i])||mxIsComplex(prhs[i])){
			mexErrMsgTxt("First three input arguments must be real single vectors.");
		}
	}
	for(unsigned int i=3;i<5;i++){
		if(!mxIsNumeric(prhs[i])||mxIsComplex(prhs[i])){
			mexErrMsgTxt("Last two input arguments must be numeric real values.");
		}
	}

	double radius,thelier,nchunks;
	
	//input
	size_t mscalar = mxGetM(prhs[3]);
	size_t nscalar = mxGetN(prhs[3]);
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
           
    	float* pointset = (float*)mxGetData(prhs[0]);
    	int dimMp = mxGetM(prhs[0]);
    	int dimNp = mxGetN(prhs[0]);
    	
    	float* query = (float*)mxGetData(prhs[1]);
    	int dimMq = mxGetM(prhs[1]);  //signallength
    	int dimNq = mxGetN(prhs[1]);  //pointdim
        
        float *vecradius = (float*)mxGetData(prhs[2]);
        int dimMr = mxGetM(prhs[2]);
        int dimNr = mxGetN(prhs[2]);

	if(dimNp!=dimNq) mexErrMsgTxt("Pointset and queryset must have points of same dimension");
	if(dimMp!=dimMq) mexErrMsgTxt("Pointset and queryset matrices must have same dimension");
    	if(dimMr!=dimMp) mexErrMsgTxt("Bad dimension in vector of distances");
	
	//output
	mwSize outdim[2];
	outdim[0]=dimMq;
	outdim[1]=1;
	plhs[0]=mxCreateNumericArray(2,outdim,mxINT32_CLASS,mxREAL);  //indexes. TODO: int32 or int64

	int* npointsrange = (int*)mxGetData(plhs[0]);
    
	//fprintf(stderr,"\nPointer: %p",npointsrange);
	int success=1;
	success=cudaFindRSAll((int*)npointsrange,(float*)pointset,(float*)query,(float*)vecradius,(int)thelier,(int)nchunks,dimNq,dimMq);
	if(!success) mexErrMsgTxt("Error detected in the GPU (check console)");


}
