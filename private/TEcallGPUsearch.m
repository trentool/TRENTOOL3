function [ncount] = TEcallGPUsearch(cfg,channelpair,ps_1,ps_2,ps_p2,ps_21,ps_12,ps_p21)

% TEcallGPUsearch(cfg,channelpair,ps_1,ps_2,ps_p2,ps_21,ps_12,ps_p21)
%
% TECALLGPUSEARCH: This function is wrapper for the GPU functions conducting
% the k-nearest neighbour and range search on trialwise embedded data. The 
% function splits up the data according to the capacities (memory and max
% no. of inputs, see comments in the code).
% 
%
%
% * REFERENCE INFORMATION
%
%   - transfer entropy
%     - The concept of TE appears in Schreiber's article,
%       "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464
%       (2000).
%     - For the estimation of probability densities needed for the TE
%       computation, the function implements the Kraskov-Stoegbauer-
%       Grassberger estimator described in Kraskov et al. "Estimating
%       mutual information", Phys. Rev. E 69 (6) 066138, (2004).
%
%   - knn algorithm
%     - Merkwirth, Parlitz & Lauterborn (2000). Fast nearest-neighbor 
%       searching for nonlinear signal processing. Phys. Rev. E, 
%       American Physical Society, 62, 2089-2097.
%
%   - ensemble method
%     - Gomez-Herrero, Wu, Rutanen, Soriano, Pipa & Vicente (2010). 
%       Assessing coupling dynamics from an ensemble of time series. 
%       arXiv preprint arXiv:1008.0539.
%
%
% * DEPENDENCIES
%     - functions 'range_search_all_gpu.mexa64' and 
%	'fnearneigh_gpu.mexa64' are used for nearest neighbour and
%	range searches (Martinez-Zarzuela, 2012)
%     - this function is called by TEsurrogatestats_ensemble.m, the data
%	entered into this function has to be embedded first
%
%
%
% * INPUT PARAMETERS
%
%   cfg             = configuration structure - it MUST contain:
%       .TheilerT   = number of temporal neighbors excluded to avoid serial
%                     correlations (Theiler correction) (for default ACT,
%                     the maximum ACT for current signalcombination over 
%                     all trials is taken)
%       .k_th       = number of neighbors for fixed mass search (controls
%                     balance of bias/statistical errors) (default = 4)
%	    .MIcalc     = 1 if mutual information (MI) should be calculated 
%                     additionally to TE, 0 if not (faster)
%       .chunk_ind  = vector with indices encoding the chunk a single
%                     pointset belongs to
%       .GPUid      = if more than one GPU device are installed, this index
%                     can be used to call the desired device (has to be an
%                     integer between 1 and no. devices)
% 	    .GPUmemsize = size of the GPU's memory in MB
%       .numthreads = number of threads that can be run in a single block
%                     on the GPU
%       .maxgriddim = maximum dimension of the grid of blocks
%     (The last three parameters are used to determine the maximum size of
%     the array that can be passed to the GPU in one call. The memory size 
%     limits the number of values that can be in one input array. The max. 
%     number of threads and maximum grid size limits the dimension of the 
%     input array.)
%
%
%   channelpair  = current channelcombination (this is needed to read
%                 the current TheilerT from the vector cfg.TheilerT
%
%   pointsets   = concatenation of individually embedded trials, original
%		  and surrogate data are stacked in the first dimension: 
%		  indices for individual chunks (all pointsets for one
%		  channelcombination are stored in 'chunk_ind' ('1' 
%		  indicates original data)
%		  dimension: [nr.trialsXdatapointsXnumpermutations dim]
%
%       ps_1	= source time series embedded
%	ps_2	= target time series embedded
%	ps_p2	= perdiction point + target time series embedded
%	ps_21	= joint target and source time series embedded
%	ps_12	= joint source and target time series embedded
%	ps_p21	= perdiction point + joint target and source time series 
%		  embedded
%
%
% * OUTPUT PARAMETERS
%
%  ncount = structure with concatenated neighbour counts for search spaces
%	    relevant for TE calculation (see function TEcalc.m)
%	.ncount_p21_p2
%	.ncount_p21_p2
%	.ncount_p21_21 
%	.ncount_p21_2 
%	.ncount_12_1 
%	.ncount_12_2 
%	.nchunks       = total nr of chunks used in the neighbour searches
%			 (equal to numpermutation+1)
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
% 
% Version 1.0 by Michael Wibral, Patricia Wollstadt, Mario Martinez-
% Zarzuela
% Frankfurt 2013

% CHANGELOG
%
% 2014-01-13 PW: corrected the estimator, 'eps' is now subtracted from the
% distance vector before range searching
% 2014-01-24 PW: I made MI calculation optional (can be switched on/off
% using cfg.MIcalc
% 2014-04-14 PW: bugfix - ncount.nchunks didn't return the correct value
% before

%% Get variables from cfg
% -------------------------------------------------------------------------

k_th      = cfg.kth_neighbors;
TheilerT  = cfg.TheilerT(channelpair); 
chunk_ind = cfg.chunk_ind;
MIcalc    = cfg.MIcalc;

% get device ID, note: Mario's code uses zero-based indexing (hence -1)!
if isfield(cfg,'GPUid')
    gpuid = cfg.GPUid-1;
else
    gpuid = 0;
end


%% Set parameters for GPU device
% -------------------------------------------------------------------------

% data that can go onto the GPU is limited in size by two parameters:
%	(1) Memory in the GPU RAM (on a Tesla C2075 this is ~5259 MB)
%	(2) Number of threads in a block and maximum grid dimension (on a 
%       Tesla C2075 this is 512 threads, 65535 maximum grid dimension)
%
% Both parameters limit the maximum number of data, that can go onto the
% GPU in one call. Input size has to be smaller in memory than (1).
% Additionally the first dimension of the input array can not exceed the 
% number of threads times the grid size.
% (Note that the GPU RAM has to hold input, output and temporary files 
% generated during the execution of the GPU functions, this is taken into
% account when determining maximum size of the input.)

mem_free = cfg.GPUmemsize;                % max. memory in MB for data input
max_dim  = cfg.numthreads*cfg.maxgriddim; % max. 1st dimension for data input
fprintf('Free memory: %.0f MB\n', mem_free);
fprintf('Max. dimension: %.0f\n', max_dim);


%% Calculate data partitioning for GPU calls
% -------------------------------------------------------------------------

chunksize    = numel(ps_p21(chunk_ind==1,:));   % get no. entries in biggest chunk (pointset p21)
chunk_dim    = size(ps_p21(chunk_ind==1,:),1);  % get size of 1st dimension
chunksperrun = floor(max_dim/chunk_dim);        % decide how many of these chunks fit on the card in one run

% check whether the intended chunksperrun exceed the max. memory on the GPU
% if so, recalculate 'chunksperrun' based on the available memory

mem_run = ((chunksperrun*chunksize*4)/1024)/1024; % memory needed per run (in MB) with current chunks per run

if mem_run*2 >= mem_free;    
    mem_free = mem_free * 1024 * 1024;                  % convert free memory to Byte
    mem_free = mem_free/2;
    mem_free_n_values = mem_free/4;                     % calculate no. single precision values that fit into free memory
    fprintf('Max. no. values (single precision): %.0f \n', mem_free_n_values);
    chunksperrun = floor(mem_free_n_values/chunksize);  % overwrite 'chunksperrun' based on memory
end

% calculate no. runs given 'chunksperrun'
nrruns	     = ceil(chunk_ind(end)/chunksperrun);	

fprintf('Max. number of chunks per run: %.0f \n',chunksperrun);
fprintf('Chunks in current data set: %.0f \n',chunk_ind(end));
fprintf('\nNumber of runs: %.0f \n',nrruns);

%% Call GPU functions: knn search + range search
% -------------------------------------------------------------------------

% all arrays are generated, if no MI calculation is requested ncount_12_1 and ncount_12_2 remain empty
ncount_p21_p2 = [];
ncount_p21_21 = []; 
ncount_p21_2  = [];
ncount_12_1   = [];
ncount_12_2   = [];

cutpoint = chunksperrun;
nchunks  = chunksperrun;

for ii=1:nrruns

	fprintf('\nRun %.0f of %.0f\n',ii,nrruns);
	
	%% get data for this run (i.e. call of GPU functions) by concatenating indiv. chunks
	
	fprintf('\t get data ...\n');
	t = tic;
	if ii==nrruns
		
		% take remaining data if this is the last call        
        
		% get number and index of first chunk for this run
		chunk_start = cutpoint-chunksperrun+1;        
		ind1 = find(chunk_ind==chunk_start,1);
		
		% get data from concatenated datasets
		pointset_p21 = ps_p21(ind1:end,:);        	
		pointset_p2  = ps_p2(ind1:end,:);
		pointset_21  = ps_21(ind1:end,:);	
		pointset_2   = ps_2(ind1:end,:);
		if MIcalc
		    pointset_1   = ps_1(ind1:end,:);
		    pointset_12  = ps_12(ind1:end,:);
		end    
		
		nchunks = (chunk_ind(end)-(chunk_start))+1;		
    
	else
        
		% get number and indices of first and last chunk for this run
		chunk_start = cutpoint-chunksperrun+1;
		chunk_end   = cutpoint;
		ind1        = find(chunk_ind==chunk_start,1);
		ind2        = find(chunk_ind==chunk_end,1,'last');
		
		pointset_p21 = ps_p21(ind1:ind2,:);	
		pointset_p2  = ps_p2(ind1:ind2,:);
		pointset_21  = ps_21(ind1:ind2,:);	
		pointset_2   = ps_2(ind1:ind2,:);
		if MIcalc
		    pointset_1   = ps_1(ind1:ind2,:);
		    pointset_12  = ps_12(ind1:ind2,:);
		end   
		
		cutpoint = cutpoint+chunksperrun;

	end
	t = toc(t); 
	fprintf('\t \t getting the data for this run took %.1f minutes\n',t/60);
	clear t;
    
	%% TE
	% k nearest neighbors search (fixed mass)
	fprintf('\t k nearest neighbour search for TE...\n');
	t = tic;
	%[index_p21, distance_p21] = fnearneigh_multigpu(single(pointset_p21),single(pointset_p21),k_th,TheilerT,nchunks,gpuid);	
	[index_p21, distance_p21] = fnearneigh_gpu(single(pointset_p21),single(pointset_p21),k_th,TheilerT,nchunks);
	t = toc(t);
	fprintf('\t \t the GPU knn search took %.1f minutes\n',t/60);
	clear t;
    
	%% n nearest neighbor range search (fixed radius)
	
	fprintf('\t range search ...\n');
        t = tic;
	radius_p21 = single(distance_p21(:,k_th));
	radius_p21 = radius_p21 - eps('single');
	clear distance_p21 index_p21;
    	
	%ncount_p21_p2 = cat(1,ncount_p21_p2,range_search_all_multigpu(single(pointset_p2),single(pointset_p2),radius_p21,TheilerT,nchunks,gpuid));
	%ncount_p21_21 = cat(1,ncount_p21_21,range_search_all_multigpu(single(pointset_21),single(pointset_21),radius_p21,TheilerT,nchunks,gpuid));
	%ncount_p21_2  = cat(1,ncount_p21_2,range_search_all_multigpu(single(pointset_2),single(pointset_2),radius_p21,TheilerT,nchunks,gpuid));
	ncount_p21_p2 = cat(1,ncount_p21_p2,range_search_all_gpu(single(pointset_p2),single(pointset_p2),radius_p21,TheilerT,nchunks));
	ncount_p21_21 = cat(1,ncount_p21_21,range_search_all_gpu(single(pointset_21),single(pointset_21),radius_p21,TheilerT,nchunks));
	ncount_p21_2  = cat(1,ncount_p21_2,range_search_all_gpu(single(pointset_2),single(pointset_2),radius_p21,TheilerT,nchunks));
	t = toc(t); 
	
	fprintf('\t \t each GPU range search took %.1f minutes (total for three searches: %.1f minutes)\n',(t/3)/60,t/60);
	clear t;
	
	
	%% MI
	if MIcalc
	
		% k nearest neighbors search (fixed mass)
		fprintf('\t k nearest neighbour search for MI...\n');
		t = tic;		
		%[index_12, distance_12]   = fnearneigh_multigpu(single(pointset_12),single(pointset_12),k_th,TheilerT,nchunks,gpuid);
		[index_12, distance_12]   = fnearneigh_gpu(single(pointset_12),single(pointset_12),k_th,TheilerT,nchunks);
		t = toc(t);
		fprintf('\t \t the GPU knn search took %.1f minutes \n',t/60);
		clear t;	
		
		%% n nearest neighbor range search (fixed radius)	
		fprintf('\t range search ...\n');
		t = tic;
		radius_12  = single(distance_12(:,k_th));
		radius_12  = radius_12 - eps('single');
		clear distance_12 index_12;
		
		%ncount_12_1 = cat(1,ncount_12_1,range_search_all_multigpu(single(pointset_1),single(pointset_1),radius_12,TheilerT,nchunks,gpuid));
		%ncount_12_2 = cat(1,ncount_12_2,range_search_all_multigpu(single(pointset_2),single(pointset_2),radius_12,TheilerT,nchunks,gpuid));	
		ncount_12_1 = cat(1,ncount_12_1,range_search_all_gpu(single(pointset_1),single(pointset_1),radius_12,TheilerT,nchunks));
		ncount_12_2 = cat(1,ncount_12_2,range_search_all_gpu(single(pointset_2),single(pointset_2),radius_12,TheilerT,nchunks));	
		t = toc(t); 
		
		fprintf('\t \t each GPU range search took %.1f minutes (total for two searches: %.1f minutes)\n',(t/2)/60,t/60);
		clear t;
	end	
	
end

fprintf('\nNeighbour count - ok \n');

% build output structure
ncount.ncount_p21_p2 = ncount_p21_p2;
ncount.ncount_p21_21 = ncount_p21_21;
ncount.ncount_p21_2  = ncount_p21_2;
ncount.ncount_12_1   = ncount_12_1;
ncount.ncount_12_2   = ncount_12_2;
ncount.nchunks	     = cfg.chunk_ind(end);

return;
