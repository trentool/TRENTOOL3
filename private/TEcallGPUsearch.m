function [ncount] = TEcallGPUsearch(cfg,channelpair,ps_1,ps_2,ps_p2,ps_21,ps_12,ps_p21)

% TEcallGPUsearch(cfg,channelpair,ps_1,ps_2,ps_p2,ps_21,ps_12,ps_p21)
%
% TECALLGPUSEARCH: This function is wrapper for the GPU functions
% conducting the k-nearest neighbour and range search on trialwise embedded
% data. The function splits up the data according to the capacities (memory
% and max no. of inputs, see comments in the code).
%
% NOTES:
%   - this version of the code calls the multigpu-functions
%   - this version doesn't require a GPU ID, the ID is set by the srm
%   - the code is written for MATLAB Version 7.9 (R2009b) and later (~ is
%     used to ignore some function outputs)
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
%     - running srmd (small resource manager) if cfg.site = 'ffm_srmc'
%     - functions 'range_search_all_multigpu.mexa64' and
%	  - 'fnearneigh_multigpu.mexa64' are used for nearest neighbour and
%       range searches (Martinez-Zarzuela, 2012)
%     - this function is called by TEsurrogatestats_ensemble.m, the data
%       passed to this function has to be embedded first
%     - TEconsoleoutput
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
%       .numthreads = number of threads that can be run in a single block
%                     on the GPU
%       .maxgriddim = maximum dimension of the grid of blocks
%     (The last three parameters are used to determine the maximum size of
%     the array that can be passed to the GPU in one call. The memory size
%     limits the number of values that can be in one input array. The max.
%     number of threads and maximum grid size limits the dimension of the
%     input array.)
%   	.verbosity   = set the verbosity of console output (see 'help
%                      TEconsoleoutput', default: 'info_minor')
%   	.site        = can be either 'ffm', 'ffm_srmc', or 'other'
%                      (default = 'other'); set to 'ffm' if TRENTOOL is
%                      executed on the cluster in the MEG Lab Frankfurt using
%                      multiple GPUs with manual parallelization via
%                      specification of the GPUid; set to 'ffm_srmc' if TRENTOOL
%                      is executed on the cluster in the MEG Lab Frankfurt using
%                      multiple GPUs with the small resource manager (srmc),
%                      the manager automatically distributes jobs over GPUs,
%                      managing the available main memory and devices
%
%   additional hardware specifications of GPU device if .site='ffm':
%      .GPUid        = device to be used for GPU computation
%   additional hardware specifications of GPU device if .site='ffm' or 'other':
%      .GPUmemsize   = the memory of the GPU in MB (e.g. cfg.GPUmemsize =
%                      4200), this information is mandatory if cfg.site is
%                      not set to 'ffm'.
%   if .site='other' no multi-GPU support is provided, all jobs will be
%   executed on the device with the lowest ID.
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
% Zarzuela, Thomas Sattler
% Frankfurt 2015

% CHANGELOG
%
% 2014-01-13 PW: corrected the estimator, 'eps' is now subtracted from the
% distance vector before range searching
% 2014-01-24 PW: I made MI calculation optional (can be switched on/off
% using cfg.MIcalc
% 2014-04-14 PW: bugfix - ncount.nchunks didn't return the correct value
% before
% 2015-05-20 PW/TS included calls to the srm
% 2015-05-27 PW: function now calls 'srmc max' to determine max. memory on
% GPU device

%% define logging levels
LOG_INFO_MAJOR   = 1;
LOG_INFO_MINOR   = 2;
LOG_DEBUG_COARSE = 3;
LOG_DEBUG_FINE   = 4;

%% Get variables from cfg
% -------------------------------------------------------------------------

k_th      = cfg.kth_neighbors;
TheilerT  = cfg.TheilerT(channelpair);
chunk_ind = cfg.chunk_ind;
MIcalc    = cfg.MIcalc;

if isfield(cfg, 'pause_time')
    pause_time = cfg.pause_time;
else
    pause_time = 20;
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

switch cfg.site
    case 'ffm_srmc'
        gpu_memsize = TEsrmc(pause_time, 'maxmem', cfg.verbosity);
    case {'other', 'ffm'}
        gpu_memsize = cfg.GPUmemsize;
    otherwise
        error('TRENTOOL ERROR: unkown site (has to be ''ffm'', ''ffm_srmc'', or ''other'').')
end

max_dim  = cfg.numthreads*cfg.maxgriddim; % max. 1st dimension for data input


%% Calculate data partitioning for GPU calls
% -------------------------------------------------------------------------

chunkelem    = numel(ps_p21(chunk_ind==1,:));   % get no. values in biggest chunk (pointset p21)
chunkdim     = size(ps_p21(chunk_ind==1,:),1);  % get size of 1st dimension, i.e. no. search points
chunksize    = (chunkelem*4)/(1024*1024);       % get size in memory (in MB)


% decide how many chunks fit on the card in one run
max_chunksperrun = min([ ...
    floor(max_dim/chunkdim) ...
    floor(gpu_memsize/chunksize) ...
    ]);

% calculate the number of runs and chunks per run
n_chunks = chunk_ind(end);
if n_chunks > max_chunksperrun
    nrruns = ceil(n_chunks/max_chunksperrun);
    chunksperrun = ceil(n_chunks/nrruns);
else
    nrruns = 1;
    chunksperrun = n_chunks;
end
mem_run = 2.5*chunksperrun*chunksize;

msg = sprintf('Chunks in current data set: %.0f (%0.4f MB per chunk, total: %.2f MB)', n_chunks, chunksize, chunksize*n_chunks);
TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
msg = sprintf('Number of runs: %.0f (%0.4f MB per run)', nrruns, mem_run);
TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);

switch cfg.site
    case 'ffm_srmc'
        gpu_id = TEsrmc(pause_time, 'request', cfg.verbosity, sprintf('gpumem:%d', ceil(mem_run)));
    case 'ffm'
        gpu_id = cfg.GPUid - 1;
	% fprintf('\n\n site is Frankfurt - using GPU %d\n\n', gpu_id);
    case 'other'
        ; % do nothing
    otherwise
        error('TRENTOOL ERROR: unkown site (has to be ''ffm'', ''ffm_srmc'', or ''other'').')
end

%% allocate memory for neighbor counts

% all arrays are generated, if no MI calculation is requested ncount_12_1 and ncount_12_2 remain empty
ncount_p21_p2 = zeros(size(ps_p21,1), 1);
ncount_p21_2 = zeros(size(ps_p21,1), 1);
ncount_p21_21 = zeros(size(ps_p21,1), 1);

if MIcalc
	ncount_12_1 = zeros(size(ps_p21,1), 1);
	ncount_12_2 = zeros(size(ps_p21,1), 1);
end

%% Call GPU functions: knn search + range search
% -------------------------------------------------------------------------

cutpoint = chunksperrun;
nchunks  = chunksperrun;

for ii=1:nrruns

	msg = sprintf('channelpair %d (u = %d) - run %d of %d',channelpair,cfg.u_in_ms(channelpair),ii,nrruns);
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);

	%% get data for this run (i.e. call of GPU functions) by concatenating indiv. chunks

    if ii==nrruns  % take remaining data if this is the last run

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

        ind2 = size(ps_p21,1);
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

	%% TE

	% k nearest neighbors search (fixed mass)
	t = tic;
    switch cfg.site
        case {'ffm', 'ffm_srmc'}
            [~, distance_p21] = fnearneigh_multigpu(single(pointset_p21),single(pointset_p21),k_th,TheilerT,nchunks,gpu_id);
        case 'other'
            [~, distance_p21] = fnearneigh_gpu(single(pointset_p21),single(pointset_p21),k_th,TheilerT,nchunks);
    end
	clear index_p21;
	t = toc(t);
	msg = sprintf('knn search for TE - %.1f minutes',t/60);
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
	clear t;

	% n nearest neighbor range search (fixed radius)
    t = tic;
	radius_p21 = single(distance_p21(:,k_th));
	radius_p21 = radius_p21 - eps('single');
	clear distance_p21;

    switch cfg.site
        case {'ffm', 'ffm_srmc'}
            ncount_p21_p2(ind1:ind2) = range_search_all_multigpu(single(pointset_p2),single(pointset_p2),radius_p21,TheilerT,nchunks,gpu_id);
            ncount_p21_21(ind1:ind2) = range_search_all_multigpu(single(pointset_21),single(pointset_21),radius_p21,TheilerT,nchunks,gpu_id);
            ncount_p21_2(ind1:ind2)  = range_search_all_multigpu(single(pointset_2),single(pointset_2),radius_p21,TheilerT,nchunks,gpu_id);
        case 'other'
            ncount_p21_p2(ind1:ind2) = range_search_all_gpu(single(pointset_p2),single(pointset_p2),radius_p21,TheilerT,nchunks);
            ncount_p21_21(ind1:ind2) = range_search_all_gpu(single(pointset_21),single(pointset_21),radius_p21,TheilerT,nchunks);
            ncount_p21_2(ind1:ind2)  = range_search_all_gpu(single(pointset_2),single(pointset_2),radius_p21,TheilerT,nchunks);
    end
    t = toc(t);

	msg = sprintf('range search for TE - %.1f minutes',t/60);
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
	clear t;


	%% MI estimation (optional)
	if MIcalc

		% k nearest neighbors search (fixed mass)
        t = tic;
        switch cfg.site
            case {'ffm', 'ffm_srmc'}
                [~, distance_12]   = fnearneigh_multigpu(single(pointset_12),single(pointset_12),k_th,TheilerT,nchunks,gpu_id);
            case 'other'
                [~, distance_12]   = fnearneigh_gpu(single(pointset_12),single(pointset_12),k_th,TheilerT,nchunks);
        end
        t = toc(t);
		msg = sprintf('knn search for MI - %.1f minutes',t/60);
        TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
		clear t;

		% n nearest neighbor range search (fixed radius)
		t = tic;
		radius_12  = single(distance_12(:,k_th));
		radius_12  = radius_12 - eps('single');
		clear distance_12 index_12;

        switch cfg.site
            case {'ffm', 'ffm_srmc'}
                ncount_12_1(ind1:ind2) = range_search_all_multigpu(single(pointset_1),single(pointset_1),radius_12,TheilerT,nchunks,gpu_id);
                ncount_12_2(ind1:ind2) = range_search_all_multigpu(single(pointset_2),single(pointset_2),radius_12,TheilerT,nchunks,gpu_id);
            case 'other'
                ncount_12_1(ind1:ind2) = range_search_all_gpu(single(pointset_1),single(pointset_1),radius_12,TheilerT,nchunks);
                ncount_12_2(ind1:ind2) = range_search_all_gpu(single(pointset_2),single(pointset_2),radius_12,TheilerT,nchunks);
        end
        t = toc(t);

		msg = sprintf('range search for MI - %.1f minutes',t/60);
        TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
		clear t;
	end

end

%fprintf('\nNeighbour count - ok \n');

switch cfg.site
    case 'ffm_srmc'
        resource = sprintf('gpumem:%d', ceil(mem_run));
        unit     = sprintf('/dev/nvidia%d', gpu_id);
        TEsrmc(pause_time, 'return', cfg.verbosity, resource, unit);
    case {'other', 'ffm'}
        ;
    otherwise
        error('TRENTOOL ERROR: unkown site (has to be ''ffm'' or ''other'').')
end

% build output structure
ncount.ncount_p21_p2 = ncount_p21_p2;
ncount.ncount_p21_21 = ncount_p21_21;
ncount.ncount_p21_2  = ncount_p21_2;
if MIcalc
    ncount.ncount_12_1   = ncount_12_1;
    ncount.ncount_12_2   = ncount_12_2;
end
ncount.nchunks	     = cfg.chunk_ind(end);

return;
