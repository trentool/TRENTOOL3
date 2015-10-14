function [TEresult]=transferentropy(cfg, data, varargin)

% TRANSFERENTROPY computes the transfer entropy (TE) among given pairs of
% channels for a sequence of trials over periods of time
%
% This function is called by the function TEsurrogatestats,
% TEconditionstatssingle,TEconditionsstatsgroup, and TEgroupstats
%
% !! varargin is only used for the shift test when called from the function
% TEsurrogatestats, TEconditionstatssingle or TEgroup_calculate !!
%
%
% * REFERENCE INFORMATION
%     - The concept of TE appears in Schreiber's article,
%       "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464
%       (2000).
%     - For the estimation of probability densities needed for the TE
%       computation, the function implements the Kraskov-Stoegbauer-
%       Grassberger estimator described in Kraskov et al. "Estimating
%       mutual information", Phys. Rev. E 69 (6) 066138, (2004).
%
% * DEPENDENCIES
%     - Package TSTOOL is used at nearest neighbors searches
%       required for the KSG estimator. (Gnu Public License)
%       http://www.dpi.physik.uni-goettingen.de/tstool/
%     - The following Matlab toolboxes:
%           - signal processing toolbox
%           - statistic toolbox
%     - The functions
%           - TEactdetect
%           - TEchannelselect 
%           - TEtrialselect
%           - TEvalues
%           - TECvalues
%           - TEconsoleoutput
%           - TEwaitbar
%
%
% * INPUT PARAPETERS
%
%   data            = Fieldtrip datastructure MUST contain:
%       .trials     = three dimensional data matrix
%       .time       = vector 1 x numtoi, the time included in the
%                     data
%       .label      = vector 1 x numlabel, labels of channelnames included
%                     in the data
%       .fsample    = value of sampling rate
%
%   The configuration MUST contain:
%   cfg.sgncmb      = list of channelpairs
%                     cell array (Nx(source, target))
%   or
%   cfg.channel     = list of channels testing all by all
%
%   and
%
%   cfg.toi         = the time range of interest (vector 1 x numtoi) in ms
%                     e.g. (time_from, time_to)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   WARNING:
%   The span of time needed for embedding is: (max(dim)-1)*max(tau)
%   The prediction time starts after this embedding time. Hence the span of
%   time defined in cfg.toi must be a good deal longer than the embedding
%   time, at least embedding time plus 150 samples or max(cfg.predicttime_u).
%
%       |<  embedding time  >|< prediction time u...
%   ----|--------------------|-----------------------------------|-->
%       |<                       cfg.toi                        >|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   cfg.dim         = embedding dimension
%   cfg.predicttime_u = points ahead for the advance vector in ms
%
%   cfg.kth_neighbors = number of neighbors for fixed mass search (controls
%                     balance of bias/statistical errors) (default = 4)
%   cfg.TheilerT    = number of temporal neighbors excluded to avoid serial
%                     correlations (Theiler correction) (default = 'ACT')
%   cfg.trialselect = ACT threshholding of trials - 'ACT' ,'range' or 'no'
%                     (default = 'ACT')
%       if you chose 'ACT' or nothing add:
%       cfg.actthrvalue = max threshold for ACT; min threshold default = 0
%       cfg.minnrtrials = minimum Nr of trials with ACT < actthrest used to
%                     calculate transfer entropy
%       if you chose 'range' add:
%       cfg.trial_from  = Inferior limit for the trials to be considered
%       cfg.trial_to    = Superior limit for the trials to be considered
%
%       cfg.TEcalctype =   'V' : self-prediction of the target signal and cross-
%                     prediction are both made from states in source and 
%                     target that precede the target state to be predicted
%                     by cfg.predicttime_u.
%                     (this option is considered deprecated)
%                     'VW' : the self-prediction time for the target is tau
%                     and cross-predictions are made from source states
%                     that precede the target state to be predicted by
%                     cfg.predicttime_u.
%                     (to solve the problem of decreasing self-prediction
%                     accuracy for large prediction times)
%                     (this option is considered deprecated)
%                     'VW_ds' : the self-prediction time for the target is tau
%                     and cross-predictions are made from source states
%                     that precede the target state to be predicted by
%                     1 sample. This estimator respects certain conditions
%                     necessary for d-separation properties  in the causal
%                     graph of source and target that are necessary for
%                     proper delay reconstruction
%                     (default='VW_ds')
%   cfg.verbosity   = set the verbosity of console output (see 'help
%                     TEconsoleoutput', default: 'info_minor')
%
%   cfg.maxlag      = the range of lags for computing the auto correlation
%                     time: from -MAXLAG to MAXLAG (default = 1000)
%   cfg.surrogatetype = 'trialshuffling','trialreverse','blockresampling',
%                     'blockreverse1','blockreverse2', 'blockreverse3', or
%                     'swapneighbor'.
%
%                       original trial:     1 2 3 4 5 6
%                       trialshufling:      2 3 4 5 6 1
%                       trialreverse:       6 5 4 3 2 1
%                       blockresampling:    4 5 6 1 2 3
%                       blockreverse1:      3 2 1 6 5 4
%                       blockreverse2:      6 5 4 1 2 3
%                       blockreverse3:      4 5 6 3 2 1
%                       swapneighbors:      2 1 4 3 6 5
%
% * OUTPUT PARAMETERS
%
%  TEresult             = Output structure
%          .TEmat       = resultmatrix including transfer entropy(TE)
%                         values
%          .MImat       = resultmatrix including mutual information (MI)
%                         values
%          .dimord      = dimensions of TEmat and MImat
%          .cfg         = configuration file used to calculate TE
%          .trials      = trial numbers selcted from raw dataset
%          .act         = ACT matrix (channelcombi x 2 trial)
%          .sgncmb      = labels of channel combinations (source -> target)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 2.0 by Raul Vicente, Michael Wibral, and Michael Lindner
% Bonn 2011
%
% CHNAGELOG:
% 2011-12-28: ML changed the internal cells and matrices (datacell, datamat,
% ACT, trials, nrtrials) to a common structure (channelcombi x ??)
% Additionally options for the fMRI data analysis are added.
% 2011-08-09: MW - at each call to TEvalues.m I added an if else statement
% to swicth to the use of TCvalues if data.TEprepare.cfg.TEcalctype is
% 'VW' - i.e. the use of the new predictor is desired by the user.
% 2012-19-07: NP - added parallel computing option, Data4Embedding is not
% working with this option yet
% 2013-04-16: PW removed all old estimators, i.e. 'VW_ds' is the only
% estimator used for TE calculation
% 2013-04-16: PW added a missing index for passing the ACT value passed to
% TEC_dsvalues.m, variable 'exec_type' was renamed to 'Theiler_ACT' and now
% takes values '1' for ACT used for Theiler T or '0' for user defined 
% TheilerT
% 2015-04-17: PW added a check for scalar delays that lets users call 
% TEsurrogatestats without calling InteractionDelayReconstruction_calculate
% 2015-10-02: PW fixed a bug in calls to TEC_dsvalues (didn't use the 
% fixed Theiler correction if requested)

LOG_INFO_MINOR = 2;
LOG_DEBUG_COARSE = 2;

if ~isfield(cfg, 'verbosity'), cfg.verbosity = 'info_minor'; end

% check data
% -------------------------------------------------------------------------
TEconsoleoutput(cfg.verbosity, 'Checking data and config', LOG_DEBUG_COARSE);

[data] = ft_checkdata(data, 'datatype','raw');

% check the data structure
if ~isfield(data, 'trial'),
    error('\nTRENTOOL ERROR: data must be in ''.trial''-structure, see help!');
end;
if ~isfield(data, 'time'),
    error('\nTRENTOOL ERROR: data contains no ''.time''-structure, see help!');
end;
if ~isfield(data, 'label'),
    error('\nTRENTOOL ERROR: data contains no ''.label''-structure, see help!');
end;
if ~isfield(data, 'fsample'),
    error('\nTRENTOOL ERROR: data contains no ''.fsample''-structure, see help!');
end;
if size(data.time,1)>size(data.time,2)
    data.time=data.time';
end



% check configuration and set defaults
% -------------------------------------------------------------------------

% if not defined set defaults
if ~isfield(cfg, 'maxlag'),       cfg.maxlag = 1000;                 end;
if ~isfield(cfg, 'trialselect'),  cfg.trialselect = 'ACT';           end;
if ~isfield(cfg, 'shuffle'),      cfg.shuffle = 'no';                end;

if ~isfield(cfg, 'embedsource') || strcmp(cfg.embedsource, 'yes')
    noSourceEmb = false;
else
    noSourceEmb = true;
    TEconsoleoutput(cfg.verbosity, 'No embedding for source time series is used', LOG_INFO_MINOR);
end

% check if channel or channelcombinations are defined
if isfield(cfg, 'channel') && isfield(cfg, 'sgncmb') ,
    error('\nTRENTOOL ERROR: specify cfg.channel OR sfg.sgncmb, see help!');
elseif isfield(cfg, 'channel') && ~isfield(cfg, 'sgncmb') ,
    channelselect = 1;
elseif ~isfield(cfg, 'channel') && isfield(cfg, 'sgncmb') ,
    channelselect = 2;
end;

% check TheilerT input
if isfield(cfg, 'TheilerT'),
    if strcmp(cfg.TheilerT, 'ACT')
        %exc_type = 1;
        Theiler_ACT = 1;
    else
        %exc_type = 2;
        Theiler_ACT = 0;
    end
end;


% check the format of input vectors
if size(cfg.toi,1)>size(cfg.toi,2)
    cfg.toi=cfg.toi';
elseif size(cfg.dim,1)>size(cfg.dim,2)
    cfg.dim=cfg.dim';
elseif size(cfg.tau,1)>size(cfg.tau,2)
    cfg.tau=cfg.tau';
elseif size(cfg.predicttime_u,1)>size(cfg.predicttime_u,2)
    cfg.predicttime_u=cfg.predicttime_u';
elseif size(cfg.kth_neighbors,1)>size(cfg.kth_neighbors,2)
    cfg.kth_neighbors=cfg.kth_neighbors';
elseif size(cfg.TheilerT,1)>size(cfg.TheilerT,2)
    cfg.TheilerT=cfg.TheilerT';
end


if nargin == 2;
    shifttest = 0;
elseif nargin == 3 && strcmp(varargin{1}, 'shifttest') ;
    shifttest = 1;
end



% get values from cfg
% -------------------------------------------------------------------------
channelcombi=data.TEprepare.channelcombi ;
channelcombilabel=data.TEprepare.channelcombilabel ;
ACT=data.TEprepare.ACT;
trials=data.TEprepare.trials;
nrtrials=data.TEprepare.nrtrials;





% read data
% -------------------------------------------------------------------------
TEconsoleoutput(cfg.verbosity, 'Reading data', LOG_DEBUG_COARSE);

% read data in to a cell {channelcombi x 2} including data matrices
% (trial x time)
% for each channelcombination only data will be read into a cell array that
% are allowed by the criteria put on their ACT values

data4TE = cell(size(channelcombi,1),2);
for cc = 1:size(channelcombi,1)
    for pp = 1:2
        datamat = zeros(nrtrials(cc,pp),size(data.trial{1},2)); % MW: check if trial{1} should be trial{2} because the valid trials of the TARGET matter
        for ii = 1:nrtrials(cc,pp) % should be 1:nrtrials(cc,2) to take the data for the source at the valid trials of the TARGET
            datamat(ii,:)=data.trial{trials{cc,pp}(ii)}(channelcombi(cc,pp),:);
        end
        data4TE{cc,pp}=datamat;
        clear datamat;
    end
end

% read data for spatial embedding for fMRI Data as '3DAsEmbed'
if isfield(data, 'Data4Embedding')
    embcell = cell(size(channelcombi,1),2);
    for cc = 1:size(channelcombi,1)
        for pp = 1:2
            embdatamat = zeros(size(Data.Data4Embedding,2),size(Data.Data4Embedding{1,1},2),size(Data.Data4Embedding{1,1},3));
            for ii = 1:size(data.trial,2)
                embdatamat(ii,:,:)=data.Data4Embedding{ii}(channelcombi(cc,pp),:,:);
            end
            embcell{cc,pp}=embdatamat;
            clear embdatamat
        end
    end
end

% get time indices from TEprepare
timeindices = data.TEprepare.timeindices;


% convert k value from ms to sampling points
dimu = data.TEprepare.u_in_samples;
%dimu = round(cfg.predicttime_u/1000*data.fsample);

if isscalar(dimu)
    dimu = repmat(dimu, 1, size(channelcombi,1));
end


% create empty result structure
TEresult=[];





% calculate number of points inside the time series used for advance and delay
% embedding and check whether there is a sufficient number of them left
% -------------------------------------------------------------------------
if strcmp(cfg.trialselect, 'ACT')
    multiplyact= min( [max(max(ACT(:,2,:))) cfg.actthrvalue] );
    mindatapoints = (timeindices(2)+1-timeindices(1))-(max(cfg.dim)-1)*max(cfg.tau)*multiplyact-max(dimu);
else
    mindatapoints = (timeindices(2)+1-timeindices(1))-(max(cfg.dim)-1)*max(cfg.tau)*max(max(ACT(:,2,:)))-max(dimu);
end


if isfield(data, 'datatype')
    if strcmp(data.datatype, 'fMRI')
        minsamples = 50;
    else
        minsamples = 150;
    end
else
    minsamples = 150;
end

if mindatapoints <= minsamples
    disp(mindatapoints)
    error('\nTRENTOOL ERROR: not enough data points left after embedding ');
end


% Change to directory containing mex files for the nearest neighbors search
% PW 10/11/2014: this is obsolete, TSTool functions come now bundled with
% TRENTOOL
%[dir_mex] = TEarch(cfg);
%cd(dir_mex);
TEarch;



% Check calculation time by calculating one TE with pessimistic values for
% tau, dim, .... etc
% -------------------------------------------------------------------------
if strcmp(cfg.calctime, 'yes')
    if isfield(data, 'Data4Embedding') && strcmp(data.TEprepare.cfg.datatype, 'fMRI')
        %fprintf('\ncalculation of time is not implementend for fMRI data.\n')
    else
        TEconsoleoutput(cfg.verbosity, 'Checking calculation time of TE', LOG_INFO_MINOR);
        
        %get time of single TE calculation (pessimistic case)
        timetest = 0;
        for ii = 1:5
            tic
	    
	    currentTau = round(cfg.tau(1)*data.TEprepare.maxact);
	    if currentTau < 1
            currentTau = 1;		% this may become zero for small ACT values -> overwrite manually
        end
        
        if Theiler_ACT            
            theiler_corr = data.TEprepare.maxact;
        else
            theiler_corr = cfg.TheilerT;
        end
        
        [te, mi] = TEC_dsvalues( ...
            squeeze(data4TE{1,1}(1,timeindices(1):timeindices(2))),squeeze(data4TE{1,2}(1,timeindices(1):timeindices(2))), ...
            cfg.dim(1), ...
            currentTau, ...
            dimu(1), ...
            cfg.kth_neighbors, ...
            theiler_corr, ...
            cfg.extracond, ...
            1, ...
            1);
        
        
            timetest = timetest + toc;
        end
        timetest = timetest/5;
        
        % time * nr of loops
        timeappr = timetest*size(channelcombi,1)*mean(mean(nrtrials))*size(cfg.TheilerT,2);
        
        
        % if unshuffled and shuffeled data are calculated double the time
        if isfield(cfg, 'permtest')
            timeappr = timeappr*2;
        elseif isfield(cfg, 'permtest') && strcmp(cfg.shifttest, 'yes')
            timeappr = timeappr*3;
        elseif isfield(cfg, 'permtest') && isfield(cfg, 'NrSubjects')
            timeappr = timeappr*2*cfg.NrSubjects;
        elseif isfield(cfg, 'permtest') && strcmp(cfg.shifttest, 'yes') && isfield(cfg, 'NrSubjects')
            timeappr = timeappr*3*cfg.NrSubjects;
        end
        
        % time in minutes OR hours and days
        timehh = floor(timeappr/60^2);                      %hours
        if timehh<1
            timemm = floor(mod((timeappr/60), 60));         %minutes
            msg = sprintf('Calculation of TE takes appr. : %d minutes!', timemm);
            TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
        else
            msg = sprintf('Calculation of TE takes appr. : %d hours (%d days)!', timehh, timehh/24);
            TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
        end
    end
end

% Start calculation of TE
% -------------------------------------------------------------------------
TEconsoleoutput(cfg.verbosity, 'Calculating transfer entropy. Please wait ...', LOG_INFO_MINOR);
TEwaitbar('init', size(channelcombi,1), cfg.verbosity);

% create zeros result matrices

TEresult.TEmat=zeros(size(channelcombi,1),max(max(nrtrials)));
TEresult.MImat=zeros(size(channelcombi,1),max(max(nrtrials)));

%--------------------------------------------------------------------------
par_state = 0;
if isfield(cfg,'TEparallel') && ft_hastoolbox('DCT');
            if isfield(cfg.TEparallel,'parON')
                if strcmp(cfg.TEparallel.parON,'yes')
                    
                    par_state = 1;
                    parallelConfig = findResource('scheduler','configuration',defaultParallelConfig);
                    max_workers = parallelConfig.ClusterSize;
                    if ~isfield(cfg.TEparallel,'workers')
                        cfg.TEparallel.workers = max_workers;
                    end
                    
                    if  matlabpool('size')== 0                    
                        
                        if cfg.TEparallel.workers <= max_workers
                           matlabpool(cfg.TEparallel.workers)
                        else
                            matlabpool(max_workers)
                        end
                        
                    else if matlabpool('size') > cfg.TEparallel.workers
                            matlabpool close;
                            matlabpool(cfg.TEparallel.workers)
                        end
                    end
                end
            end
end
%--------------------------------------------------------------------------   
% loops for scanning channels with different parameter values for TE


if ~par_state  % non-parallel part
    for channelpair = 1:size(channelcombi,1)
        TEwaitbar('update', channelpair, cfg.verbosity);
        
        for t4t = 1:nrtrials(channelpair,2)
            
            % OLD trial1 = trials{channelpair,2}(t4t);
            % trial1 (and trial2) will be used as indices into the datamats
            % that are inside data4TE, these are indexed by the indices of the
            % trial numbers, not the trial numbers!
            
            trial1 =t4t;
            
            if strcmp(cfg.shuffle, 'no')
                
                trial2 = trial1;
                timespan = timeindices(1):timeindices(2);
                
            elseif strcmp(cfg.shuffle, 'yes')
                
                switch cfg.surrogatetype,
                    case 'trialshuffling'
                        if mod(t4t,nrtrials(channelpair,2)) == 0
                            trial2 = 1;
                        else
                            trial2 = t4t+1;
                        end
                        timespan = timeindices(1):timeindices(2);
                        
                    case 'blockresampling'
                        trial2 = t4t;
                        cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                        timespan =[cutpoint:timeindices(2),timeindices(1):cutpoint-1];
                        
                    case 'trialreverse'
                        trial2 = t4t;
                        timespan = timeindices(1):timeindices(2);
                        flipdim(timespan,2)
                        
                    case 'blockreverse1'
                        trial2 = t4t;
                        cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                        timespan =flipdim([cutpoint:timeindices(2),timeindices(1):cutpoint-1],2);
                        
                    case 'blockreverse2'
                        trial2 = t4t;
                        cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                        timespan =[flipdim(cutpoint:timeindices(2),2),timeindices(1):cutpoint-1];
                        
                    case 'blockreverse3'
                        trial2 = t4t;
                        cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                        timespan =[cutpoint:timeindices(2),flipdim(timeindices(1):cutpoint-1,2)];
                    case 'swapneighbors'
                        if mod(t4t,2)==0
                            trial2 = t4t-1;
                        else
                            trial2 = t4t+1;
                        end
                        timespan = timeindices(1):timeindices(2);
                    otherwise
                        error('TRENTOOL ERROR: Unknown surrogate type!')
                end
                
            end
            
            
            
            if shifttest == 1
                if strcmp(cfg.shifttype, 'onesample')
                    a=squeeze(data4TE{channelpair,1}(trial1,timeindices(1)+1:timeindices(2)));
                    b=squeeze(data4TE{channelpair,2}(trial2,timeindices(1):timeindices(2)-1));
                elseif strcmp(cfg.shifttype, 'predicttime')
                    a=squeeze(data4TE{channelpair,1}(trial1,timeindices(1)+dimu(channelpair):timeindices(2)));
                    b=squeeze(data4TE{channelpair,2}(trial2,timeindices(1):timeindices(2)-dimu(channelpair)));
                end
                
            else
                
                a=squeeze(data4TE{channelpair,1}(trial1,timespan));
                b=squeeze(data4TE{channelpair,2}(trial2,timespan));
                if isfield(data, 'Data4Embedding') && strcmp(data.TEprepare.cfg.datatype, 'fMRI')
                    a_e=squeeze(embcell{channelpair,1}(trial1,:,timespan));
                    b_e=squeeze(embcell{channelpair,2}(trial2,:,timespan));
                end
            end
            t = tic;
            
            currentTau = round(cfg.tau(channelpair)*ACT(channelpair,2,trials{channelpair,2}(t4t)));
            if currentTau < 1
                currentTau = 1;		% this may become zero for small ACT values -> overwrite manually
            end
            
            if Theiler_ACT
                theiler_corr = ACT(channelpair,2,trials{channelpair,2}(t4t));
            else
                theiler_corr = cfg.TheilerT;
            end
                
            if noSourceEmb
                [te, mi] = TEC_dsvalues_noSourceEmb( ...
                    a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.kth_neighbors,theiler_corr,cfg.extracond,t4t,channelpair);
            else
                [te, mi] = TEC_dsvalues( ...
                    a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.kth_neighbors,theiler_corr,cfg.extracond,t4t,channelpair);
            end
            
            t = toc(t);
            TEresult.TEmat(channelpair,t4t)=te;
            TEresult.MImat(channelpair,t4t)=mi;
            
        end
        
    end
    
else % parallel part
    
    for channelpair = 1:size(channelcombi,1)
        
        TEwaitbar('update', channelpair, cfg.verbosity);
        
        aux_TEmat = zeros(1,nrtrials(channelpair,2));
        aux_MImat = zeros(1,nrtrials(channelpair,2));
        
        parfor t4t = 1:nrtrials(channelpair,2)
            
            % OLD trial1 = trials{channelpair,2}(t4t);
            % trial1 (and trial2) will be used as indices into the datamats
            % that are inside data4TE, these are indexed by the indices of the
            % trial numbers, not the trial numbers!
            
            trial1 =t4t;
            
            if strcmp(cfg.shuffle, 'no')
                trial2 = trial1;
                timespan = timeindices(1):timeindices(2);
            elseif strcmp(cfg.shuffle, 'yes')
                if strcmp(cfg.surrogatetype, 'trialshuffling')
                    if mod(t4t,nrtrials(channelpair,2)) == 0
                        trial2 = 1;
                    else
                        trial2 = t4t+1;
                    end
                    timespan = timeindices(1):timeindices(2);
                    
                elseif strcmp(cfg.surrogatetype, 'blockresampling')
                    trial2 = t4t;
                    cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                    timespan =[cutpoint:timeindices(2),timeindices(1):cutpoint-1];
                    
                elseif strcmp(cfg.surrogatetype, 'trialreverse')
                    trial2 = t4t;
                    timespan = timeindices(1):timeindices(2);
                    flipdim(timespan,2)
                    
                elseif strcmp(cfg.surrogatetype, 'blockreverse1')
                    trial2 = t4t;
                    cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                    timespan =flipdim([cutpoint:timeindices(2),timeindices(1):cutpoint-1],2);
                    
                elseif strcmp(cfg.surrogatetype, 'blockreverse2')
                    trial2 = t4t;
                    cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                    timespan =[flipdim(cutpoint:timeindices(2),2),timeindices(1):cutpoint-1];
                    
                elseif strcmp(cfg.surrogatetype, 'blockreverse3')
                    trial2 = t4t;
                    cutpoint = round( (timeindices(2)-timeindices(1)+1) * rand(1));
                    timespan =[cutpoint:timeindices(2),flipdim(timeindices(1):cutpoint-1,2)];
                elseif strcmp(cfg.surrogatetype, 'swapneighbors')
                    if mod(t4t,2)==0
                        trial2 = t4t-1;
                    else
                        trial2 = t4t+1;
                    end
                    timespan = timeindices(1):timeindices(2);
                end
                
            end
            
            
            
            if shifttest == 1
                if strcmp(cfg.shifttype, 'onesample')
                    a=squeeze(data4TE{channelpair,1}(trial1,timeindices(1)+1:timeindices(2)));
                    b=squeeze(data4TE{channelpair,2}(trial2,timeindices(1):timeindices(2)-1));
                elseif strcmp(cfg.shifttype, 'predicttime')
                    a=squeeze(data4TE{channelpair,1}(trial1,timeindices(1)+dimu(channelpair):timeindices(2)));
                    b=squeeze(data4TE{channelpair,2}(trial2,timeindices(1):timeindices(2)-dimu(channelpair)));
                end
                
            else
                
                a=squeeze(data4TE{channelpair,1}(trial1,timespan));
                b=squeeze(data4TE{channelpair,2}(trial2,timespan));
                
            end
            
            currentTau = round(cfg.tau(channelpair)*ACT(channelpair,2,trials{channelpair,2}(t4t)));
            if currentTau < 1
                currentTau = 1;		% this may become zero for small ACT values -> overwrite manually
            end
            
            if Theiler_ACT
                theiler_corr = ACT(channelpair,2,trials{channelpair,2}(t4t));
            else
                theiler_corr = cfg.TheilerT;
            end
                
            if noSourceEmb
                [te, mi] = TEC_dsvalues_noSourceEmb( ...
                    a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.kth_neighbors,theiler_corr,cfg.extracond,t4t,channelpair);
            else
                [te, mi] = TEC_dsvalues( ...
                    a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.kth_neighbors,theiler_corr,cfg.extracond,t4t,channelpair);
            end
            
            aux_TEmat(t4t) = te;
            aux_MImat(t4t) = mi;
            
        end
        
        TEresult.TEmat(channelpair,1:length(aux_TEmat))=aux_TEmat;
        TEresult.MImat(channelpair,1:length(aux_MImat))=aux_MImat;
        
    end
end


% results of unshuffled data
% -------------------------------------------------------------------------

if shifttest ==1
    TEresult.shifttest='yes';
end

if strcmp(cfg.shuffle, 'no')
    
    TEresult.act=ACT;
    TEresult.trials = trials; % save used trials in Result matrix
    TEresult.dimord = 'chanpair_trial';
    TEresult.cfg = cfg;
    if channelselect == 1
        TEresult.label=cfg.channel;
        TEresult.sgncmb=channelcombilabel;
    elseif channelselect == 2
        TEresult.sgncmb=cfg.sgncmb;
    end
    
    
end


%% Returning to the working directory
%cd(working_directory)

return;

