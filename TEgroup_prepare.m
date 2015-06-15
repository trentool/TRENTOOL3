function TEgroup_prepare(cfg,fileCell)

% TEGROUP_PREPARE:  This function prepares data for group analysis by
% finding a common embedding dimension for all data sets in the group.
% TEgroup_prepare loops over all combinations of interaction delays u,
% ragwitz dimensions dim and ragwitz delays tau, provided by the user. The
% maximum embedding dimension found over all combinations of parameters
% will be used as embedding dimension for all files within the group
% analysis.
% TEgroup_prepare furthermore checks the data for equality with respect to
% no. channels, channel combination and time of interest (toi). If any of
% these parameters are not equal over subjects, the function will give an
% error message.
%
% Note, other than most TRENTOOL functions, TEgroup_prepare takes file
% names (not paths!) instead of data structures as inputs ('fileCell'),
% see below for detailed information regarding the function's input and
% output. Results are appended to the provided files (field
% '.groupprepare') and saved.
%
%
% * REFERENCE INFORMATION
%
%   - transfer entropy
%     - The concept of TE appears in Schreiber's article,
%       "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464
%
% * DEPENDENCIES
%     - Package TSTOOL is used at nearest neighbors searches
%       required for the KSG estimator. (Gnu Public License)
%       http://www.dpi.physik.uni-goettingen.de/tstool/
%     - The following Matlab toolboxes:
%         - statistic toolbox
%     - The functions
%         - transferentropy
%         - TEcmc
%         - TEperm
%         - TEvalues
%         - TEwait
%
%
%
% * INPUT PARAMETERS
%
% fileCell = Cell including the names (not paths!) of all files for the
% 	      group analyses (note that the function has to be called from
% 	      within the folder, that contains the data)
%
% cfg: The configuration MUST contain (for further details on these
%      parameters see 'help TEprepare'):
%
%       .outputpath   = save path for prepared data sets
%       .ragdim       = range of embedding dimensions to scan vector
%                       (1xnumdim)
%       .ragtaurange  = vector (1x2) of min and max embedding delays (in
%                       multiples of the autocorrelation decay time)
%       .ragtausteps  = number of equidistant steps in ragtaurange
%                       (min 5) (default = 10)
%       .flagNei = 'Range' or 'Mass' type of neighbor search
%       .sizeNei = Radius or mass for the neighbor search according to
%                  flagNeighborhood
%       .repPred = repPred represents the number of points for which the
%                  prediction is performed (it has to be smaller than
%                  length(timeSeries)-(dimEmb-1)*tauEmb-u)
%
%	.predicttimemin_u    = minimum interaction delay u to be scanned
%	.predicttimemax_u    = maximum interaction delay u to be scanned
% 	.predicttimestepsize = step size of the resolution of u
%
%       .sgncmb      = list of channelpairs
%                      cell array (Nx(source, target))
%  or
%       .channel     = list of channels - testing will be done all-by-all
%
%       .Path2TSTOOL = Path to the folder including the TSTOOL package
%	               (this is not required if ensemble method is chosen)
%       .toi         = the time range of interest (vector 1 x 2) in seconds
%                      e.g. (time_from, time_to) (units: seconds)
%
% The data provided within the files referenced by 'fileCell' should have
% the following format:
%
%   data           = Fieldtrip raw data structure - it MUST contain:
%       .trial     = cell array (nr of channels x nr of samples) containing
%                    the data for each trial
%       .time      = cell (1xnr of samples) containing the time indices for
%                    each trial (in seconds)
%       .label     = cell (1xnr of channels), containing the labels
%                    (strings) of channels included in the data
%       .fsample   = value of sampling rate (in Hertz)
%
%
%
% OUTPUT PARAMETERS
% TEgroup_prepare does not return an output. The group dimension as well as
% the u range used for parameter estimation is appended to the file as a
% field '.groupprepared' and saved. TEprepare will check the input data for
% this field and overwirte any conflicting user input (see help TEprepare).
%
%
%
%
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 2.0 by Michael Lindner, Raul Vicente, Michael Wibral and Patricia
% Wollstadt
% Frankfurt 2010
%

% TEPREPARE this function checks the input data and parameter for
% completeness and correctness. Further, it optimizes the embedding
% parameters and adds a substructure to the data, which is nesseccary for
% the further functions.
% TEPREPARE has to be performed on all datasets first!!!
%
%
%
% * INPUT PARAMETERS
%
%
%   cfg: The configuration MUST contain:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   WARNING:
%   The span of time needed for embedding is: (max(dim)-1)*max(tau)
%   The prediction time starts after this embedding time. Hence the span of
%   time defined in cfg.toi must be a good deal longer than the embedding
%   time, leastwise a multiple of the prediction time (nrk).
%
%       |<  embedding time  >|< prediction time ...
%   ----|--------------------|-----------------------------------|-->
%       |<                       cfg.toi                        >|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  cfg.TEcalctype =   'VW_ds' : the self-prediction time for the target is tau
%                     and cross-predictions are made from source states
%                     that precede the target state to be predicted by
%                     1 sample. This estimator respects certain conditions
%                     necessary for d-separation properties  in the causal
%                     graph of source and target that are necessary for
%                     proper delay reconstruction
%                     (default='VW_ds')
%                     THE FOLLOWING ESTIMATORS ARE NO LONGER SUPPORTED:
%                     'V' : self-prediction of the target signal and cross-
%                     prediction are both made from states in source and
%                     target that precede the target state to be predicted
%                     by cfg.predicttime_u.
%                     'VW' : the self-prediction time for the target is tau
%                     and cross-predictions are made from source states
%                     that precede the target state to be predicted by
%                     cfg.predicttime_u.
%                     (to solve the problem of decreasing self-prediction
%                     accuracy for large prediction times)
%
%  cfg.predicttime_u = time ahead for the advance prediction (scalar, in
%                      ms) -- this is the effective delay of the acition
%                      from source to target that is assumed
%
%  cfg.ensemblemethod = 'yes' to use ensemble-method (Gomez-Herrero, 2013)
%                       (default = 'no')
%
%  cfg.outputpath     = path to a folder, where group prepared data is
%                       saved
%
%
%  if you choose 'ensemblemethod = 'yes'':
%       The size of the GPU memory has to be provided to
%       TEsurrogatestats_ensemble.m. Note, that using the Faes-Method
%       is mandatory when using the GPU. Also no shifttest can be computed
%       and the surrogatetype has to be set to 'trialperm' (see also
%       help TEsurrogatestats_ensemble).
%
%
%
%
%   cfg.kth_neighbors = number of neighbors for fixed mass search (controls
%                     balance of bias/statistical errors) (default = 4)
%   cfg.TheilerT    = number of temporal neighbors excluded to avoid serial
%                     correlations (Theiler correction) (default = ACT)
%
%
%  cfg.trialselect = ACT threshholding of trials - 'ACT' ,'range' or 'no'
%                    (default = 'ACT'; for fMRI default = 'no')
%      if you chose 'ACT' (or nothing):
%      cfg.actthrvalue = max threshold for ACT; min threshold
%      cfg.minnrtrials = minimum Nr of trials with ACT < actthrest used to
%                    calculate transfer entropy
%      if you chose 'range':
%      cfg.trial_from  = Inferior limit for the trials to be considered
%      cfg.trial_to    = Superior limit for the trials to be considered
%  cfg.maxlag      = the range of lags for computing the auto correlation
%                    time: from -MAXLAG to MAXLAG (default = 1000)
%
%
% in case of parallel computing
%
%  cfg.TEparallel.parON = 'yes' for paralell computing tool.
%
%  cfg.TEparallel.workers = number of workers for parallel computing, if it
%       is bigger than the default matlab configuration then the maximum
%       workers will equal the default matlab configuration
%
%

% CHANGELOG:
%
% 2013/10/10: PW changed function to the TRENTOOL 3 workflow
%
% 2014/06/23: PW: the function now only performs one TEprepare per data set, using
%		the biggest requested u - TEprepare doesn't differ for different values 
%		for u (the u becomes relevant only in TE estimation)


%% Remember the working directory
working_directory = pwd;

%% check input
% -------------------------------------------------------------------------

% check if outputpath was provided
if ~isfield(cfg,'outputpath');
    error('TRENTOOL: Please provide an output path for prepared data!');
elseif ~isdir(cfg.outputpath)
    error('TRENTOOL: The output path you provided doesn''t seem to exist!');
end

% check if u vector was provided, construct u vector
if ~isfield(cfg,'predicttimemin_u') || ~isfield(cfg,'predicttimemax_u') || ~isfield(cfg,'predicttimestepsize')
    error('TRENTOOL: Please provide a minimum and maximum u value as well as the step size!')
else
    cfg.predicttime_u = cfg.predicttimemax_u;
end

% check if user requested ragwitz as the optimize method for embedding
% dimensions (cao not supported)
if ~strcmp(cfg.optimizemethod,'ragwitz')
    error('TRENTOOL: Group comparison only work for optimize method ''ragwitz''!')
end

max_dim_scan =  max(cfg.ragdim);
N = length(fileCell);
fprintf('You provided %.0f data sets\n', N);

%% create output structures
% -------------------------------------------------------------------------

dim_all      = nan(N,1);
nrtrials_all = nan(N,1);

%% loop over data: run TEprepare for each data set
% -------------------------------------------------------------------------

%dataprepAll = cell(length(predicttimevec_u),N);

for currentSubject = 1:N
    
    % load data
    varinfile = who('-file',fileCell{currentSubject});
    load(fileCell{currentSubject});
    
    % store loaded data within the variable 'data' and clear
    x = strcat('data =',varinfile{1},';');
    eval(x)
    if ~strcmp(varinfile{1},'data')
        y=strcat( ['clear ',varinfile{1} ]);
        eval(y)
    end
    clear x y varinfile
    
    
    fprintf('\tPreparing data for subject %d ...\n', currentSubject);        
    dataprepared      = TEprepare(cfg,data);
    
    % get relevant data from TEprepare
    dim_all(currentSubject)      = dataprepared.TEprepare.optdim;
    nrtrials_all(currentSubject) = min(min(dataprepared.TEprepare.nrtrials));
      
    % check if parameters/data are equal across subjects
    if currentSubject == 1; % if its the first subject remember all relevant parameters
        
        % time of interest
        toifrom      = dataprepared.TEprepare.timeindices(1);
        toito        = dataprepared.TEprepare.timeindices(2);
        
        % channels/channelcombis
        channelcombi = dataprepared.TEprepare.channelcombi;
        channels     = dataprepared.label;
                
    else                    % compare parameters of all subjects to first subject to check for equality
        if toifrom ~= dataprepared.TEprepare.timeindices(1) || toito   ~= dataprepared.TEprepare.timeindices(2);
            error('TRENTOOL: Time of interest differs between subjects!');
        elseif ~isequal(channelcombi,dataprepared.TEprepare.channelcombi);
            error('TRENTOOL: Channel combinations differ between subjects!');
        elseif ~isequal(ismember(channels,dataprepared.label),ismember(dataprepared.label,channels));
            error('TRENTOOL: Channels (channel label) differ between subjects!');
        end;
    end
    
    % break if a subject needs the max. possible embedding dimension
    if dim_all(currentSubject) == max_dim_scan;
        break;
    end
        
    clear dataprepared;    
end

%% evaluate parameters returned by TEprepare and append to group data
% -------------------------------------------------------------------------

% find maximum dimension over all u and subjects
max_dim = max(dim_all);
fprintf('The maximum embedding dimension was %.0f\n\n', max_dim);
fprintf('Appending group prepare information to data sets ...\n');

% create code to check whether data sets went through TEprepare_group
% within the same run
TEsetRandStream;

numpart = round(47+rand(1,10)*10);
charpart = char(round(65+rand(1,10)*26));
code = [numpart(1:5),charpart(1:5),numpart(6:10),charpart(6:10)];

% prepare output structure
groupprepare = [];
groupprepare.max_dim          = max_dim;
groupprepare.min_nrtrials     = min(nrtrials_all);
groupprepare.code             = code;
groupprepare.predicttimevec_u = cfg.predicttimemin_u:cfg.predicttimestepsize:cfg.predicttimemax_u;

% % for debugging
% save(fullfile(cfg.outputpath,'groupprepare.mat'),'groupprepare');

for currentSubject = 1:N
    
    fprintf('Saving data for %s\n',fileCell{currentSubject});
    
    % load data
    varinfile = who('-file',fileCell{currentSubject});
    load(fileCell{currentSubject});
    
    % append group output to data and save to desired location
    appendstring = strcat(varinfile{1}, '.groupprepare = groupprepare;');
    eval(appendstring)
    savestring = strcat('save(''',fullfile(cfg.outputpath,fileCell{currentSubject}),''',''',varinfile{1},''');');
    eval(savestring)
    clear x y varinfile
    
end

cd(working_directory);
