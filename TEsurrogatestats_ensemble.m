function TEpermtest=TEsurrogatestats_ensemble(cfg,data)

% TESURROGATESTATS_ENSEMBLE: This function calculates the transfer entropy 
% values and performs a test of transfer entropy from experimental data 
% sets against surrogate data.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!      The function TEprepare has to be run on the data first!      !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% You can call this function directly as follows:
%         TEsurrogatestats(cfg, data)
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
%   - permutation test
%       - Maris & Oostenveld (2007). Nonparametric statistical testing of
%         EEG- and MEG-data. J. of Neuroscience Methods, 164, 177-190.
%
%   - ensemble method
%       - Gomez-Herrero, Wu, Rutanen, Soriano, Pipa & Vicente (2010). 
%         Assessing coupling dynamics from an ensemble of time series. 
%         arXiv preprint arXiv:1008.0539.
%   - Faes method
%       - Faes, Nollo & Porta (2013). Compensated transfer entropy as a 
%         tool for reliably estimating information Transfer in 
%         physiological time series. Entropy, 15, 198-219.
%
%
% * DEPENDENCIES
%     - Functions 'range_search_all_multigpu.mexa64' and 
%	  'fnearneigh_multigpu.mexa64' are used for nearest neighbour and
%	  range searches (Wollstadt, 2014)
%     - The following Matlab toolboxes:
%         - signal processing toolbox
%         - statistic toolbox
%     - The functions
%         - TEprepare
%	  - TEembedding
%	  - TEcallGPUsearch
%	  - TEcalc
%	  - TEpvalue
%
% * INPUT PARAMETERS
%
%   data            = Fieldtrip raw data structure - it MUST contain:
%       .trials     = cell array (1xnr of trials) containing the data for
%                     each trial
%       .time       = cell (1xnr of trials) containing the time indices for
%                     each trial
%       .label      = cell (1xnr of channels), containing the labels of 
%                     channels included in the data
%       .fsample    = value of sampling rate (in Hertz)
%       .TEprepare  = structure added by TEprepare 
%
% AND
%
%   cfg: The configuration MUST contain:
%
%  cfg.optdimusage  = 'maxdim' to use maximum of optimal dimensions over
%                      all channels for all channels, or 'indivdim' to use
%                      the individual optimal dimension for each channel.
%                      In case of using ragwitz criterion also the optimal
%                      embedding delay tau per channelcombi is used.
%
%   cfg.dim         = Value(s) for embedding dimension. In case of using
%                     cfg.optdimusage = 'maxdim' this has to be a scalar
%                     value. In case of cfg.optdimusage = 'indivdim' this 
%                     has to be a vector of the size (channelcombi x 1). 
%                     If not specified, the optimal dimension(s) found in 
%                     TEprepare will be used, which is the recommended 
%                     option!
%   cfg.tau         = embedding delay in units of act (x*act). If not
%                     specified (recommended option), the tau is used as
%                     followed:
%                     Depending optimizemethod in TEprepare:
%                           'ragwitz' = optimal tau found via ragwitz 
%                                       critrion
%                           'cao'     = cfg.tau given by user in TEprepare
%                     If not specified, the optimal embedding delay found  
%                     in TEprepare will be used, which is the recommended 
%                     option!
%   cfg.alpha       = significance level for statisatical permutation test 
%                     and correction for multiple comparison 
%                     (default = 0.05)
%   cfg.MIcalc	    = tells TRENTOOL to also calculate mutual information 
%		      (MI) for the current data set if set to 1. If set to
%		      0, MI will not be calculated, this makes the 
%		      calculation faster and requires less memory 
%		      (default = 1).
%   cfg.site        = can be set to 'ffm' if TRENTOOL is executed on the
%                     cluster in the MEG Lab Frankfurt (default = 'other')
%
%  hardware specifications of your GPU device:
%   cfg.GPUmemsize  = the memory of the GPU in MB (e.g. cfg.GPUmemsize =
%                     4200), this information is mandatory if cfg.site is
%                     not set to 'ffm'.
%   cfg.numthreads  = max. number of threads that can be run in one block
%                     on the GPU (default = 512)
%   cfg.maxgriddim  = max. grid dimension of the GPU (default = 65535)
%   cfg.GPUid       = tells TRENTOOL which GPU device to use if multiple
%                     devices are installed, set to 1 if only one device is
%                     installed (default = 1)
%
%   cfg.surrogatetype = If ensemble method is chosen for TE calculation, 
%                       'trialshuffling' is the only possible option for  
%                       the creation of surrogate data. Parameter will be  
%                       set to 'trialshuffling' if no other option is 
%                       chosen and will throw an error if any different 
%                       option is provided.
%                       'trialshuffling' will lead to a permutation of  
%                       trials times the number of permutations provided 
%                       for the permutation test, i.e.:
%
%                       original trial:    1 2 3 4 5 6
%                       permtrial1:        2 4 1 3 6 5
%                       permtrial2:        6 1 5 2 3 4
%                       permtrial3:        ...
%  
%   cfg.extracond   = perform conditioning in tansfer entropy formula on
%                     additional variables. This option is mandatory for 
%                     the ensemble method. 
%                     Values:
%                     'Faes_Method' - include the future sample of the 
%                     source at the prediction time into the state vector
%                     of the past of of the target to condition on it.
%                     In principle, this removes any volume conduction effect.
%                     (Faes L. et al, Phys Rev E, 2011)
%                     'Battaglia_Method' - NOT implemented yet. This
%                     methods conditions on the "mean activity" of the
%                     system, i.e. a global system state.This will be done
%                     by creating a channel carrying that signal which will
%                     then be used as an addional entry for the past state
%                     of the source 
%
%   cfg.shifttest   = has to be set to 'no' or not defined if ensemble 
%                     method is chosen for TE calculation as no shift test 
%                     can be computed.
%                     Accordingly, parameters 'shifttesttype' and 
%                     'shifttype' will be ignored.
%
%   cfg.numpermutation = nr of permutations in permutation test
%                       (default = 500)
%   cfg.permstatstype  = no permutation test is conducted when using ensemble 
%                        method, this parameter will be ignored
%   cfg.tail        = when using the ensemble method this parameter is 
%		      ignored as a one-tailed test is conducted always
%                     permutation tests) (default in TEsurrogatestats= 1)
%   cfg.correctm    = correction method used for correction of the multiple
%                     comparison problem - False discovery rate 'FDR' or
%                     Bonferroni correction 'BONF' (default = 'FDR')
%   cfg.fileidout   = string for the first part of the output filename.
%   cfg.embedsource = string ('yes'/'no') that may be used to switch off the
%		      embedding for the source time series.
%
%
%
% * OUTPUT PARAMETERS
%
%
%  TEpermtest
%            .TEpermvalues  = matrix with size:
%                             (channelpair,value)
%                           The last dimension "value" includes:
%                           1 - p_values of the statistic within the
%                               distribution given by the permutations
%                           2 - 1 (0), if the statistics is significant at
%                               the prescribed alpha level (or not)
%                           3 - 1 (0), if the statistics is significant
%                               after correction for multiple comparisons
%                               (or not)
%                           4 - difference between empirical TE value and
%                               median of the surrogate data distribution
%                           5 - 0, volume conduction (is set to 0 as the 
%                               Faes Method is mandatory for using the 
%                               ensemble method, see Faes, 2013)
%
%            .dimord     = dimensions of TEpermvalues
%            .cfg        = configuration file used to calculate TE and
%                          permtest
%            .sgncmb     = labels of channel combinations (source ->
%                          target)
%            .numpermutation = number of permutations
%            .ACT        = structure including
%                .act    = ACT matrix (channelcombi x 2 x trial)
%                .label  = label of channels in ACT matrix
%            .nr2cmc     = number of tests to correct for multiple
%                          comparisons
%            .TEprepare  = results of the function TEprepare from the
%                          data
%
% AND 
%
%  TEresult             = Equivalent to the utput structure of the function tranferentropy
%          .TEmat       = resultvector including transfer entropy(TE)
%                         values (one per channelpair)
%          .MImat       = resultvector including mutual information (MI)
%                         values (one per channelpair)
%          .dimord      = 'channelpair'; the dimensions of TEmat 
%                         and MImat
%          .cfg         = configuration file used to calculate TE
%          .trials      = trial numbers selected from raw dataset
%          .act         = ACT matrix (channelcombi x 2 x trial)
%          .sgncmb      = labels of channel combinations (source -> target)
%          .TEprepare   = results of the function TEprepare from the
%                         data
%   the field .instantaneousmixing will never be added to data if ensemble method is chosen
%   since no shift test is conducted
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
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral, Patricia
% Wollstadt
% Frankfurt 2013
%

% CHANGELOG
% 
% 2014-24-13 PW: I made MI calculation optional (can be switched on/off
% using cfg.MIcalc
% 2014-04-14 PW: function now uses TEsetRandStream.m
% 2014-04-14 PW: it is now checked, if tau is zero after multiplication with
% the ACT and rounding -> if yes, set to 1 manually
% 2014-04-14 PW: changed surrogate type from trialperm to trialshuffling to
% beconsistent with the CPU method


%% Remember the working directory
working_directory1 = pwd;

%% check data
% -------------------------------------------------------------------------
fprintf('Check data and config');

% check if TEprepare was performed
if ~isfield(data, 'TEprepare'),
    fprintf('\n')
    error('TRENTOOL ERROR: The function TEprepare must be performed on the data, see help!');
end;

% check data using checkdata from Fieldtrip
[data] = ft_checkdata(data, 'datatype','raw');

% check the data structure
if ~isfield(data, 'trial'),
    fprintf('\n')
    error('TRENTOOL ERROR: data must be in ''.trial''-structure, see help!');
end;
if ~isfield(data, 'time'),
    fprintf('\n')
    error('TRENTOOL ERROR: data contains no ''.time''-structure, see help!');
end;
if ~isfield(data, 'label'),
    fprintf('\n')
    error('TRENTOOL ERROR: data contains no ''.label''-structure, see help!');
end;
if ~isfield(data, 'fsample'),
    fprintf('\n')
    error('TRENTOOL ERROR: data contains no ''.fsample''-structure, see help!');
end;
if size(data.time,1)>size(data.time,2)
    data.time=data.time';
end

%% compare new cfg and cfg from TEprepare if equal fields exist
% -------------------------------------------------------------------------

doublefields = 0;
cfgTEprepare = data.TEprepare.cfg;

cfgfields = fieldnames(cfgTEprepare);
cfgfields2 = fieldnames(cfg);

for ii = 1:size(cfgfields,1);
    for jj = 1:size(cfgfields2,1);
        if strcmp(cfgfields{ii},cfgfields2{jj})
            doublefields = doublefields + 1;
        end
    end
end

clear cfgTEprepare

if doublefields  > 0
    fprintf('\n')
    error('TRENTOOL ERROR: Illegal attempt to overwrite entry generated by or used for TEprepare! Change cfg or rerun TEprepare. (see help)')
end


% add structures and values of data.TEprepare.cfg to cfg
names1 = fieldnames(data.TEprepare.cfg);
nr1 = size(names1,1);
for ii = 1:nr1
    eval(strcat('cfg.',names1{ii},' = getfield(data.TEprepare.cfg, {1}, names1{ii});'))
end

cfg.u_in_ms = data.TEprepare.u_in_ms;

%% check general and GPU configuration and set defaults
% -------------------------------------------------------------------------

% if not defined set defaults
if ~isfield(cfg, 'alpha'),          cfg.alpha = 0.05;                end;
if ~isfield(cfg, 'correctm'),       cfg.correctm = 'FDR';            end;
if ~isfield(cfg, 'surrogatetype'),  cfg.surrogatetype = 'trialshuffling'; end;
if ~isfield(cfg, 'embedsource'),    cfg.embedsource = 'yes';         end;

if isfield(cfg, 'tail') && cfg.tail ~= 1
    cfg = rmfield(cfg, 'tail');
    warning('You provided "cfg.tail" with tail not equal one, this will be ignored, when ensemble method is chosen, see help. \n');
end;

if ~isfield(cfg, 'fileidout'),
    fprintf('\n')
    error('TRENTOOL ERROR: cfg.fileidout must be defined, see help!');
end;

% check if Frankfurt specific options are requested
if ~isfield(cfg, 'site')
    cfg.site = 'other';
elseif ~(strcmp(cfg.site, 'ffm') || strcmp(cfg.site, 'other'))
    error('TRENTOOL ERROR: Unkown site.')
end

% check if the user provided the GPU's memory
if ~strcmp(cfg.site, 'ffm') && ~isfield(cfg,'GPUmemsize');
    error('TRENTOOL ERROR: You requested ensemble method with GPU calculation. Please provide your GPUs memory in "cfg.GPUmemsize". See help!');
end;

% check if the user provided the GPU's no. threads and grid dimension
% this info limits the maximum size of the input array to the knn- and
% range-search
if ~isfield(cfg,'numthreads'); cfg.numthreads = 512; end;    
if ~isfield(cfg,'maxgriddim'); cfg.maxgriddim = 65535; end;

% check if the user provided faes method or explicitly set 'extracond'
% to 'no' -> 'Faes_Method' has to be used when using the ensemble method
% to avoid volume conduction
% faes method is needed as a shift test is not possible when using the
% ensemble method
if ~isfield(cfg,'extracond') || ~strcmp(cfg.extracond,'Faes_Method');
    error('TRENTOOL ERROR: You requested ensemble method with GPU calculation. Please provide "Faes_Method" in "cfg.extracond". See help!');
end;

if isfield(cfg, 'shifttest')
    if ~strcmp(cfg.shifttest,'no');
        error('TRENTOOL ERROR: You requested ensemble method with GPU calculation AND a shifttest. A shifttest is not possible when using the ensemble method, use the Faes method instead, see help.');
    end
end

% check if the correct permutation method for surrogate generation is
% chosen
if ~isfield (cfg,'surrogatetype') || ~strcmp(cfg.surrogatetype,'trialshuffling')
    error('TRENTOOL ERROR:  If ensemble method is used, "cfg.surrogatetype" has to be set to "trialshuffling", see help.');
end

% check optimizemethod
if ~isfield(cfg, 'optdimusage'),  
    fprintf('\n')
    error('TRENTOOL ERROR: cfg.optdimusage is not defined, see help!')
else
    if strcmp(cfg.optdimusage, 'maxdim') == 0 && strcmp(cfg.optdimusage, 'indivdim') == 0 
        fprintf('\n')
        error(['TRENTOOL ERROR: ',cfg.optdimusage,' is a wrong input for cfg.optdimusage , see help!'])
    end
end;

% check if MI calculation is requested
if ~isfield(cfg,'MIcalc')
    % by default, switch MI calculation on 
    fprintf('\n')
    fprintf('TRENTOOL will also calculate mutual information (MI).')
    MIcalc = 1;
    cfg.MIcalc = 1;
else
    if cfg.MIcalc == 1
        MIcalc = 1;
    else
        MIcalc = 0;
    end
end
    
% check if local TE calculation is requested
if ~isfield(cfg,'TELcalc')
    % by default, switch TEL calculation off (requires a lot of disc space)     
    cfg.TELcalc = 0;
else 
	if cfg.TELcalc
		warning('TRENTOOL: You''re calculating local TE, this may require a lot of disc space when saving files!');
	end
end

% check dim 
if ~isfield(cfg, 'dim')
    if strcmp(cfg.optdimusage, 'indivdim')
        cfg.dim = data.TEprepare.optdimmat;
%         cfg.optdimusage = cfg.optdimusage;
    else
        cfg.dim(1:size(data.TEprepare.optdimmat,1),1) = data.TEprepare.optdim;
%         cfg.optdimusage = cfg.optdimusage;
    end
else
    if strcmp(cfg.optdimusage, 'indivdim')
        if size(cfg.dim,1) ~= size(data.TEprepare.channelcombi,1)
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.dim has to be in that size: (channelcombi x 1), see help!')
        elseif size(cfg.dim,2)>1
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.dim has to be in that size: (channelcombi x 1), see help!')
        end
    else
        if size(cfg.dim,1)>1 && size(cfg.dim,2)>1
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.dim must include a scalar, see help!');
        end
        if cfg.dim < data.TEprepare.optdim
            fprintf('\n')
            fprintf('TRENTOOL WARNING: specified embedding dimension (cfg.dim) is smaller then the optimal dimension from TEprepare.')
        elseif cfg.dim > data.TEprepare.optdim
            fprintf('\n')
            fprintf('TRENTOOL WARNING: specified embedding dimension (cfg.dim) is bigger then the optimal dimension from TEprepare.')
        end
    end
end;

% check tau
if ~isfield(cfg, 'tau')
    if strcmp(data.TEprepare.cfg.optimizemethod, 'ragwitz') 
        if strcmp(cfg.optdimusage, 'indivdim')
            cfg.tau = data.TEprepare.opttaumat;
        else
            cfg.tau(1:size(data.TEprepare.channelcombi,1)) = data.TEprepare.opttau;
        end
    elseif strcmp(data.TEprepare.cfg.optimizemethod, 'cao') 
        cfg.tau(1:size(data.TEprepare.channelcombi,1)) = data.TEprepare.cfg.caotau;
    end
    
else
    if strcmp(cfg.optdimusage, 'indivdim') && strcmp(data.TEprepare.cfg.optimizemethod, 'ragwitz') 
        if size(cfg.tau,1) ~= size(data.TEprepare.channelcombi,1)
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.tau has to be in that size: (channelconmbi x 1), see help!')
        elseif size(cfg.tau,2)>1
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.tau has to be in that size: (channelconmbi x 1), see help!')
        end
    else
        if size(cfg.tau,1)>1 && size(cfg.tau,2)>1
            fprintf('\n')
            error('TRENTOOL ERROR: cfg.tau must include a scalar, see help!');
        end
    end
    
end

% select trials
trials                = data.TEprepare.trials;
nrtrials              = data.TEprepare.nrtrials;
cfg.permtest.trials   = trials;
cfg.permtest.nrtrials = nrtrials;
    
    
% check TE parameter
if isempty(cfg.predicttime_u), error('TRENTOOL ERROR: specify cfg.predicttime_u, see help!');  end;

if ~isfield(cfg, 'kth_neighbors'),  cfg.kth_neighbors = 4;  end;

if ~isfield(cfg, 'TheilerT'),       cfg.TheilerT = 'ACT';   end;
if strcmp(cfg.TheilerT, 'ACT');
    
    % get TheilerT if cfg.TheilerT == 'ACT' -> max ACT for each sgncombi over trials
    % ACT: [#signalcombis 2 nrtrials]
    cfg.TheilerT = max(squeeze(data.TEprepare.ACT(:,2,:)),[],2);
        
else
    if size(cfg.TheilerT,1)>1 || size(cfg.TheilerT,2)>1
        fprintf('\n')
        error('TRENTOOL ERROR: cfg.TheilerT must include a scalar, see help!');
    end
end


% check the format of input vectors
if size(cfg.toi,1)>size(cfg.toi,2)
    cfg.toi=cfg.toi';
elseif size(cfg.predicttime_u,1)>size(cfg.predicttime_u,2)
    cfg.predicttime_u=cfg.predicttime_u';
elseif size(cfg.kth_neighbors,1)>1 || size(cfg.kth_neighbors,2)>1
    fprintf('\n')
    error('TRENTOOL ERROR: cfg.dim must include a scalar, see help!');
end


fprintf(' - ok');

%% check nr of permutations 
% -------------------------------------------------------------------------
fprintf('\nChecking number of permutations');

%nr2cmc=size(data.TEprepare.channelcombilabel,1)*size(cfg.predicttime_u,2);
nr2cmc=size(data.TEprepare.channelcombilabel,1);

findDelay = 0;

if ~isfield(cfg, 'numpermutation'),
    %cfg.numpermutation = 190100; % for p<0.01 with a possible bonferroni correcetion of 100
    cfg.numpermutation = ceil(1/(cfg.alpha/nr2cmc));
    fprintf('TRENTOOL: You didn''t specify a number of permutations. It was set to %d (1/(alpha/no_channelcombis)).', cfg.numpermutation);
end

if strcmp(cfg.numpermutation, 'findDelay');
    cfg.numpermutation = 0;
    findDelay = 1;
else 

    if cfg.numpermutation < ceil(1/cfg.alpha)
    	fprintf('\n')
    	error('TRENTOOL ERROR: cfg.numpermutation too small (< 1/alpha)!');
    elseif cfg.numpermutation < ceil(1/(cfg.alpha/nr2cmc))
       fprintf('\n###############################################\n# WARNING: Nr of permutations not sufficient for correction for multiple comparisons! #\n#######################################################################################\n'); 
    elseif max(nrtrials)>31 && cfg.numpermutation > 2^31
        fprintf('\n')
        error('TRENTOOL ERROR: cfg.numpermutation too huge (> 2^31)!');
    elseif max(nrtrials)>31 && cfg.numpermutation > 2^min(min(nrtrials)) % nrtrials is now 2-D!
        fprintf('\n')
        error('TRENTOOL ERROR: cfg.numpermutation too huge (> 2^n_trials)!');
    end

end

fprintf(' - ok\n');

%% get channels, ACT and trials for TEperm and embedding
% ------------------------------------------------------------------------

channelcombi      = data.TEprepare.channelcombi ;
channelcombilabel = data.TEprepare.channelcombilabel ;
ACT               = data.TEprepare.ACT;
timeindices       = data.TEprepare.timeindices;
dimu              = data.TEprepare.u_in_samples;
numpermutation    = cfg.numpermutation;

% remember the TEprepare structure to later add it to the output
TEpreparestruct   = data.TEprepare;

%% read data
% -------------------------------------------------------------------------
fprintf('Read data');

% read data in to a cell {channelcombi x 2} including data matrices
% (trial x time)
% for each channelcombination only data will be read into a cell array that
% are allowed by the criteria put on their ACT values

data4embedding = cell(size(channelcombi,1),2);
for cc = 1:size(channelcombi,1)
    for pp = 1:2
        datamat = zeros(nrtrials(cc,pp),size(data.trial{1},2)); % MW: check if trial{1} should be trial{2} because the valid trials of the TARGET matter
        for ii = 1:nrtrials(cc,pp) % should be 1:nrtrials(cc,2) to take the data for the source at the valid trials of the TARGET
            datamat(ii,:)=data.trial{trials{cc,pp}(ii)}(channelcombi(cc,pp),:);
        end
        data4embedding{cc,pp}=datamat;
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

fprintf(' - ok\n');

%% check no. samples
% -------------------------------------------------------------------------

% calculate number of points inside the time series used for advance and delay
% embedding and check whether there is a sufficient number of them left

if strcmp(cfg.trialselect, 'ACT')
    multiplyact= min( [max(max(ACT(:,2,:))) cfg.actthrvalue] );
    mindatapoints = (timeindices(2)+1-timeindices(1))-(max(cfg.dim)-1)*max(cfg.tau)*multiplyact-max(dimu);
else
    mindatapoints = (timeindices(2)+1-timeindices(1))-(max(cfg.dim)-1)*max(cfg.tau)*max(max(ACT(:,2,:)))-max(dimu);
end

fprintf('Min. sample points left after embedding: %.0f ',mindatapoints);

if isfield(data, 'datatype')
    if strcmp(data.datatype, 'fMRI')
        minsamples = 50/min(min(nrtrials));
    else
        minsamples = 150/min(min(nrtrials));
    end
else    
    % determine minsamples over trials, as samples are later pooled over 
    % trials for neighbor searches
    minsamples = 150/min(min(nrtrials)); 
end

fprintf('(min. required: %.0f.)',ceil(minsamples));

% compare samples in analysis window against minimum feasible no. sample points
if mindatapoints <= minsamples    
    error('\nTRENTOOL ERROR: not enough data points left after embedding.');
end

fprintf(' - ok');

%% clear input data to save memory
% -------------------------------------------------------------------------

% data is no longer needed at this point
% this serves memory as embedding all datapoints at once over trials 
% is quite memory intensive
clear data;


%% loop over channelcombinations, embed original data and create surrogates
% -------------------------------------------------------------------------
timeindices = TEpreparestruct.timeindices; 

%fprintf('\nStart embedding original and surrogate data for individual channelpairs');

% prepare output structures
TEpermvalues = zeros(size(channelcombi,1),5);
TEmat        = zeros(size(channelcombi,1),1);
MImat        = zeros(size(channelcombi,1),1);
TEmat_sur    = zeros(size(channelcombi,1),numpermutation);
if cfg.TELcalc
    TELmat 	     = cell(1,size(channelcombi,1));
else
    TELmat = [];
end

TEsetRandStream;

for channelpair = 1:size(channelcombi,1)
	
    fprintf('\nEmbedding original data for channelpair %.0f of %.0f',channelpair,size(channelcombi,1));		
    %fprintf('\n\tEmbedding orginal data');
	
    % prepare data strucuters for embedding	
	
	pointsets_concat_2   = [];
	pointsets_concat_p2  = [];
	pointsets_concat_21  = [];	
	pointsets_concat_p21 = [];
	
	% if no MI calculation is requested, these remain empty and are passed to TEcallGPUsearch
	pointsets_concat_1   = [];
	pointsets_concat_12  = [];

	
	% embed original data per trial
	for t4t = 1:nrtrials(channelpair,2)
		
        timespan = timeindices(1):timeindices(2);
        trial1 = t4t; 
        trial2 = trial1;
		a=squeeze(data4embedding{channelpair,1}(trial1,timespan));
		b=squeeze(data4embedding{channelpair,2}(trial2,timespan));
		
		
		% TEembedding(a,b,cfg.dim(channelpair),round(cfg.tau(channelpair)*ACT(channelpair,trials{channelpair,2}(t4t))),dimu,cfg.extracond);
		
		currentTau = round(cfg.tau(channelpair)*ACT(channelpair,2,trials{channelpair,2}(t4t)));
		if currentTau < 1
			currentTau = 1;		% this may become zero for small ACT values -> overwrite manually
		end
	
		% check if source needs to be embedded too
		if strcmp(cfg.embedsource,'no')
			pointset = TEembedding_noSourceEmb(a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.extracond);
		else
			pointset = TEembedding(a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.extracond);
		end

		pointsets_concat_2 = cat(1,pointsets_concat_2,pointset.pointset_2);
		pointsets_concat_p2 = cat(1,pointsets_concat_p2,pointset.pointset_p2);
		pointsets_concat_21 = cat(1,pointsets_concat_21,pointset.pointset_21);		
		pointsets_concat_p21 = cat(1,pointsets_concat_p21,pointset.pointset_p21);
		
		if MIcalc
		    pointsets_concat_1 = cat(1,pointsets_concat_1,pointset.pointset_1);
		    pointsets_concat_12 = cat(1,pointsets_concat_12,pointset.pointset_12);
		end
		
		clear pointset;
    end
    
    % preallocate memory for surrogate data embedding
    % data is casted to single precision to save memory
    % 'chunksize' is the number of embedded points summed over all trials
    chunksize  = size(pointsets_concat_2,1);      
    pointsets_concat_2   = cat(1,single(pointsets_concat_2),zeros(chunksize*numpermutation,size(pointsets_concat_2,2)));
    pointsets_concat_p2  = cat(1,single(pointsets_concat_p2),zeros(chunksize*numpermutation,size(pointsets_concat_p2,2)));
    pointsets_concat_21  = cat(1,single(pointsets_concat_21),zeros(chunksize*numpermutation,size(pointsets_concat_21,2)));    
    pointsets_concat_p21 = cat(1,single(pointsets_concat_p21),zeros(chunksize*numpermutation,size(pointsets_concat_p21,2)));
    if MIcalc
        pointsets_concat_1   = cat(1,single(pointsets_concat_1),zeros(chunksize*numpermutation,size(pointsets_concat_1,2)));
        pointsets_concat_12  = cat(1,single(pointsets_concat_12),zeros(chunksize*numpermutation,size(pointsets_concat_12,2)));
    end


    % carry a indices vector to identify individual chunks later   
    chunk_ind = zeros(chunksize*(numpermutation+1),1);
    chunk_ind(1:chunksize) = 1;
    cutpoint = chunksize*2;
    
    fprintf('\t - ok\n');
	
	% shuffle orig data and embed shuffled data per trial  
    if numpermutation > 0
        ft_progress('init', 'text','Generating surrogate data sets ...')
    end
    
    for ii = 1:numpermutation
        ft_progress(ii/numpermutation, '\tdata set %d of %d', ii, numpermutation);
        
        % get permutation for trials in second/target channel
        channel2_shuffle = randperm(nrtrials(channelpair,2)); %\\ TODO check randstream
        
        % build auxiliary data structures to collect embedded surrogate data over trials       
        pointsets_concat_2_aux   = [];
        pointsets_concat_p2_aux  = [];
        pointsets_concat_21_aux  = [];        
        pointsets_concat_p21_aux = [];
	if MIcalc
	    pointsets_concat_1_aux   = [];
	    pointsets_concat_12_aux  = [];
	end
        
        % embed shuffled data per trial
        for t4t = 1:nrtrials(channelpair,2)
            
            % get trials
            trial1 = t4t;
            trial2 = channel2_shuffle(t4t);
            timespan = timeindices(1):timeindices(2);
            
            a=squeeze(data4embedding{channelpair,1}(trial1,timespan));
            b=squeeze(data4embedding{channelpair,2}(trial2,timespan));
            
            % TEembedding(a,b,cfg.dim(channelpair),round(cfg.tau(channelpair)*ACT(channelpair,trials{channelpair,2}(t4t))),dimu,cfg.extracond);
	    
	    currentTau = round(cfg.tau(channelpair)*ACT(channelpair,2,trials{channelpair,2}(t4t)));
		if currentTau < 1
			currentTau = 1;		% this may become zero for small ACT values -> overwrite manually
        end
	    
        % check if source needs to be embedded too
		if strcmp(cfg.embedsource,'no')
			pointset = TEembedding_noSourceEmb(a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.extracond);
		else
			pointset = TEembedding(a,b,cfg.dim(channelpair),currentTau,dimu(channelpair),cfg.extracond);
		end
                        
            % concatenate individually embedded trials in the first dimension
            pointsets_concat_2_aux   = cat(1,pointsets_concat_2_aux,pointset.pointset_2);
            pointsets_concat_p2_aux  = cat(1,pointsets_concat_p2_aux,pointset.pointset_p2);
            pointsets_concat_21_aux  = cat(1,pointsets_concat_21_aux,pointset.pointset_21);
            pointsets_concat_p21_aux = cat(1,pointsets_concat_p21_aux,pointset.pointset_p21);
	    if MIcalc
	        pointsets_concat_1_aux   = cat(1,pointsets_concat_1_aux,pointset.pointset_1);
            pointsets_concat_12_aux  = cat(1,pointsets_concat_12_aux,pointset.pointset_12);            
	    end
            clear pointset;
        end
        
        % stack embedded surrogate data in the 1st dimension
        pointsets_concat_2(cutpoint-chunksize+1:cutpoint,:)   = pointsets_concat_2_aux;
        pointsets_concat_p2(cutpoint-chunksize+1:cutpoint,:)  = pointsets_concat_p2_aux;
        pointsets_concat_21(cutpoint-chunksize+1:cutpoint,:)  = pointsets_concat_21_aux;        
        pointsets_concat_p21(cutpoint-chunksize+1:cutpoint,:) = pointsets_concat_p21_aux;
	if MIcalc
            pointsets_concat_1(cutpoint-chunksize+1:cutpoint,:)   = pointsets_concat_1_aux;
            pointsets_concat_12(cutpoint-chunksize+1:cutpoint,:)  = pointsets_concat_12_aux;
	end
	
        chunk_ind(cutpoint-chunksize+1:cutpoint)=ii+1;
        cutpoint = cutpoint + chunksize;
        
        clear *_aux;
        
    end
        
    if numpermutation > 0; 
        ft_progress('close');
        fprintf('\t - ok\n'); 
    end;
    
    fprintf('\nStarting GPU neighbour count ...\n');
    
    % remember indices of individual chunks
    cfg.chunk_ind = chunk_ind;    
       
	% get point counts from GPU functions
	[ncount] =  TEcallGPUsearch(cfg,channelpair,pointsets_concat_1,pointsets_concat_2, ...
        pointsets_concat_p2, pointsets_concat_21,pointsets_concat_12,pointsets_concat_p21);

    clear pointsets*; cfg = rmfield(cfg,'chunk_ind');
	
	% calculate TE and MI for orig and surrogate data (returns arrays of size nchunks = numpermutation + 1)
	[te mi tel] = TEcalc(cfg,ncount);
    clear ncount;    
	
    % get values for original data and add it to output structure
	TEmat(channelpair)       = te(1); 
	MImat(channelpair)       = mi(1);
    TEmat_sur(channelpair,:) = te(2:end);
	if cfg.TELcalc
		TELmat{channelpair} = tel;
    end	
    
    fprintf('\nCalculating Transfer Entropy')
	[p TE_diff] = TEpvalue(te,numpermutation);
    
	% do statistical comparison (is TE_orig an extreme value with respect to the TE values of the surrogate data?)
	% add results to output structure
    % 	1 - p_values of the statistic within the surrogate data
    %   2 - 1 (0), if the statistics is significant at the prescribed alpha level (or not)
    %   3 - 1 (0), if the statistics is significant after correction for multiple comparisons (or not)
    %   4 - difference between empirical TE value and median of the surrogate data distribution
    %   5 - 0, volume conduction (is set to 0 as the Faes Method is
    %       mandatory for using the ensemble method, see Faes, 2011)
	TEpermvalues(channelpair,1) = p;                                         % p value
	TEpermvalues(channelpair,2) = TEpermvalues(channelpair,1) < cfg.alpha;	 % statistical significance at alpha level
    % TEpermvalues(channelpair,3) is set below (significance after cmc)
	TEpermvalues(channelpair,4) = TE_diff;                                   % difference between empirical TE and median of surrogate distribution
	TEpermvalues(channelpair,5) = 0;                                         % Faes method is mandatory for ensemble method, takes care of volume conduction    
    fprintf(' - ok\n')
end



%% correct for multiple comparisons
% -------------------------------------------------------------------------

if ~findDelay
    fprintf('\nCorrection for multiple comparison...')
    pvalues                 = TEpermvalues(:,1);
    nrinstmix               = 0; % instanteneous mixing is not tested when ensemble method is used
    [significance,correctm] = TEcmc(pvalues, cfg.correctm, cfg.alpha, nrinstmix);
    TEpermvalues(:,3)       = significance;
    cfg.correctm            = correctm;    % the correction method might change, if only a small no. channels is analyzed, this should be updated
    fprintf(' - ok\n');
end

%% prepare output structure
% -------------------------------------------------------------------------

% add results to TEpermtest
TEpermtest = [];
TEpermtest.TEpermvalues   = TEpermvalues; clear TEpermvalues;
TEpermtest.dimord         = 'chanpair_value';
TEpermtest.cfg            = cfg;
TEpermtest.ACT.actvalue   = TEpreparestruct.ACT;
TEpermtest.sgncmb         = channelcombilabel;
TEpermtest.numpermutation = cfg.numpermutation;
TEpermtest.TEprepare      = TEpreparestruct;
TEpermtest.nr2cmc         = nr2cmc;
TEpermtest.TEmat          = TEmat;
TEpermtest.MImat          = MImat; 
TEpermtest.TEmat_sur      = TEmat_sur; 
TEpermtest.TELmat         = TELmat;

% add results to TEresult
TEresult           = [];
TEresult.TEmat     = TEmat; clear TEmat;
TEresult.MImat     = MImat; clear MImat;
TEresult.TEmat_sur = TEmat_sur; clear TEmat_sur;
TEresult.act       = ACT;
TEresult.dimord    = 'chanpair';
TEresult.cfg       = cfg;
TEresult.sgncmb    = channelcombilabel;
TEresult.TEprepare = TEpreparestruct; clear TEpreparestruct;

%% save results
% -------------------------------------------------------------------------


% fprintf('\nSaving ...')
% fprintf('\n\tresults of TE estimation')
% save(strcat(cfg.fileidout,'_time',num2str(cfg.toi(1)),'-',num2str(cfg.toi(2)),'s_TE_output.mat'), 'TEresult','-v7.3');
% fprintf(' - ok');
% fprintf('\n\tresults of permutation test')
% save(strcat(cfg.fileidout,'_time',num2str(cfg.toi(1)),'-',num2str(cfg.toi(2)),'s_TEpermtest_output.mat'), 'TEpermtest','-v7.3');
% fprintf(' - ok');


%% Returning to the working directory
cd(working_directory1)

return;

