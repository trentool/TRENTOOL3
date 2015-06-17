function TEconditionstatssingle(cfg,varargin)

% TECONDITIONSTATS: This function calculates the transferentropy values and
% performs a permutation test on two transfer entropy data sets of two
% conditions.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!    The function TEprepare has to be run on all datasets first!    !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% You can call this function directly as follows:
%         TEconditionstats(cfg, data1, data2)
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
%
% * INPUT PARAMETERS
%
%   data1 - datan   = datasets e.g. from subjects and conditions
%                     Fieldtrip raw data structures which MUST contain:
%       .trials     = three dimensional datamatrix
%       .time       = vector 1 x numtoi, the time included in the
%                     data
%       .label      = vector 1 x numlabe, labels of channelnames included
%                     in the data
%       .fsample    = value of sampling rate
%
% AND
%
%   cfg: The configuration MUST contain:
%
%   cfg.dim         = Value(s) for embedding dimension. In case of using
%                     cfg.optdimusage = 'maxdim' this has to be a scalar
%                     value. In case of cfg.optdimusage = 'indivdim' this 
%                     has to be a cell of the size {number of subjects x1} 
%                     including vectors of the size (channelcombi x 1). 
%                     If not specified, the optimal dimension(s) found in 
%                     TEprepare will be used, which is the recommended 
%                     option!
%   cfg.tau         = embedding delay in units of act (x*act). If not
%                     specified (recommended option), the tau is used as
%                     followed:
%                     In case of optimizemethod in TEprepare:
%                           'ragwitz' = optimal tau found via ragwitz
%                                       critrion
%                           'cao'     = cfg.tau given by user in TEprepare
%
%   cfg.alpha       = significance level for statisatical shift test,
%                     permutation test and correction for multiple
%                     comparison (default = 0.05)
%   cfg.shifttest   = perform shift test to identify instantaneous mixing
%                     between the signal pairs. Values: 'yes' or 'no'
%                     (default = 'yes')
%                     This shift test is important for EEG and MEG data,
%                     because linear mixing is always present in the data.
%                     In case of instantaneous mixing transfer entropy
%                     should not be calculated for the affected
%                     channelpairs with the corresponding parameter sets,
%                     because it could result in false positive results.
%                     Hence the TE values for these cases will be set to
%                     NaN and the corresponding p-values of the permutation
%                     test to 1.
%   cfg.shifttesttype = The shift test can be calculated for the direction
%                     TE value of original data > TE values of shifted data
%                     (value = 'TE>TEshift') or for the other direction
%                     (value = 'TEshift>TE'). In this case the alpha is
%                     set to 0.1. (default = 'TE>TEshift')
%   cfg.shifttype     = Shifting the data 'onesample' or the length of the
%                      'predicttime' (default = 'predicttime')
%   cfg.numpermutation = nr of permutations in permutation test
%                     (default = 190100)
%   cfg.permstatstype  = 'mean' to use the distribution of the mean
%                     differences and 'depsamplesT' or
%                    'indepsamplesT' for distribution of the
%                     t-values. (default = 'mean')
%   cfg.tail        = 1 tail or 2 tailed test of significance (for the
%                     permutation tests) (default = 2)
%   cfg.correctm    = correction method used for correction of the multiple
%                     comparison problem - False discovery rate 'FDR' or
%                     Bonferroni correction 'BONF' (default = 'FDR')
%   cfg.fileidout   = string for the first part of the output filename.
%
%
%
% OUTPUT PARAMETERS
%
%
%  TEpermtestcondition
%            .TEpermvalues  = matrix with size:
%                             (channelpair x value)
%                           The last dimension "value" includes:
%                           1 - p_values of the statistic within the
%                               distribution given by the permutations
%                           2 - 1 (0), if the statistics is significant at
%                               the prescribed alpha level (or not)
%                           3 - 1 (0), if the statistics is significant
%                               after correction for mulitple comparisons
%                               (or not)
%                           4 - 1 (0), mean difference or tvalue of mean
%                               difference depending on cfg.permstatstype
%                           5 - 1 (0), if instantaneous mixing (volume
%                               conduction) exists (or not)
%            .dimord        = dimensions of TEpermvalues
%            .cfg           = configuration file used to calculate TE and
%                             permtest
%            .sgncmb        = labels of channel combinations (source ->
%                             target)
%            .numpermutation = number of permutations
%            .ACT           = structure including
%                .actvalue  = ACT matrix (channelcombi x 2 x trial)
%            .TEprepare     = results of the function TEprepare fron the
%                             data
%
%
% AND
%
%
%  TEresult             = Output structure of the function tranferentropy
%          .TEmat       = resultmatrix including transfer entropy(TE)
%                         values. (channelpairs x trial)
%          .MImat       = resultmatrix including mutual information (MI)
%                         values. (channelpairs x trial)
%          .dimord      = 'channelpair_u_trial'; the dimensions of TEmat
%                         and MImat
%          .cfg         = configuration file used to calculate TE
%          .trials      = trial numbers selected from raw dataset
%          .act         = ACT matrix (channelcombi x 2 x trial)
%          .sgncmb      = labels of channel combinations (source -> target)
%          .TEprepare   = results of the function TEprepare from the
%                         data
%   if instantaneous mixing is found in the data, then another field will
%   be added:
%          .instantaneousmixing = matrix (channelcombi) which indicates were
%          the instantaneous mixings were found (1) or not (0).
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 2.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Bonn 2011
%

%% Remember the working directory
working_directory1 = pwd;


warning('!!!!!!!! THIS FUNCTION IS DEPRECATED - USE TEgroup_prepare and TEgroup_stats instead !!!!!!!!');
pause(5);


%% check data
% -------------------------------------------------------------------------
fprintf('\nCheck data and config');

if length(varargin) > 2
    error('\nTRENTOOL ERROR: Only two datasets can be compared. More than two datasets are given as input. (see help)')
end

% check parameters of the input data for equality

DataCell2Comp{2}=[];
for ii = 1:2
    DataCell2Comp{ii}=varargin{ii};
end
groupres = TEpreparegroup(DataCell2Comp);


% compare new cfg and cfg from TEprepare if equal entries exist
% -------------------------------------------------------------------------

doublefields = 0;
cfgTEprepare = varargin{1}.TEprepare.cfg;
if isfield(cfg, 'Path2TSTOOL')
    cfgTEprepare = rmfield(cfgTEprepare, 'PathTSTOOL');
end

cfgfields = fieldnames(cfgTEprepare);
cfgfields2 = fieldnames(cfg);


for ii = 1:size(cfgfields,1);
    for jj = 1:size(cfgfields2,1);
        if strcmp(cfgfields{ii},cfgfields2{jj})
            doublefields = doublefields + 1;
        end
    end
end

if doublefields  > 0
    fprintf('\n')
    error('TRENTOOL ERROR: Illegal attempt to overwrite entry generated by or used for TEprepare! Change cfg or rerun TEprepare. (see help)')
end


% add structures and values of data.TEprepare.cfg to cfg
nr1 = size(cfgfields,1);
for ii = 1:nr1
    eval(strcat('cfg.',cfgfields{ii},' = getfield(cfgTEprepare, {1}, cfgfields{ii});'));
end



% check configuration and set defaults
% -------------------------------------------------------------------------

cfg.nrdatasets = size(varargin,2);

% if not defined set defaults
if ~isfield(cfg, 'alpha'),          cfg.alpha = 0.05;           end;
if ~isfield(cfg, 'correctm'),       cfg.correctm = 'FDR';       end;
if ~isfield(cfg, 'tail'),           cfg.tail = 2;               end;


if ~isfield(cfg, 'permstatstype'),  cfg.permstatstype = 'mean'; end;
if strcmp(cfg.permstatstype , 'mean') == 0 && strcmp(cfg.permstatstype , 'indepsamplesT') == 0 && strcmp(cfg.permstatstype , 'depsamplesT') == 0
    fprintf('\n');
    error('TRENTOOL ERROR: wrong cfg.permstatstype - use ''mean'' ''depsamplesT'' or ''indepsamplesT'', see help!');
end

if ~isfield(cfg, 'shifttest'),  cfg.shifttest = 'yes'; end;
if strcmp(cfg.shifttest , 'yes') == 0 && strcmp(cfg.shifttest , 'no') == 0
    fprintf('\n');
    error('TRENTOOL ERROR: wrong cfg.shifttest - use ''yes'' or ''no'', see help!');
end
if strcmp(cfg.shifttest , 'yes')
    if ~isfield(cfg, 'shifttype'),    cfg.shifttype = 'predicttime'; end;
    if ~isfield(cfg, 'shifttesttype'),  cfg.shifttesttype = 'TE>TEshift'; end;
    if strcmp(cfg.shifttesttype , 'TE>TEshift') == 0 && strcmp(cfg.shifttesttype , 'TEshift>TE') == 0
        fprintf('\n');
        error('TRENTOOL ERROR: wrong cfg.shifttesttype - use ''TE>TEshift'' or ''TEshift>TE'', see help!');
    end
end

if ~isfield(cfg, 'fileidout'),
    fprintf('\n');
    error('TRENTOOL ERROR: cfg.fileidout must be defined, see help!');
end;


% check dim
if ~isfield(cfg, 'dim'), cfg.dim=groupres.groupoptdim;  end;
if cfg.dim < groupres.groupoptdim
    fprintf('\n');
    fprintf('TRENTOOL WARNING: specified embedding dimension (cfg.dim) is smaller then the optimal dimension from TEprepare.')
elseif cfg.dim > groupres.groupoptdim
    fprintf('\n');
    fprintf('TRENTOOL WARNING: specified embedding dimension (cfg.dim) is bigger then the optimal dimension from TEprepare.')
end
% cfg.dim to vector
cfg.dim(1:groupres.nrchannelcombis,1) = cfg.dim;


% check tau
if strcmp(varargin{1}.TEprepare.cfg.optimizemethod, 'ragwitz')
    cfg.tau(1:size(varargin{1}.TEprepare.channelcombi,1)) = groupres.groupopttau;
elseif strcmp(varargin{1}.TEprepare.cfg.optimizemethod, 'cao')
    cfg.tau(1:size(varargin{1}.TEprepare.channelcombi,1)) = varargin{1}.TEprepare.cfg.caotau;
end

% check TE parameter
if isempty(cfg.predicttime_u), 
    fprintf('\n');
    error('TRENTOOL ERROR: cfg.predicttime_u missing in .TEprepare!');  
end;

if ~isfield(cfg, 'kth_neighbors'),  cfg.kth_neighbors = 4;  end;

if ~isfield(cfg, 'TheilerT'),       cfg.TheilerT = 'ACT';   end;
if ~strcmp(cfg.TheilerT, 'ACT');
    if size(cfg.TheilerT,1)>1 || size(cfg.TheilerT,2)>1
        fprintf('\n');
        error('TRENTOOL ERROR: cfg.TheilerT must include a scalar, see help!');
    end
end


% check the format of input vectors
if size(cfg.toi,1)>size(cfg.toi,2)
    cfg.toi=cfg.toi';
elseif size(cfg.kth_neighbors,1)>1 || size(cfg.kth_neighbors,2)>1
    fprintf('\n');
    error('TRENTOOL ERROR: cfg.kth_neighbors must include a scalar, see help!');
end

fprintf(' - ok');



%% Loop over datasets - calculate TE and shifttest for each separately
% -------------------------------------------------------------------------

% create empty result cell
TEresult{cfg.nrdatasets}=[];

% start loop
for subject = 1:cfg.nrdatasets

    fprintf(strcat(['\nDataset ',num2str(subject),' of ',num2str(cfg.nrdatasets),'\n\n']))



    %% get channels, ACT and trials from the cfg.TEprepare
    % ------------------------------------------------------------------------

    fprintf('\nGet channels, trials, and ACT info from TEprepare')

    % select channels
    channelcombi=varargin{subject}.TEprepare.channelcombi ;
    cfg.permtest.channelcombi = channelcombi;
    cfg.permtest.channelcombilabel = varargin{subject}.TEprepare.channelcombilabel ;
    
    cfg.permtest.ACT{subject}=varargin{subject}.TEprepare.ACT;

    % select trials
    trials=varargin{subject}.TEprepare.trials;
    nrtrials=varargin{subject}.TEprepare.nrtrials;
    cfg.permtest.trials{subject}=trials;
    cfg.permtest.nrtrials{subject}=nrtrials;

    


    fprintf(' - ok');

    %% check nr of permutations
    % -------------------------------------------------------------------------
    fprintf('\nChecking number of permutations');

    nr2cmc=size(channelcombi,1);

    if ~isfield(cfg, 'numpermutation'),
        cfg.numpermutation = 190100; % for p<0.01 with a possible bonferroni correcetion of 100
    elseif cfg.numpermutation < ceil(1/cfg.alpha)
        fprintf('\n');
        error('TRENTOOL ERROR: cfg.numpermutation too small!');
    else
        if nrtrials>31
            if cfg.numpermutation > 2^31
                fprintf('\n');
                error('\nTRENTOOL ERROR: cfg.numpermutation too huge!');
            end
        else
            if cfg.numpermutation > 2^min(min(nrtrials))
                fprintf('\n');
                error('\nTRENTOOL ERROR: cfg.numpermutation too huge!');
            end
        end
        if cfg.numpermutation < ceil(1/(cfg.alpha/nr2cmc))
            fprintf('\n#######################################################################################\n# WARNING: Nr of permutations not sufficient for correction for multiple comparisons! #\n#######################################################################################\n');
        end
    end

    fprintf(' - ok\n');


    %% start calculating TE
    % -------------------------------------------------------------------------

    if subject == 1
        cfg.calctime = 'yes';
    else
        cfg.calctime = 'no';
    end

    % for unshuffled data
    % ----------------------
    fprintf('\nStart calculating transfer entropy for unshuffled data');
    cfg.shuffle = 'no';
    [TEresult{subject}] = transferentropy(cfg,varargin{subject});
    TEresult{subject}.TEprepare = varargin{subject}.TEprepare;


    cfg.calctime = 'no';

    % for shifted data
    % ----------------------
    % TEshift is created inside transferentropy.m as a reduced version of
    % TEresult without certain fields. TEshift is never written to disk/file
    % to avoid later confusion. Please save TEshift yourself if necessary.
    if strcmp(cfg.shifttest, 'yes')
        fprintf('\nStart calculating transfer entropy for shifted data');
        cfg.shuffle = 'no';
        [TEshift] = transferentropy(cfg,varargin{subject},'shifttest');

        % permutation test for shift test
        fprintf('\nStart permutation tests for shift test');
        permstatstype = cfg.permstatstype;
        cfg.permstatstype = 'indepsamplesT';
        tailtype = cfg.tail;
        cfg.tail = 1;
        if strcmp(cfg.shifttesttype, 'TE>TEshift')
            alpha = cfg.alpha;
            cfg.alpha = 0.05;
            TEpermshift = TEperm(cfg,TEresult{subject},TEshift);
            cfg.alpha = alpha;
        elseif strcmp(cfg.shifttesttype, 'TEshift>TE')
            alpha = cfg.alpha;
            cfg.alpha = 0.1;
            TEpermshift = TEperm(cfg,TEshift,TEresult{subject});
            cfg.alpha = alpha;
        end
        cfg.permstatstype = permstatstype;
        cfg.tail=tailtype;

        % analyze shift test
        fprintf('\nanalyze shift test\n');
        if strcmp(cfg.shifttesttype, 'TE>TEshift')
            indexinstmix = find(TEpermshift.TEpermvalues(:,2)==0);
            if size(indexinstmix,1) == 0
                fprintf('No instantaneous mixing found!\n')
            else
                fprintf(strcat(num2str(size(indexinstmix,1)),' instantaneous mixings found!\nFor these cases TEvalues of all trials are set to NaN!\n'))
                mask=repmat((TEpermshift.TEpermvalues(:,2)-1)*-1, [1 1 size(TEresult{subject}.TEmat,2)]);
                TEresult{subject}.TEmat(mask==1) = NaN;
                TEresult{subject}.MImat(mask==1) = NaN;
                clear mask;
                TEresult{subject}.instantaneousmixing = (TEpermshift.TEpermvalues(:,2)-1)*-1;
            end
        elseif strcmp(cfg.shifttesttype, 'TEshift>TE')
            indexinstmix = find(TEpermshift.TEpermvalues(:,2)==1);
            if size(indexinstmix,1) == 0
                fprintf('No instantaneous mixing found!\n')
            else
                fprintf(strcat(num2str(size(indexinstmix,1)),' instantaneous mixings found!\nFor these cases TEvalues of all trials are set to NaN!\n'))
                mask=repmat(TEpermshift.TEpermvalues(:,2), [1 1 size(TEresult{subject}.TEmat,2)]);
                TEresult{subject}.TEmat(mask==1) = NaN;
                TEresult{subject}.MImat(mask==1) = NaN;
                clear mask;
                TEresult{subject}.instantaneousmixing = TEpermshift.TEpermvalues(:,2);
            end
        end

        clear TEpermshift TEshift
    end



    cfg = rmfield(cfg, 'calctime');

end



%% permutation tests
% -------------------------------------------------------------------------
fprintf('\nStart permutation tests');

%TEpermtestcondition=[];

TEpermtestcondition = TEperm(cfg,TEresult{1},TEresult{2});

TEpermtestcondition.dimord = 'chanpair_value';
TEpermtestcondition.cfg = cfg;
% TEpermtestcondition.ACT.label = varargin{subject}.label;
%TEpermtestcondition.ACT.actvalue = ACT;
TEpermtestcondition.sgncmb = TEresult{subject}.sgncmb;
TEpermtestcondition.numpermutation = cfg.numpermutation;
TEpermtestcondition.nrdatasets = cfg.nrdatasets;
TEpermtestcondition.TEprepare = varargin{subject}.TEprepare;

fprintf('\nCalculation ready\n')


%% save results
% -------------------------------------------------------------------------
fprintf('\nSaving ...')
fprintf('\nResults of TE')
savename1 = strcat(cfg.fileidout,'_time',num2str(cfg.toi(1)),'-',num2str(cfg.toi(2)),'ms_TE_output.mat');
save(savename1, 'TEresult','-v7.3');
fprintf(' - ok');
fprintf('\nResults of permutation test')
save(strcat(cfg.fileidout,'_time',num2str(cfg.toi(1)),'-',num2str(cfg.toi(2)),'ms_TEpermtest_output.mat'), 'TEpermtestcondition','-v7.3');
fprintf(' - ok');


%% Returning to the working directory
cd(working_directory1)

fprintf('\n\nThank you for using this transfer entropy tool!\n')

return;

