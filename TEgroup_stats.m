function TEgroup_stats(cfg,filesTEpermtest)

% TEGROUP_STATS: The function uses a permutation test to compare TE values
% between two groups of subjects or between two conditions. Data have to be
% prepared by TEgroup_prepare and have to be analyzed using
% InteractionDelayReconstruction_calculate and 
% InteractionDelayReconstruction_analyze. 
%
% !!!!!!!! DATA HAS TO BE PREPARED AND ANALYZED FOR GROUP STATICIS !!!!!!!!
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
% * INPUT PARAMETERS
%
%   cfg.design      = matrix containing a row with subject number and a row
%                     with independent variable representing the order of
%                     the data input.
%                       example:
%                       datasets:    1 2 3 4 5 1 2 3 4 5
%                       conditions:  1 1 1 1 1 2 2 2 2 2
%   cfg.uvar        = row in cfg.design which contains the dataset number
%                     (in the example: 1)
%   cfg.ivar        = row in cfg.design which contains the independent
%                     variable (in the example: 2)
%   cfg.alpha       = significance level for statistical shift test,
%                     permutation test and correction for multiple
%                     comparison (default = 0.05)
%   cfg.numpermutation = nr of permutations in permutation test
%                     (default = 190100)
%   cfg.permstatstype  = 'mean' to use the distribution of the mean
%                     differences and 'depsamplesT' or
%                     'indepsamplesT' for distribution of the
%                     t-values. (default = 'mean')
%   cfg.tail        = '1' tail or '2' tailed test of significance (for the
%                     permutation tests) (default = 2)
%   cfg.correctm    = correction method used for correction of the multiple
%                     comparison problem - False discovery rate 'FDR' or
%                     Bonferroni correction 'BONF' (default = 'FDR')
%   cfg.fileidout   = string for the output path and the first part of the 
%                     output filename. TEgroup_stats will write two files
%                     with suffix 'TE_output.mat' (containing the raw 
%                     values used in the statistical test) and
%                     'TEpermtestgroup_output.mat' (containing the output
%                     of the statistical test).
%
%
% * OUTPUT PARAMETERS (ouput is saved to the location specified in
%                      cfg.fileidout)
%
%  TEpermtestgroup (saved in file '*TE_output.mat')
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
%                           4 - test statistic mean difference or tvalue 
%                               of mean difference depending on 
%                               cfg.permstatstype
%                           5 - 1 (0), if instantaneous mixing/volume
%                               conduction exists (or not)
%			    6 - mean estimated interaction delay u from the
%				group with the larger TE value
%
%            .dimord     = dimensions of TEpermvalues
%            .cfg        = configuration file used for group statistics
%            .sgncmb     = labels of channel combinations (source ->
%                          target)
%            .numpermutation = number of permutations
%            .nrdatasets     = number of datasets that entered the analysis
%            .TEgroupprepare = results of the function TEgroup_prepare 
%                              from the data
%
%  Raw data used for statistical testing (saved in file '*TE_output.mat')
%   TEresult1    = raw TE values for group/condition 1
%   TEresult2    = raw TE values for group/condition 2
%   TEresultmean = mean over raw values per group/condition
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2010
%

% CHANGELOG:
% 2013-06-25 PW: changed ival/uval to ivar/uvar in the function's help
% 2013-06-25 PW: nrchannels was not written into cfg structure, but is used
% later to construct the output structures
% 2013-10-20 PW: TEpermtestgroup.TEpermvalues now returns the mean interaction
% delay u for channel combinations with significant differences (mean value 
% from the group with the higher TE value)

%% Remember the working directory
working_directory1 = pwd;

%% load Data
% -------------------------------------------------------------------------

nrdatasets    = length(filesTEpermtest);
cfg.maxtrials = 0;                          % get the max. no. trials in any data set
allTEpermtest = cell(1,nrdatasets);         % holds all data for stat test
u_values = [];

for ll = 1:nrdatasets
    
    % load files containing 'TEpermteset'
    varinfile =  who('-file',filesTEpermtest{ll});
    load(filesTEpermtest{ll},varinfile{1});
    allTEpermtest{ll} = eval(varinfile{1});
    
    % check data for TEprepare/TEproupprepare structure
    if ~isfield(allTEpermtest{ll}, 'groupprepare')
        error(strcat('TRENTOOL ERROR: Field .groupprepare is missing in dataset ',filesTEpermtest{ll}));  
    end
        
    % check code 
    if allTEpermtest{1}.groupprepare.code ~= allTEpermtest{ll}.groupprepare.code
        error('TRENTOOL ERROR: Only datasets, which were prepared together with TEgroup_prepare can be analysed together! Code is not identical')
    end   
    
    % get the max no. trials in any of the data sets that enter the test
    if max(max(allTEpermtest{1}.TEprepare.nrtrials))>cfg.maxtrials
        cfg.maxtrials = max(max(allTEpermtest{1}.TEprepare.nrtrials));
    end

    u_values = cat(2,u_values,allTEpermtest{ll}.TEpermvalues(:,6));
        
end
clear ll 

% get number of channel combinations
nrchannelcombi = size(allTEpermtest{1}.sgncmb,1);


%% check cfg
% -------------------------------------------------------------------------

if ~isfield(cfg, 'alpha'),          cfg.alpha = 0.05;           end;
if ~isfield(cfg, 'correctm'),       cfg.correctm = 'FDR';       end;
if ~isfield(cfg, 'tail'),           cfg.tail = 2;               end;

if ~isfield(cfg, 'design'),
    error('\nTRENTOOL ERROR: cfg.design must be defined, see help!');
end;
if ~isfield(cfg, 'ivar'),
    error('\nTRENTOOL ERROR: cfg.ivar must be defined, see help!');
end;
if ~isfield(cfg, 'uvar'),
    error('\nTRENTOOL ERROR: cfg.uvar must be defined, see help!');
end;
if size(cfg.design,2) ~= nrdatasets
    error('\nTRENTOOL ERROR: cfg.design must have the same number of columns than number of datasets exist, see help!');
end

if ~isfield(cfg, 'permstatstype'),  cfg.permstatstype = 'mean'; end;
if ~strcmp(cfg.permstatstype , 'mean') && ~strcmp(cfg.permstatstype , 'indepsamplesT') && ~strcmp(cfg.permstatstype , 'depsamplesT')
    error('\nTRENTOOL ERROR: wrong cfg.permstatstype - use ''mean'' ''depsamplesT'' or ''indepsamplesT'', see help!');
end

if ~isfield(cfg, 'fileidout'),
    error('\nTRENTOOL ERROR: cfg.fileidout must be defined, see help!');
end;


%% check nr of permutations
% -------------------------------------------------------------------------
fprintf('\n\nChecking number of permutations');


% cfg.permtest.channelcombi = channelcombi;
% cfg.permtest.channelcombilabel = data.TEprepare.channelcombilabel ;


nr2cmc = nrchannelcombi*length(allTEpermtest{1}.groupprepare.predicttimevec_u);

if ~isfield(cfg, 'numpermutation'),
    cfg.numpermutation = 190100; % for p<0.01 with a possible bonferroni correcetion of 100
elseif cfg.numpermutation < ceil(1/cfg.alpha)
    error(strcat('\nTRENTOOL ERROR: cfg.numpermutation too small - Nr of permutations must be at least :',num2str(numpermutation),' !'));
else
    if cfg.maxtrials>31
        if cfg.numpermutation > 2^31
            error(strcat('\nTRENTOOL ERROR: cfg.numpermutation too huge - Nr of permutations must be at least :',num2str(numpermutation),' !'));
        end
    else
        if cfg.numpermutation > 2^cfg.maxtrials
            error(strcat('\nTRENTOOL ERROR: cfg.numpermutation too huge - Nr of permutations must be at least :',num2str(numpermutation),' !'));
        end
    end
    if cfg.numpermutation < ceil(1/(cfg.alpha/nr2cmc))
       fprintf('\n#######################################################################################\n# WARNING: Nr of permutations not sufficient for correction for multiple comparisons! #\n#######################################################################################\n'); 
    end
end

fprintf(' - ok\n');



%% Define conditions
% -------------------------------------------------------------------------
fprintf('\nDefine conditions:\n');

conds    = squeeze(cfg.design(cfg.ivar,:));
condtype = unique(conds);
nrconds  = length(condtype);
if nrconds ~=2
    fprintf(strcat(['TRENTOOL ERROR: You defined ',num2str(nrconds),' conditions.']))
    error('\nTRENTOOL ERROR: You have to define two conditions in cfg.design, see help!');
end

units    = squeeze(cfg.design(cfg.uvar,:));
unittype = unique(units);
nrunits  = length(unittype);

condindex1 = find(conds == condtype(1));
condindex2 = find(conds == condtype(2));

fprintf('Total no. datasets: %d\n', nrunits)
fprintf('\tCondition 1 -> ')
fprintf('%d ', condindex1 )
fprintf('\n\tCondition 2 -> ')
fprintf('%d ', condindex2)
fprintf('\n')

u_mean1 = mean(u_values(:,condindex1),2);
u_mean2 = mean(u_values(:,condindex2),2);
clear u_values;

%% Create TEresultmean (mean over trials per subject)
% -------------------------------------------------------------------------

% create empty matrices
TEresultmean.TEmat = zeros(nrchannelcombi,nrdatasets);
%TEresultmean.MImat = zeros(nrchannelcombi,nrdatasets);

% TE value per subject and significant channel combi (mean over trials)
for subject = 1:nrdatasets
    
    % find all significant TE values
    % sign_ind = allTEpermtest{subject}.TEmat(:,CMC_ind);
    
    TEresultmean.TEmat(:,subject) = squeeze(mean(allTEpermtest{subject}.TEmat, 2));
    %TEresultmean.MImat(:,subject) = squeeze(mean(MImat{subject}, 2));
end


%% permutation tests
% -------------------------------------------------------------------------

% split TEresultmean matrices for permutation tests
fprintf('Prepare condition matrices');

TEresult1.TEmat = TEresultmean.TEmat(:,condindex1);
%TEresult1.MImat = TEresultmean.MImat(:,condindex1);

TEresult2.TEmat = TEresultmean.TEmat(:,condindex2);
%TEresult2.MImat = TEresultmean.MImat(:,condindex2);

fprintf(' - ok');

% perform permutation test
TEpermtestgroup = TEperm(cfg,TEresult1,TEresult2);

% add mean u value from group with bigger TE value
ind_u_1 = TEpermtestgroup.TEpermvalues(:,2) & TEpermtestgroup.TEpermvalues(:,4)>0;
ind_u_2 = TEpermtestgroup.TEpermvalues(:,2) & TEpermtestgroup.TEpermvalues(:,4)<0;
u_vec = sum([u_mean1.*ind_u_1 u_mean2.*ind_u_2],2);
TEpermtestgroup.TEpermvalues = [TEpermtestgroup.TEpermvalues u_vec];

% add info to output structure
TEpermtestgroup.dimord         = 'chanpair_value'; 
TEpermtestgroup.cfg            = cfg;
TEpermtestgroup.sgncmb         = allTEpermtest{1}.sgncmb;
TEpermtestgroup.numpermutation = cfg.numpermutation;
TEpermtestgroup.nrdatasets     = nrdatasets;
%TEpermtestgroup.TEprepare      = allTEpermtest{1}.TEprepare;
TEpermtestgroup.TEgroupprepare = allTEpermtest{1}.groupprepare;
% keep some information from TEprepare (only fields that are relevant to all subjects)
TEpermtestgroup.TEgroupprepare.ensemblemethod    = allTEpermtest{1}.TEprepare.ensemblemethod;
TEpermtestgroup.TEgroupprepare.cfg               = allTEpermtest{1}.TEprepare.cfg;
TEpermtestgroup.TEgroupprepare.channelcombi      = allTEpermtest{1}.TEprepare.channelcombi;
TEpermtestgroup.TEgroupprepare.channelcombilabel = allTEpermtest{1}.TEprepare.channelcombilabel;
TEpermtestgroup.TEgroupprepare.channellabel      = allTEpermtest{1}.TEprepare.channellabel;
TEpermtestgroup.TEgroupprepare.timeindices       = allTEpermtest{1}.TEprepare.timeindices;
TEpermtestgroup.TEgroupprepare.optdim            = allTEpermtest{1}.TEprepare.optdim;

fprintf('\nCalculation finished\n')


%% save results
% -------------------------------------------------------------------------
toi = allTEpermtest{1}.cfg.toi;

fprintf('\nSaving ...')

% this saves data for each group, that enters the group statistics
fprintf('\n\tdata used in permutation test')
savename1 = strcat(cfg.fileidout,'_time',num2str(toi(1)),'-',num2str(toi(2)),'s_TEpermtestgroup_data.mat');  
save(savename1, 'TEresultmean','TEresult1','TEresult2','-v7.3');
fprintf(' - ok');

% this saves the output of the group statistics
fprintf('\n\tresults of permutation test')
save(strcat(cfg.fileidout,'_time',num2str(toi(1)),'-',num2str(toi(2)),'s_TEpermtestgroup_output.mat'), 'TEpermtestgroup','-v7.3');
fprintf(' - ok\n\n');


%% Returning to the working directory
cd(working_directory1)


end
