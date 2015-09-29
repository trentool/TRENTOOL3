function TEgroup_conditionstatssingle(cfg, data1, data2)

% TECONDITIONSTATS: Performs a permutation test on results of transfer 
% entropy analysis for a single subject under two conditions (data1 and 
% data2). TEgroup_prepare has to be run on all data sets and TE values have 
% to be estimated using InteractionDelayReconstruction_calculate (or 
% TEsurrogatestats/_ensemble.m).
%
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!! The function TEgroup_prepare has to be run on all datasets first! !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% You can call this function directly as follows:
%         TEgroup_conditionstatssingle(cfg, data1, data2)
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
%     - The following Matlab toolboxes:
%         - statistic toolbox
%     - The functions
%         - TEperm
%
%
% * INPUT PARAMETERS
%
%   data1, data2 = results structures for two single subjects returned by 
%                  InteractionDelayReconstruction_calculate,
%                  TEsurrogatestats*; structures have to contain:
%       .TEmat        = no. channels x trials, matrix of TE values for all 
%                       channel conbinations and trials
%       .groupprepare = field with results from running TEgroup_prepare on
%                       all data sets
%       .sgncmb       = labels of channel combinations 
%
% AND
%
%   cfg: 
%
%   The configuration MUST contain:
%   cfg.fileidout   = string for the output path and the first part of the 
%                     output filename. The function will write a result 
%                     file with the suffix 'TEpermtestcondsingle_output.mat' 
%                     to disk  (containing the output of the statistical 
%                     test).
%
%   Optional arguments:
%
%   cfg.tail        = 1 tail or 2 tailed test of significance (for the
%                     permutation tests) (default = 2)
%   cfg.permstatstype  = 'mean' to use the distribution of the mean
%                     differences and 'depsamplesT' or
%                    'indepsamplesT' for distribution of the
%                     t-values. (default = 'mean')
%   cfg.alpha       = significance level for statisatical shift test,
%                     permutation test and correction for multiple
%                     comparison (default = 0.05)
%   cfg.correctm    = correction method used for correction of the multiple
%                     comparison problem - False discovery rate 'FDR' or
%                     Bonferroni correction 'BONF' (default = 'FDR')
%
% OUTPUT PARAMETERS
%
%  TEpermtestcondsingle
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
%                           6 - reconstructed interaction delay of the
%                               condition with stronger TE for that channel
%                               combination
%            .dimord        = dimensions of TEpermvalues
%            .cfg           = configuration file used to conduct the 
%                             permutation test
%            .sgncmb        = labels of channel combinations (source ->
%                             target)
%            .numpermutation = number of permutations
%            .TEgroupprepare = results of the function TEprepare from the
%                              data
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
% Version 3.0 by Patricia Wollstadt, Michael Lindner, Raul Vicente, 
% Michael Wibral
% Frankfurt 2015
%


%% check data
% -------------------------------------------------------------------------
fprintf('\nCheck data and config');


% check data for TEprepare/TEproupprepare structure
if ~isfield(data1, 'groupprepare') || ~isfield(data2, 'groupprepare')
    error('TRENTOOL ERROR: Field .groupprepare is missing in one of the datasets!');
end

% check code
if data1.groupprepare.code ~= data2.groupprepare.code
    error('TRENTOOL ERROR: Only datasets, which were prepared together with TEgroup_prepare can be analysed together! Code is not identical')
end

% get the max no. trials in any of the data sets that enter the test
cfg.maxtrials = max([max(data1.TEprepare.nrtrials(:,2)) max(data2.TEprepare.nrtrials(:,2))]);



%% compare new cfg and cfg from TEprepare if equal entries exist
% -------------------------------------------------------------------------

doublefields = 0;
cfgTEprepare = data1.TEprepare.cfg;

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

%% check cfg
% -------------------------------------------------------------------------

if ~isfield(cfg, 'alpha'),          cfg.alpha = 0.05;           end;
if ~isfield(cfg, 'correctm'),       cfg.correctm = 'FDR';       end;
if ~isfield(cfg, 'tail'),           cfg.tail = 2;               end;


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


nr2cmc = size(data1.TEpermvalues,1);

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


%% calculate statistics
% -------------------------------------------------------------------------

% perform permutation test
TEpermtestcondsingle = TEperm(cfg,data1,data2);

% add mean u value from group with bigger TE value
ind_u_1 = TEpermtestcondsingle.TEpermvalues(:,2) & TEpermtestcondsingle.TEpermvalues(:,4)>0;
ind_u_2 = TEpermtestcondsingle.TEpermvalues(:,2) & TEpermtestcondsingle.TEpermvalues(:,4)<0;
u_vec = sum([...
    data1.TEpermvalues(:,6) .* ind_u_1 ...
    data2.TEpermvalues(:,6) .* ind_u_2],2);
TEpermtestcondsingle.TEpermvalues = [TEpermtestcondsingle.TEpermvalues u_vec];

% add info to output structure
TEpermtestcondsingle.dimord         = 'chanpair_value'; 
TEpermtestcondsingle.cfg            = cfg;
TEpermtestcondsingle.sgncmb         = data1.sgncmb;
TEpermtestcondsingle.numpermutation = cfg.numpermutation;
TEpermtestcondsingle.TEgroupprepare = data1.groupprepare;

% keep some information from TEprepare (only fields that are relevant to all subjects)
TEpermtestcondsingle.TEgroupprepare.ensemblemethod    = data1.TEprepare.ensemblemethod;
TEpermtestcondsingle.TEgroupprepare.cfg               = data1.TEprepare.cfg;
TEpermtestcondsingle.TEgroupprepare.channelcombi      = data1.TEprepare.channelcombi;
TEpermtestcondsingle.TEgroupprepare.channelcombilabel = data1.TEprepare.channelcombilabel;
TEpermtestcondsingle.TEgroupprepare.channellabel      = data1.TEprepare.channellabel;
TEpermtestcondsingle.TEgroupprepare.timeindices       = data1.TEprepare.timeindices;
TEpermtestcondsingle.TEgroupprepare.optdim            = data1.TEprepare.optdim;

fprintf('\nCalculation finished\n')


%% save results (output of TEperm)
% -------------------------------------------------------------------------
toi = data1.cfg.toi;

fprintf('\nSaving ...')
fprintf('\n\tresults of permutation test')
save(...
    strcat(cfg.fileidout,'_time',num2str(toi(1)),'-',num2str(toi(2)),'s_TEpermtestcondsingle_output.mat'), ...
    'TEpermtestcondsingle','-v7.3');
fprintf(' - ok\n\n');



end


