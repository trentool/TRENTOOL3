function TEpermtest = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data)

% InteractionDelayReconstruction_calculate
% Graph Timing calculation based on TE with scanned prediction times u
% 
% This function has to be supplied with the folowing INPUT:
%
% data      -    the data in the fieldtrip raw data format; see TEprepare for more
%           details
% cfgTEP    -    a configuration structure that has all the fields required for
%           TEprepare; see TEprepare for more details. The only difference is that instead of
%           cfg.predictime_u which is a single number, cfg.predicttimemin_u
%           cfg.predicttimemax_u,cfg.predicttimestepsize
%           have to be supplied, indicating the minimum and maximum prediction time of
%           inteterest and the stepsize of the resolution
% cfgTSS    -   a configuration structure that has all the fields required
%           for TEsurrogatestats
%
% The OUTPUT TEpermtest is a structure containing estimated TE values and 
% results for permutation testing using the optimal interaction delay for 
% each channel combination.


% CHANGELOG
% 2013-25-01 PW: added GPU branch/ensemblemethod
% 2013-27-06 PW: added check for group analysis
% 2014-27-11 PW: changed interaction delay reconstruction: u is now optimized w/o surrogate
%		 testing, surrogate testing is done in a second step using the reconstructed
%		 interaction delay
% 2013-27-06 PW: TEpermtest is saved to disk by this function


%% check input
if nargin ~= 3
    error('TRENTOOL: Please provide ''cfgTEP'', ''cfgTESS'' and ''data'' as input!');
end

%% define logging levels
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;

if ~isfield(cfgTEP, 'verbosity'), cfgTEP.verbosity = 'info_minor'; end;

%% checks and parameter preparations

% check if data was prepared for later group statistics
groupanalysis = isfield(data,'groupprepare');
if groupanalysis
    
    groupprepare = data.groupprepare;
    
    msg = ['Data was prepared for group statistics; the embedding ' ...
        'dimension will be set to a common value for all datasets!'];
    TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MINOR);
    if isfield(cfgTEP,'predicttimemax_u') || isfield(cfgTEP,'predicttimemin_u') || isfield(cfgTEP,'predicttimestepsize')
        fprintf('\n')
        warning(['TRENTOOL WARNING: Any parameter regarding the ' ...
            'prediction time u, provided in cfgTEP will be ' ...
            'overwritten by parameters from TEgroup_prepare.'])
    end
    
    % set prediction time vector
    predicttimevec_u = data.groupprepare.predicttimevec_u;
    
    % set ragwitz dimension   
    fprintf('\n')
    warning(['TRENTOOL WARNING: Any parameter regarding the ' ...
            'ragwitz dimension, provided in cfgTEP/cfgTESS will be ' ...
            'overwritten by parameters from TEgroup_prepare.'])
    msg = sprintf('Setting ''cfgTEP.ragdim'' to maximum over all subjects (dim = %.0f)',data.groupprepare.max_dim);
    TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MINOR);
    msg = 'Setting ''cfgTESS.optdimusage'' to ''maxdim'' ';
    TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MINOR);
    cfgTEP.ragdim        = data.groupprepare.max_dim;
    cfgTESS.optdimusage  = 'maxdim';
else
    
    % cfg.predicttimevec_u supplied ?
    if ~isfield(cfgTEP,'predicttimemax_u')
        fprintf('\n')
        error('TRENTOOL ERROR: No cfgTEP.predicttimemax_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end
    if ~isfield(cfgTEP,'predicttimemin_u')
        fprintf('\n')
        error('TRENTOOL ERROR: No cfgTEP.predicttimemin_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end
    if ~isfield(cfgTEP,'predicttimestepsize')
        fprintf('\n')
        error('TRENTOOL ERROR: No cfgTEP.predicttimestepsize specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end

    predicttimevec_u=cfgTEP.predicttimemin_u:cfgTEP.predicttimestepsize:cfgTEP.predicttimemax_u;

    if ~(predicttimevec_u(end)==cfgTEP.predicttimemax_u)
        predicttimevec_u(end+1)=cfgTEP.predicttimemax_u; % make sure the last intended ptredictiontime is also investigated
    end
    % all other checks are left to the subsidiary functions
end


t_total = tic;

%% TEprepare part
msg = '################### PREPARING DATA FOR TE ANALYSIS';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

cfgTEP.predicttime_u = cfgTEP.predicttimemax_u;  % fix config
dataprep = TEprepare(cfgTEP,data);
clear data;

%% find optimal interaction delays
msg = '################### OPTIMIZING INFORMATION TRANSFER DELAY';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

[dataprep, TEmat] = TEfindDelay(predicttimevec_u,cfgTESS,dataprep);
cfgTESS.embedsource = 'yes';

%% calulate statistics with optimal u for individual channels
msg = '################### ESTIMATING TRANSFER ENTROPY WITH OPTIMIZED PARAMETERS';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

cfgTESS.fileidout=strcat(cfgTESS.fileidout,'_RAG4_TGA_opt_u');

% branch here for GPU calculation
if strcmp(dataprep.TEprepare.ensemblemethod,'yes')
    TEpermtest=TEsurrogatestats_ensemble(cfgTESS,dataprep);
else
    TEpermtest=TEsurrogatestats(cfgTESS,dataprep);
end

% add opt u vector to TEpermvalues matrix
TEpermtest.TEpermvalues = [TEpermtest.TEpermvalues TEpermtest.TEprepare.u_in_ms];
TEpermtest.TEbyU        = TEmat;

if groupanalysis
    TEpermtest.groupprepare = groupprepare;
end

%% save results

msg = 'Saving results of TE estimation and surrogate testing';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MINOR);
save(strcat(cfgTESS.fileidout,'_time',num2str(cfgTEP.toi(1)),'-',num2str(cfgTEP.toi(2)),'s_TEpermtest_output.mat'), ...
    'TEpermtest','-v7.3');

%%

t=toc(t_total);
msg = sprintf( ...
    'Thank you for using this transfer entropy tool!\n\nTRANSFER ENTROPY CALCULATION ENDED: %s \nCALCULATION TOOK %.0f MINUTES (%.0f SECONDS)', ...
    datestr(now), t/60, t);
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);
