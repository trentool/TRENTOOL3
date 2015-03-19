function TGA_results=InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data)

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
% The OUTPUT TGA_results results is a cell array of TEpermtest results,
% each cell containing the results for one prediction u 


% CHANGELOG
% 2013-25-01 PW: added GPU branch/ensemblemethod
% 2013-27-06 PW: added check for group analysis
% 2014-27-11 PW: changed interaction delay reconstruction: u is now optimized w/o surrogate
%		 testing, surrogate testing is done in a second step using the reconstructed
%		 interaction delay

%%%
if nargin ~= 3
    error('TRENTOOL: Please provide ''cfgTEP'', ''cfgTESS'' and ''data'' as input!');
end


%%% checks and parameter preparations

% check if data was prepared for later group statistics
groupanalysis = isfield(data,'groupprepare');
if groupanalysis
    
    groupprepare = data.groupprepare;
    
    fprintf('\n\nData was prepared for group statistics!\n')
    fprintf('The embedding dimension will be set to a common value for all datasets!\n\n')
    if isfield(cfgTEP,'predicttimemax_u') || isfield(cfgTEP,'predicttimemin_u') || isfield(cfgTEP,'predicttimestepsize')
        warning('Any parameter regarding the prediction time u, provided in cfgTEP will be overwritten by parameters from TEgroup_prepare.')
    end
    
    % set prediction time vector
    predicttimevec_u = data.groupprepare.predicttimevec_u;
    
    % set ragwitz dimension
    warning('Any parameter regarding the ragwitz dimension, provided in cfgTEP and cfgTESS will be overwritten by parameters from TEgroup_prepare.')
    fprintf('\tSetting ''cfgTEP.ragdim'' to maximum over all subjects (dim = %.0f)...\n',data.groupprepare.max_dim);
    fprintf('\tSetting ''cfgTESS.optdimusage'' to ''maxdim''...\n\n');
    cfgTEP.ragdim        = data.groupprepare.max_dim;
    cfgTESS.optdimusage  = 'maxdim';
else
    
    % cfg.predicttimevec_u supplied ?
    if ~isfield(cfgTEP,'predicttimemax_u')
    error(' No cfgTEP.predicttimemax_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end
    if ~isfield(cfgTEP,'predicttimemin_u')
        error(' No cfgTEP.predicttimemin_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end
    if ~isfield(cfgTEP,'predicttimestepsize')
        error(' No cfgTEP.predicttimestepsize specified - see HELP InteractionDelayReconstruction_calculate for more information');
    end

    predicttimevec_u=cfgTEP.predicttimemin_u:cfgTEP.predicttimestepsize:cfgTEP.predicttimemax_u;

    if ~(predicttimevec_u(end)==cfgTEP.predicttimemax_u)
        predicttimevec_u(end+1)=cfgTEP.predicttimemax_u; % make sure the last intended ptredictiontime is also investigated
    end
    % all other checks are left to the subsidiary functions
end

% construct output datastructure
TGA_results = cell(1,max(size(predicttimevec_u)));

t_total = tic;

%% TEprepare part
cfgTEP.predicttime_u = cfgTEP.predicttimemax_u;  % fix config
dataprep = TEprepare(cfgTEP,data);
clear data;

%% find optimal interaction delays
dataprep = TEfindDelay(predicttimevec_u,cfgTESS,dataprep);
cfgTESS.embedsource = 'yes';

%% calulate statistics with optimal u for individual channels
fprintf('Calling TEsurrogatestats using optimal delays u \n');

cfgTESS.fileidout=strcat(cfgTESS.fileidout,'_RAG4_TGA_opt_u');

% branch here for GPU calculation
if strcmp(dataprep.TEprepare.ensemblemethod,'yes')
    TGA_results=TEsurrogatestats_ensemble(cfgTESS,dataprep);
else
    TGA_results=TEsurrogatestats(cfgTESS,dataprep);
end

% add opt u vector to TEpermvalues matrix
TGA_results.TEpermvalues = [TGA_results.TEpermvalues TGA_results.TEprepare.u_in_ms];

if groupanalysis
    TGA_results.groupprepare = groupprepare;
end

t=toc(t_total);

fprintf('\n\n\nCALCULATION ENDED: %s \n', datestr(now));
fprintf('CALCULATION TOOK %.0f MINUTES (%.0f SECONDS)\n\n', t/60,t);
