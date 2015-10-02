function [data, TEmat]= TEfindDelay(predicttimevec_u,cfgTESS,data)

% FUNCTION TEFINDDELAY
%
% Reconstructs optimal interaction delay u from a vector of assumed delays.
% The function estimates TE for all channels and all interaction delays and
% returns the optimal interaction delay as the value for u that maximizes
% the TE. 
%
% INPUT
%   predicttimevector_u - vector with interaction delays to be scanned
%   cfgTESS             - cfg structure for calls to TEsurrogatestats or
%                         TEsurrogatestats_ensemble
%   data   		- data that has already undergone preparation by TEprepare
%
% OUTPUT
%   data - data, where the prediction time in fields of data.TEprepare are
%	   replaced by a vector of the optimal interaction delay for each
%	   channel combination. The following fields are changed by the
%	   function:
%		data.TEprepare.u_in_ms
%		data.TEprepare.u_in_samples
%		data.TEprepare.cfg.predicttime_u
%   TEmat - array with size [n_sgncmb x u] raw TE values for each signal 
%      combination and assumed delay u, returned by TEfindmaxte.m
%
% PW 27/11/2014

%% define logging levels
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;

verbosity = data.TEprepare.cfg.verbosity;

%%
cfgTESS.numpermutation = 'findDelay';
TGA_results = cell(1, length(predicttimevec_u));
predicttimevec_u_samples = round(predicttimevec_u/1000*data.fsample);
fileidout = cfgTESS.fileidout;
n_channelcombis = size(data.TEprepare.channelcombilabel,1);

for uu=1:max(size(predicttimevec_u))
    
    msg = sprintf('Estimating TE for u = %.0f ms', predicttimevec_u(uu));
    TEconsoleoutput(data.TEprepare.cfg.verbosity, msg, LOG_INFO_MINOR);
    
    data.TEprepare.cfg.predicttime_u = repmat(predicttimevec_u(uu), n_channelcombis, 1);
    data.TEprepare.u_in_samples      = repmat(predicttimevec_u_samples(uu), n_channelcombis, 1);
    data.TEprepare.u_in_ms           = repmat(predicttimevec_u(uu), n_channelcombis, 1);
    
    % update fileidout to include information on u    
    cfgTESS.fileidout=strcat(fileidout,'_RAG4_TGA_u_',num2str(predicttimevec_u(uu)));
        
    % branch here for GPU calculation   
    data.TEprepare.cfg.verbosity = 'none';
    if strcmp(data.TEprepare.ensemblemethod,'yes')
        TGA_results{uu}=TEsurrogatestats_ensemble(cfgTESS,data);
    else
        TGA_results{uu}=TEsurrogatestats(cfgTESS,data);
    end
    data.TEprepare.cfg.verbosity = verbosity;
    
    if isfield(data,'groupprepare')
        TGA_results{uu}.groupprepare = data.groupprepare;
    end

end


[opt_u_vec TEmat]= TEfindmaxte(TGA_results);

data.TEprepare.u_in_ms = opt_u_vec;
data.TEprepare.u_in_samples = round(opt_u_vec/1000*data.fsample);
data.TEprepare.cfg.predicttime_u = opt_u_vec;
