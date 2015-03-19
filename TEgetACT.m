function ACT = TEgetACT(cfg, data)

% TEGETACT returns the autocorrelation time (ACT) for all channel
% combinations specified in the cfg structure. 
%
% INPUT
%   cfg         configuration structure 
%       .maxlag = maximum lag for the calculation of the ACT (in samples, 
%                 default=1000)
%       .toi    = time of interest (same units as the time vectors in the
%                 data structure, usually seconds, 1x2 vector with
%                 beginning and end of the analysis window)                 
%       .channel OR .sgncmb = list of channels or signal combinations to be
%                             analyzed, if a list of channels is provided,
%                             the function builds all possible channel
%                             combinations from the list.
%       
%   data        Fieldtrip raw data structure - it MUST contain:
%       .trial    = cell array (nr of channels x nr of samples) containing
%                   the data for each trial
%       .time     = cell (1xnr of samples) containing the time indices for
%                   each trial (in seconds)
%       .label    = cell (1xnr of channels), containing the labels
%                   (strings) of channels included in the data
%       .fsample  = value of sampling rate (in Hertz)
%
%
% OUTPUT
%   ACT     = matrix with ACT values in samples (channel combinations x 2 x
%             no. trials, where the two columns in the second dimension 
%             contain ACT values for the source and target process of the
%             respective channel combination)
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Patricia Wollstadt, Michael Wibral, Michael Lindner, Raul
% Vicente
% Frankfurt 2014

% CHANGELOG
%

%% check input

if ~isfield(cfg, 'maxlag'), cfg.maxlag = 1000;                        end;
if ~isfield(cfg, 'toi')   , error('TRENTOOL: Please specify a toi!'); end;


%% get channel combinations

% check if channel or channelcombinations are defined
if ~isfield(cfg, 'channel') && ~isfield(cfg, 'sgncmb') ,
    fprintf('\n')
    error('TRENTOOL ERROR: specify cfg.channel OR cfg.sgncmb, see help!');
elseif isfield(cfg, 'channel') && isfield(cfg, 'sgncmb') ,
    fprintf('\n')
    error('TRENTOOL ERROR: specify cfg.channel OR cfg.sgncmb, see help!');
elseif isfield(cfg, 'channel') && ~isfield(cfg, 'sgncmb') ,
    if size(cfg.channel,2)>size(cfg.channel,1)
        cfg.channel=cfg.channel';
    end
    channelselect = 1;
    % a warning because of some issue if only a subselection of
    % channels enters the analysis
    if max(size(cfg.channel))<size(data.trial{1},1) % If there are less channels
        fprintf('\nTRENTOOL WARNING: your are specifying a subselection of channels \n - please use cfg.sgncmb to specify channelcombinations directly'); 
    end
elseif ~isfield(cfg, 'channel') && isfield(cfg, 'sgncmb') ,
    if size(cfg.sgncmb) ~= 2
        fprintf('\n')
        error('TRENTOOL ERROR: cfg.sgncmb has wrong dimensions, see help!');
    end
    channelselect = 2;
end;

[channelcombi,channelcombilabel] = TEchannelselect(cfg, data, channelselect);

%% prepare data cell

timeindices=zeros(1,2);
for ii = 1:size(cfg.toi,2)
    [col]=nearest(data.time{1}, cfg.toi(ii));
    timeindices(ii)=col;
end

datacell = cell(size(channelcombi,1),2);
for cc = 1:size(channelcombi,1)
    for pp = 1:2
        datamat = zeros(size(data.trial,2),size(data.trial{1},2));
        for ii = 1:size(data.trial,2)
            datamat(ii,:)=data.trial{ii}(channelcombi(cc,pp),:);
        end
        datacell{cc,pp}=datamat;
        clear datamat;
    end
end

%% calculate ACT

ACT = TEactdetect(datacell,cfg.maxlag,timeindices);