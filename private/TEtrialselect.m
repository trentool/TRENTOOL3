function [trials,nrtrials]=TEtrialselect(cfg,datacell,ACT,channelcombi)

% TEtrialselect calculates the indices in the trialdimension of the data
% based on the alternatives ('ACT','range','all') given in the cfg.
%
% This function is called by the functions transferentropystats or
% transferentropy.
%
%
% * OUTPUT PARAMETERS
%   trials   = cell containing the indices of trials per channel used for
%              calculation of transfer entropy. (channelcombi x n)
%   nrtrials = matrix including the number of trials used per
%              channelcombination (channelcombi x 2).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2009
%
% CHANGELOG:
% 2011-12-28: ML changed the internal cells and matrices (datacell, datamat,
% ACT, trials, nrtrials) to a common structure (channelcombi x ??)

% note: channelcombi either refers to
% a ROW vector of all used channels
% when called from Ragwitz or Cao criteria calculation
% a nx2 matrix of channel index combinations
% when called from TEprepare outside of ragwitz or Cao
% (the cao method for embedding optimisation is deprecated)





%% Remember the working directory
working_directory = pwd;

% select trials per channel
% ---------------------------------------------------

% create zeros result matrics
nrtrials=zeros(size(channelcombi, 1),1);
trials=cell(size(channelcombi, 1), 2);

if strcmp(cfg.trialselect, 'ACT')
    % select trials by ACT criteria
    actthreshold = [0 cfg.actthrvalue];
    badchancmb=zeros(size(channelcombi, 1),1);


    for ii = 1:size(channelcombi, 1) % here channelcombi is a nx2 matrix of channel indices from all combinations
        for pp = 1:2
            % select trials from target channel reaching ACT criteria
            trialselect = find(ACT(ii,2,:)>actthreshold(1) & ACT(ii,2,:)<=actthreshold(2));
            nrtrials(ii,pp)=size(trialselect,1);
            if nrtrials(ii,pp)<cfg.minnrtrials              %detect bad channel combinations
                badchancmb(ii)=1;
            end
        trials{ii, pp}=trialselect;
        end
    end

    % error message if bad channel exist
    if sum(badchancmb) > 0
         fprintf(['\nbad channel combinations: ',num2str(sum(badchancmb)),'\n']);
%         for ii = 1 : size(badchancmb,1)
%             if badchancmb == 1
%                 fprintf(strcat(cell2mat(cfg.sgncmb(ii,:)),'\n'));
%             end
%         end
        fprintf('\n')
        error('TRENTOOL ERROR: less than mininmum nr of trials reached ACT threshold criteria - to solve the problem change the ACT threshold (cfg.actthrvalue) or remove the bad channel combinations');
    end

elseif strcmp(cfg.trialselect, 'range')
    % select trials from range

    for ii = 1:size(channelcombi, 1)
        for pp = 1:2
            trials{ii,pp} = cfg.trial_from:cfg.trial_to;
            nrtrials(ii,pp) = cfg.trial_to+1-cfg.trial_from;
        end
    end


elseif strcmp(cfg.trialselect, 'no')
    % use all trials
    for ii = 1:size(channelcombi, 1)
        for pp = 1:2
            trials{ii,pp} = 1:size(datacell{ii,pp},1);
            nrtrials(ii,pp) = size(datacell{ii,pp},1);
        end
    end

end;


%% Returning to the working directory
cd(working_directory)

end