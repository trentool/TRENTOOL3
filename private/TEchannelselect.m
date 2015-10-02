function  [channelcombi,channelcombilabel]= TEchannelselect(cfg, data, channelselect)

% TEchannelselect calculates the indices of the used channels in
% in the channeldimension in the data based on the list of channels or 
% channelcombinations given in the
% cfg.
%
% This function is called by the functions transferentropystats or
% transferentropy.
%
% OUTPUT:   channelcombi      = nx2 matrix of indices of the channel
%                               combinations.
%           channelcombilabel = nx2 cell array of labels  of the channel
%                               combinations.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 2.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2012
%

%% Remember the working directory
working_directory = pwd;

% Select Channels
% -----------------------------------------------

%read channel labels from data
if size(data.label,1)==1
    data.label = data.label';
end

allchannels=data.label;


if channelselect == 1   % given channel list
    usedchannels=zeros(size(cfg.channel,1));
    channelcmp=cfg.channel;
    usedcounter = 0;
    for jj=1:size(allchannels,1)
        for kk=1:size(channelcmp,1)
            if strcmp(allchannels{jj},channelcmp{kk})
                usedchannels(kk)=jj;
                usedcounter = usedcounter + 1;
            end
        end
    end
    
%     channelcombi = nan((size(usedchannels,1)^2-size(usedchannels,1)),2);
%     channelcombilabel = nan((size(usedchannels,1)^2-size(usedchannels,1)),2);
    
    combicount = 1;
    for ii = 1:size(usedchannels,1)
        for jj = 1:size(usedchannels,1)
            if ii~=jj
                channelcombi(combicount,:)=[usedchannels(ii), usedchannels(jj)];
                channelcombilabel(combicount,:)=[data.label(usedchannels(ii)), data.label(usedchannels(jj))];
                combicount = combicount + 1;
            end
        end
    end    
    
    if usedcounter ~= size(usedchannels,1)
        fprintf('\n')
        error('TRENTOOL ERROR: mismatch between cfg.channel and data.label - check for typos');
    end
elseif channelselect == 2   % given channel combinations
    usedchannels=zeros(size(cfg.sgncmb)); % naming not optimal
    for cc = 1:2
        channelcmp=cfg.sgncmb(:,cc);
        usedcounter = 0;
        for jj=1:size(allchannels,1)
            for kk=1:size(channelcmp,1)
                if strcmp(allchannels{jj},channelcmp{kk})
                    usedchannels(kk,cc)=jj;
                    usedcounter = usedcounter + 1;
                end
            end
        end
        channelcombilabel = cfg.sgncmb;
        if usedcounter ~= size(usedchannels,1)
            fprintf('\n')
            error('TRENTOOL ERROR: mismatch between cfg.sgncmb and data.label - check for typos');
        end
    end
    channelcombi=usedchannels;
end


%% Returning to the working directory
cd(working_directory)

