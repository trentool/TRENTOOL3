function [ACT]= TEactdetect(datacell, maxlag, timeindices)

% TEactdetect estimates the auto-correlation time (ACT) over the
% lag range [-maxlags:maxlags] of the target channel of the datacell
%
% This function is called by the function transferentropystats and
% transferentropy.
%
% OUTPUT: ACT matrix (channelcombi x 2 x trial)
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
% CHANGELOG:
% 2011-12-28: ML changed the internal cells and matrices (datacell, datamat,
% ACT, trials, nrtrials) to a common structure (channelcombi x ??)


%% Remember the working directory
working_directory = pwd;

ACT=zeros(size(datacell,1),size(datacell,2),size(datacell{1,1},1));


thresh = exp(-1);
counter = 1;

fprintf('\t');
for cc = 1:size(datacell,1) % the number of channel combinations
    for pp = 1:size(datacell,2)  %%%% ML: TODO; MW note: number of channels in a combination (??)
        for trial = 1:size(datacell{cc,pp},1) % the number of trials in a particular combination
            TEwait(counter);            
            % correct data for mean
            data_cut = squeeze(datacell{cc,pp}(trial,timeindices(1):timeindices(2)));
            data_corr = data_cut-mean(data_cut);
            
            % calculate ACT
            c = TEautocorr(data_corr',maxlag);
%             % start debug
%             pause on
%             plot(c)
%             pause
%             % end debug
            d = c(maxlag+1:end);
            auxlag = 0:maxlag;
            % size(min(auxlag(find(d<thresh)))) % debugging
            if ~isempty(min(auxlag(find(d<thresh))));
                ACT(cc,pp,trial) = min(auxlag(find(d<thresh))); % ACT in samples (!)
            elseif isempty(min(auxlag(find(d<thresh))));
                %warning('ACT not determined using cfg.maxlag vaule. Please increase.')
                ACT(cc,pp,trial) = inf; %undetectable ACT is set to infitiy -> ACT threshhold in function transferentropy ACT<threshold
            end
            counter = counter + 1;
        end
    end
end

%% Returning to the working directory
cd(working_directory)



function [c]=TEautocorr(data,maxlag)
% calculating the autocorrelation of the data

% size of data
[M,N] = size(data);

% find nearest power
m = 2*M-1;
[F,E] = log2(abs(m));

% Check if m is an exact power of 2.
if ~isempty(F) && F == 0.5
    E = E-1;
end

% check infinity
checkinfinity = ~isfinite(F);
E(checkinfinity) = F(checkinfinity);


% Compute autocorrelation via FFT
fft_data = fft(data,2^E);
c = ifft(abs(fft_data).^2);


% Force real data
c = real(c);

% define lags for moving
lags = -maxlag:maxlag;
if maxlag >= M
    c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
else
    c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
end

% normalize data
c = c./c(maxlag+1);

