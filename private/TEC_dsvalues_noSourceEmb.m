function [te,mi] = TEC_dsvalues_noSourceEmb(ts_1,ts_2,dim,tau,u,k_th,TheilerT, extracond,trialnum,channelpair)

% function [te,mi] = TEC_dsvalues_noSourceEmb(ts_1,ts_2,dim,tau,u,k_th,TheilerT, extracond,trialnum,channelpair)
%
% TRANSFERENTROPYVALUES computes the transfer entropy (TE) among a given
% pair of time series. source (ts_1) -> target (ts_2). The function doesn't
% use an embedding for the source.
%
% This function is called by the transferentropy.
%
% REFERENCE INFORMATION
%   - The concept of TE appears in Schreiber's article,
%     "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464 (2000).
%   - For the estimation of probability densities needed for the TE
%     computation, the function implements the Kraskov-Stoegbauer-Grassberger
%     estimator described in Kraskov et al. "Estimating mutual information",
%     Phys. Rev. E 69 (6) 066138, (2004).
%
% * DEPENDENCIES
%   - Package TSTOOL is used at nearest neighbors searches
%     required for the KSG estimator.
%
% INPUT PARAMETERS
%   - cfg       = configuration structure
%   - ts_1      = time series 1
%   - ts_2      = time series 2 (ts_2 should be of equal length than ts_1)
%   - dim       = embedding dimension
%   - tau       = embedding delay in number of sampled points AND points 
%                 ahead for the advance vector in number of sampled points 
%                 (this is different from TEvalues)
%   - u         = points ahead for the advance vector in number of sampled
%                 points (from ts_1 to ts_2 prediction point only,
%                 this is different from TEvalues)
%   - k_th      = number of neighbors for fixed mass search (controls 
%                 balance of bias/statistical errors)
%   - TheilerT  = number of temporal neighbors excluded to avoid serial 
%                 correlations (Theiler correction)
%   - extracond = string indicating any extra conditioning to be done
%
% developers want to have access to pinotsets created here, per trial and channelpair
% hence inputs also includes: 
%   - trialnum  = the current trial number
%   - channelpair = the current channel pair
%
% OUTPUT PARAMETERS
%   - te = transfer entropy time series 1 -> time series 2
%   - mi = mutual information
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
% Frankfurt 2012

%% Changelog
%
% 2014/11/11: PW, TSTool functions are now called "in bulk" for all points, 
%             i.e. calls to neighbour searches are no longer conducted for
%             each point separately

%% Developper note:
% this code takes u as an interaction delay from ts_1 to ts_2; 
% and tau as the interval for which to predict ts_2 
% DO NOT feed this function with the difference u-tau for the input u




%% Preprocessing of the data to be read by TSTOOL %%

% Z-scoring the time series
z_data_1 = zscore(ts_1);
z_data_2 = zscore(ts_2);

%% Creating the embedding vectors

% Computing effective lengths

T = length(ts_1);     % length of full time series
M = T-(dim-1)*tau;    % (only used in next line) number of points inside the time series available for delay embedding
L = M-max(u,tau)-1;   % number of points inside the time series available for advance and delay embedding
FirstPredictionP=T-L; % all referencing of embeddings is done wrt the prediction points
%WOI = 1:L;            % Indices of the prediction (and their corresponding emb. vectors)
                      % the first prediction is made for FirstPredictionP


% Initialization of embedding vectors
pointset_1 = zeros(L,1);
pointset_2 = zeros(L,dim);
pointset_p2 = zeros(L,dim+1);
pointset_21 = zeros(L,dim+1);
pointset_12 = zeros(L,dim+1);
pointset_p21 = zeros(L,dim+2);

% Embedding vectors names containing '2' indicate the 'target' time series

for ii = 1:L           % Marginal distributions
    pointset_1(ii) = z_data_1(ii+FirstPredictionP-u); 
    for jj = 1:dim
        % jump ahead to FirstPredictionP, then back by u (for ts_1) or tau
        % (for ts_2)        
        pointset_2(ii,jj) = z_data_2(ii+FirstPredictionP-1-(jj-1)*tau);
    end
end

for ii = 1:L           % Join distributions of the two time series
    
    pointset_12(ii,1) = z_data_1(ii+FirstPredictionP-u);
    
    for jj = 1:dim
        
            % vector with all the source embeddings first, then the targets
            % (do we really need this?)
            pointset_12(ii,jj+1) = z_data_2(ii+FirstPredictionP-1-(jj-1)*tau);
            % vector with all the target embeddings first, then the sources
            pointset_21(ii,jj) = z_data_2(ii+FirstPredictionP-1-(jj-1)*tau);
    end
    
    pointset_21(ii,dim+1) = z_data_1(ii+FirstPredictionP-u);
    
end

for ii = 1:L           % Join distributions of join marginal and own future states
    for jj = 1:2*dim+1
        if jj == 1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP);
        elseif jj > 1 && jj <= dim+1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP-1-(jj-2)*tau);
        else
            pointset_p21(ii,jj) = z_data_1(ii+FirstPredictionP-u-(jj-dim-2)*tau);
        end
    end
end

for ii = 1:L           % Join distributions of join marginal and own future states
    for jj = 1:dim+2
        if jj == 1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP);
        elseif jj > 1 && jj <= dim+1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP-1-(jj-2)*tau);
        else
            pointset_p21(ii,jj) = z_data_1(ii+FirstPredictionP-u);
        end
    end
end

% Augment conditioning time series '2' by the sample of the source 
% time series (1) that occurs  at the time of the prediction point in the
% target time series, i.e. the
% source equivalent of 'p'. This goes back to an idea of Luca Faes
% (reference?)
if  strcmp(extracond, 'Faes_Method')
 % augment the pointset (2) for the past of the target timeseries by the sample
 % of the source (1) at the prediction point (ii+FirstPredictionP)
 % disp('***** using faes method***')
 % do this for all elevant pointsets that contain (2)
 % do it in their after their last dimension i.e. at (end+1)
 %
 % preallocate additional memory
 pointset_2(:,end+1)   = zeros(L,1);
 pointset_21(:,end+1)  = zeros(L,1);
 pointset_p21(:,end+1) = zeros(L,1);
 pointset_p2(:,end+1)  = zeros(L,1);
 
 % fill the new last rows of the embedding state matrices with source values
 % at the prediction time (ii+FirstPredictionP)
 for ii = 1:L 
 pointset_2(ii,end) = z_data_1(ii+FirstPredictionP); % augment marginal target
 pointset_21(ii,end)= z_data_1(ii+FirstPredictionP); % augment joint of source and target
 pointset_p21(ii,end)= z_data_1(ii+FirstPredictionP); % augment joint of source and target and predictee
 pointset_p2(ii,end)= z_data_1(ii+FirstPredictionP); % augment marginal target plus predictee
 end
 % TODO: fix these changes for MI computation...
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for developers uncomment below
%save(strcat('/data/home1/pwollsta/CUDA/test_RS-data/Pointsets/pointset_sgncmb_',num2str(channelpair),'_trial_',num2str(trialnum),'.mat'), ...
%    'pointset_12','pointset_2','pointset_p2','pointset_p21','pointset_21')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nearest neighbors search (fixed mass)

% Preprocessing for nearest neighbor searches
% makes lookups for the points of interest faster
atria_1 = nn_prepare(pointset_1,'maximum');
atria_2 = nn_prepare(pointset_2,'maximum');
atria_p2 = nn_prepare(pointset_p2,'maximum');
atria_12 = nn_prepare(pointset_12,'maximum');
atria_21 = nn_prepare(pointset_21,'maximum');
atria_p21 = nn_prepare(pointset_p21,'maximum');

% Finding the k_th nearest neighbor
[index_p21, distance_p21] = nn_search(pointset_p21,atria_p21,1:L,k_th,TheilerT);
[index_12, distance_12] = nn_search(pointset_12,atria_12,1:L,k_th,TheilerT);
clear index_p21 index_12

%% Nearest neighbor search (fixed radius)

ncount_p21_p2 = range_search(pointset_p2,atria_p2,1:L,distance_p21(1:L,k_th)-eps,TheilerT);
ncount_p21_21 = range_search(pointset_21,atria_21,1:L,distance_p21(1:L,k_th)-eps,TheilerT);
ncount_p21_2  = range_search(pointset_2,atria_2,1:L,distance_p21(1:L,k_th)-eps,TheilerT);

ncount_12_1 = range_search(pointset_1,atria_1,1:L,distance_12(1:L,k_th)-eps,TheilerT);
ncount_12_2 = range_search(pointset_2,atria_2,1:L,distance_12(1:L,k_th)-eps,TheilerT);

%% Transfer entropy
te = psi(k_th)+mean(psi(ncount_p21_2+1)-psi(ncount_p21_p2+1)-psi(ncount_p21_21+1));

%% Mutual Information
mi = psi(k_th)+psi(L)-mean(psi(ncount_12_1+1)+psi(ncount_12_2+1));
% MW: Quick hack until I have some more time 
if  strcmp(extracond, 'Faes_Method')
    mi = NaN;
end

return;
