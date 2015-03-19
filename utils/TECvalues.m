function [te,mi] = TECvalues(ts_1,ts_2,dim,tau,u,k_th,TheilerT)

% TRANSFERENTROPYVALUES computes the transfer entropy (TE) among a given
% pair of time series. source (ts_1) -> target (ts_2)
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
WOI = 1:L;            % Indices of the prediction (and their corresponding emb. vectors)
                      % the first prediction is made for FirstPredictionP


% Initialization of embedding vectors
pointset_1 = zeros(L,dim);
pointset_2 = zeros(L,dim);
pointset_p2 = zeros(L,dim+1);
pointset_21 = zeros(L,2*dim);
pointset_12 = zeros(L,2*dim);
pointset_p21 = zeros(L,2*dim+1);

% Embedding vectors names containing '2' indicate the 'target' time series

for ii = 1:L           % Marginal distributions
    for jj = 1:dim
        % jump ahead to FirstPredictionP, then back by u (for ts_1) or tau
        % (for ts_2)
        pointset_1(ii,jj) = z_data_1(ii+FirstPredictionP-u-(jj-1)*tau); 
        pointset_2(ii,jj) = z_data_2(ii+FirstPredictionP-tau-(jj-1)*tau);
    end
end


for ii = 1:L           % Join distributions of marginal and own future state
    for jj = 1:dim+1
        if jj == 1 % take the prediction point
            pointset_p2(ii,jj) = z_data_2(ii+FirstPredictionP);
        else % take the prediction point, go back by tau, then embedd,
             % but take care of the offset in jj
            pointset_p2(ii,jj) = z_data_2(ii+FirstPredictionP-tau-(jj-2)*tau);
        end
    end
end

for ii = 1:L           % Join distributions of the two time series
    for jj = 1:2*dim
        if jj <= dim
            % vector with all the source embeddings first, then the targets
            % (do we really need this?)
            pointset_12(ii,jj) = z_data_1(ii+FirstPredictionP-u-(jj-1)*tau );
            % vector with all the target embeddings first, then the sources
            pointset_21(ii,jj) = z_data_2(ii+FirstPredictionP-tau-(jj-1)*tau);
        else
            pointset_12(ii,jj) = z_data_2(ii+FirstPredictionP-tau-(jj-dim-1)*tau);
            pointset_21(ii,jj) = z_data_1(ii+FirstPredictionP-u-(jj-dim-1)*tau );
        end
    end
end

for ii = 1:L           % Join distributions of join marginal and own future states
    for jj = 1:2*dim+1
        if jj == 1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP);
        elseif jj > 1 && jj <= dim+1
            pointset_p21(ii,jj) = z_data_2(ii+FirstPredictionP-tau-(jj-2)*tau);
        else
            pointset_p21(ii,jj) = z_data_1(ii+FirstPredictionP-u-(jj-dim-2)*tau);
        end
    end
end


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
[index_p21, distance_p21] = nn_search(pointset_p21,atria_p21,WOI,k_th,TheilerT);
[index_12, distance_12] = nn_search(pointset_12,atria_12,WOI,k_th,TheilerT);

%% Nearest neighbor search (fixed radius)

ncount_p21_p2 = zeros(L,1);
ncount_p21_21 = zeros(L,1);
ncount_p21_2 = zeros(L,1);
ncount_12_1 = zeros(L,1);
ncount_12_2 = zeros(L,1);

for i=1:L
    [count_p21_p2, neighbors_p21_p2] = range_search(pointset_p2,atria_p2,i,distance_p21(i,k_th)-eps,TheilerT);
    [count_p21_21, neighbors_p21_21] = range_search(pointset_21,atria_21,i,distance_p21(i,k_th)-eps,TheilerT);
    [count_p21_2,  neighbors_p21_2]  = range_search(pointset_2,atria_2,i,distance_p21(i,k_th)-eps,TheilerT);
    ncount_p21_p2(i) = count_p21_p2;
    ncount_p21_21(i) = count_p21_21;
    ncount_p21_2(i)  = count_p21_2;
end


for i=1:L
    [count_12_1, neighbors_12_1] = range_search(pointset_1,atria_1,i,distance_12(i,k_th)-eps,TheilerT);
    [count_12_2, neighbors_12_2] = range_search(pointset_2,atria_2,i,distance_12(i,k_th)-eps,TheilerT);
    ncount_12_1(i) = count_12_1;
    ncount_12_2(i) = count_12_2;
end


%% Transfer entropy
te = psi(k_th)+mean(psi(ncount_p21_2+1)-psi(ncount_p21_p2+1)-psi(ncount_p21_21+1));

%% Mutual Information
mi = psi(k_th)+psi(L)-mean(psi(ncount_12_1+1)+psi(ncount_12_2+1));


return;
