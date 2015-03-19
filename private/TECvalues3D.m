function [te,mi] = TECvalues3D(ts_1,Emb_1,ts_2,Emb_2,tau,u,k_th,TheilerT)

% TECVALUES computes the transfer entropy (TE) among a given
% pair of time series. source (ts_1) -> target (ts_2)
% This function is based on a modifued predictor where the autoprediction
% of the target future is tau points ahead and the influence of the source
% is checked u points ahead
%
% This function is called by transferentropy.m.
%
% The function has been modified to account for the problem of decreasing
% autoprediction for large prediction times (compare to TEvalues). This has
% been fixed by simply embedding such that the prediction time is actually
% tau but the source signal is shifted back k units of time. In the calling
% code the value of cfg.predictime_u should be passed as the value for k
%
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
%   - tau       = embedding delay in number of sampled points
%   - u         = points ahead for the advance vector for the source signal
%                 in number of sampled points in addition to the general 
%                 advance by tau
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
% Bonn 2012


%% Preprocessing of the data to be read by TSTOOL %%

% Z-scoring the time series
z_data_1 = zscore(ts_1);
z_data_2 = zscore(ts_2);
z_Embdata_1 = zscore(Emb_1);
z_Embdata_2 = zscore(Emb_2);

dim = size(Emb_1,1);
%% Creating the embedding vectors

% Computing effective lengths

T = length(ts_1);         % length of full time series
M = T-tau-u;              % number of points inside the time series ready for delay embedding
L = M-u;                  % number of points inside the time series ready for advance and delay embedding
WOI = 1:L;                % Window of interest


% Initialization of embedding vectors
pointset_1 = zeros(L,dim+1);
pointset_2 = zeros(L,dim+1);
pointset_p2 = zeros(L,dim+2);
pointset_21 = zeros(L,2*(dim+1));
pointset_12 = zeros(L,2*(dim+1));
pointset_p21 = zeros(L,(2*dim+1)+1); %???

% Embedding vectors
% Marginal distributions
for ii = 1:L           % loop over time samples   
    pointset_1(ii,1) = z_data_1(ii); % source
    pointset_2(ii,1) = z_data_2(ii+u); % target
    
    for jj = 1:dim     % loop over embedding dimensions
        pointset_1(ii,jj+1) = z_Embdata_1(jj,ii); % source
        pointset_2(ii,jj+1) = z_Embdata_2(jj,ii+u); % target
    end
end

% Joint distributions of marginal and own future state
for ii = 1:L          % loop over time samples     
    pointset_p2(ii,1) = z_data_2(ii+u+tau); % target
    pointset_p2(ii,2) = z_data_2(ii+u); % target
    for jj = 1:dim     % loop over embedding dimensions
        pointset_p2(ii,jj+2) = z_Embdata_2(jj,ii+u); % target
    end
    
end


pointset_12(:,1:dim+1)=pointset_1;
pointset_12(:,dim+2:end)=pointset_2;
pointset_21(:,1:dim+1)=pointset_2;
pointset_21(:,dim+2:end)=pointset_1;

% Joint distributions of the two time series - no future state
% for ii = 1:L           
%     for jj = 1:2*dim
%         if jj <= dim                      
%             pointset_12(ii,jj) = z_data_1(ii+(dim-1)*tau-(jj-1)*tau); % source
%             pointset_21(ii,jj) = z_data_2(ii+(dim-1)*tau-(jj-1)*tau+k); % target
%         else
%             pointset_12(ii,jj) = z_data_2(ii+(dim-1)*tau-(jj-dim-1)*tau+k); % target
%             pointset_21(ii,jj) = z_data_1(ii+(dim-1)*tau-(jj-dim-1)*tau); % source
%         end
%     end
% end

% Joint distributions of joint marginal and own future states
for ii = 1:L      % loop over time samples     
    
    pointset_p21(ii,1)=z_data_2(ii+tau);
    pointset_p21(ii,1:end)= pointset_21;
    
end
% for ii = 1:L      % loop over time samples     
%     for jj = 1:2*dim+1 % loop over joint dimensions
%         if jj == 1
%             pointset_p21(ii,jj) = z_data_2(ii+(dim-1)*tau+tau); % target - tau used as prediction time - OK
%         elseif jj > 1 && jj <= dim+1
%             pointset_p21(ii,jj) = z_data_2(ii+(dim-1)*tau-(jj-2)*tau+k); % target (+k ??? is this correct)
%         else
%             pointset_p21(ii,jj) = z_data_1(ii+(dim-1)*tau-(jj-dim-2)*tau); % source
%         end
%     end
% end


%% Nearest neighbors search (fixed mass)




% Preprocessing for nearest neighbor searches
atria_1 = nn_prepare(pointset_1,'maximum');
atria_2 = nn_prepare(pointset_2,'maximum');
atria_p2 = nn_prepare(pointset_p2,'maximum');
atria_12 = nn_prepare(pointset_12,'maximum');
atria_21 = nn_prepare(pointset_21,'maximum');
atria_p21 = nn_prepare(pointset_p21,'maximum');

% Finding the k_th nearest neighbor
% [index_2, distance_2] = nn_search(pointset_2,atria_2,WOI,k_th,TheilerT);
% [index_p2, distance_p2] = nn_search(pointset_p2,atria_p2,WOI,k_th,TheilerT);
% [index_21, distance_21] = nn_search(pointset_21,atria_21,WOI,k_th,TheilerT);
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
