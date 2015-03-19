function pointset = TEembedding_noSourceEmb(ts_1,ts_2,dim,tau,u,extracond)

% TEEMBEDDING_NOSOURCEEMB: This function embeds the provided data and 
% returns the embedded pointsets without embedding the source time series. 
% The output is intendet for later use in GPU TE-calculation. Not 
% embedding the source time course may help to determine the interaction 
% delay u in some applications
% 
%
%
% * REFERENCE INFORMATION
%     - The concept of TE appears in Schreiber's article,
%       "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464
%       (2000).
%     - For the estimation of probability densities needed for the TE
%       computation, the function implements the Kraskov-Stoegbauer-
%       Grassberger estimator described in Kraskov et al. "Estimating
%       mutual information", Phys. Rev. E 69 (6) 066138, (2004).
%
%
% * INPUT PARAMETERS
%
%  ts_1 = source time series
%  ts_2 = target time series (should be of equal length as ts_1)
%  dim  = embedding dimension
%  tau  = embedding delay in number of sampled points AND points 
%         ahead for the advance vector in number of sampled points 
%  u    = points ahead for the advance vector in number of sampled
%         points (from ts_1 to ts_2 prediction point only)
%
%  extracond = perform conditioning in te tansfer entropy formula on
%                  additional variables. Values:
%                     'Faes_Method' - include the future sample of the 
%                     source at the prediction time into the state vector
%                     of the past of of the target to condition on it.
%                     In principle, this removes any volume conduction effect.
%                     (Faes L. et al, Phys Rev E, 2011)
%                     'Battaglia_Method' - NOT implemented yet. This
%                     methods conditions on the "mean activity" of the
%                     system, i.e. a global system state.This will be done
%                     by creating a channel carrying that signal which will
%                     then be used as an addional entry for the past state
%                     of the source 
%
% * OUTPUT PARAMETERS
%
%   pointsets   = embedded individual trials
%		  dimension: [datapoints embedding_dimension]
%       pointset_1   = source time series embedded
%	pointset_2   = target time series embedded
%	pointset_p2  = perdiction point + target time series embedded
%	pointset_21  = joint target and source time series embedded
%	pointset_12  = joint source and target time series embedded
%	pointset_p21 = perdiction point + joint target and source time series 
%		       embedded
%
%
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
% 
%
% Version 2.0 by Michael Lindner, Raul Vicente, Michael Wibral, Patricia
% Wollstadt
% Frankfurt 2013

% CHANGELOG


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

% build result structure and return
pointset.pointset_1   = pointset_1;
pointset.pointset_2   = pointset_2;
pointset.pointset_p2  = pointset_p2;
pointset.pointset_21  = pointset_21;
pointset.pointset_12  = pointset_12;
pointset.pointset_p21 = pointset_p21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for developers uncomment below
%save(strcat('/data/home1/pwollsta/CUDA/test_RS-data/Pointsets/pointset_sgncmb_',num2str(channelpair),'_trial_',num2str(trialnum),'.mat'), ...
%    'pointset_12','pointset_2','pointset_p2','pointset_p21','pointset_21')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%