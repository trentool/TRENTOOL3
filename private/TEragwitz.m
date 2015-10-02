function [mre] = TEragwitz(cfg,timeSeries,pPoints,u,flagNei,sizeNei,dimEmb,tauEmb,Theiler,dimMax,tauMax)

%   This function returns a nonlinear predictor based on the method of
%   analogues while taking into account the possible manifold topology of the data
%
%   INPUTS:
%   cfg                 -->  The cfg structure given to the caller function
%                            of this function
%   timeSeries          -->  Time series to predict
%   trainingPoints      -->  Points for the training set
%   pPoints             -->  Number of points we want to predict (it has to
%                            be smaller than length(timeSeries)-(dimEmb-1)*tauEmb-u-1)
%   u                   -->  Prediction horizon
%   flagNei             -->  Radius or Mass nearest neighbor search
%   sizeNei             -->  Radius or Mass for the neighbor search according to flagNei
%   dimEmb              -->  Dimensionality for the Takens embedding
%   tauEmb              -->  Delay for the Takens embedding
%   Theiler             -->  Theiler correction excludes from the nearest neighbors search the Theiler samples closer in time
%   dimMax              -->  maximum embedding dimension
%   tauMax              -->  maximum tau (in samples: max of cfg.tau * max
%                            of ACT)
%
%   OUTPUTS:
%   mre                 -->  Mean relative error
%
%
%   References:
%   (1) "Markov models from data by simple nonlinear time series
%   predictors in delay embedding spaces", by Ragwitz and Kantz, Physical
%   Review E, vol 65, 056201 (2002).
%   (2) "A Global Geometric Framework for Nonlinear Dimensionality
%   Reduction", by Joshua B. Tenenbaum, Vin de Silva, John C. Langford,
%   Science, vol 290, pp 2319-2322, (2000).
%
%
%   by Juhan Aru and Raul Vicente 07/10/2010, Frankfurt am Main
%   latest code revision 10/19/10 by Michael Lindner, Frankfurt am Main


%% Preprocessing

% Z-scoring the time series
data = zscore(timeSeries);
% data = (timeSeries-mean(timeSeries))/std(timeSeries);  %changed by vp to avoid the use of the stat toolbox


%% Takens embedding

% Computing effective lengths
T = length(data);           % length of full time series
L = T-(dimMax-1)*tauMax-u;  % number of points used from the time series ready for delay embedding


% Initialization of the embedding vector
pointset = zeros(L,dimEmb);

% Embedding of the training points


for ii = 1:L
    for jj = 1:dimEmb
        pointset(ii,jj) = data(ii+(dimMax-1)*tauMax-(jj-1)*tauEmb);
    end
end

%% Nearest neighboor search and local manifold predictor

% Change to TSTOOL directory for the search of neighbors
%working_dir = pwd;
%[dir_mex] = TEarch(cfg);
%cd(dir_mex);
TEarch;

% Preprocessing for nearest neighbor searches
atria = nn_prepare(pointset,'maximum');

%%     Local manifold predictor

lp = zeros(1,pPoints);
actual = zeros(1,pPoints);

switch flagNei
    
    case 'Range'
        
        % Finding the nearest neighbor of the predictees in the range of sizeNei
        [count, neighbors] = range_search(pointset,atria,1:pPoints,sizeNei-eps,Theiler);
        
        % Finding the index of the forward time projections of the neighbors
        % of the set of predictees
        for pp = 1:pPoints
            
            % Local predictor - these predcitions are u+(dimMax-1)*tauMax
            % times steps ahead for all compared embedding so that the
            % differences in prediction really come from the embedding and
            % not from chnages in prediction horizon
            lp(pp) = sum(timeSeries(neighbors{pp,1}+u+(dimMax-1)*tauMax))./count(pp);
            
            % Actual value - these predcitions are u+(dimMax-1)*tauMax
            % times steps ahead -> see explanation above
            actual(pp) = timeSeries(pp+(dimMax-1)*tauMax+u);
            
        end
        
    case 'Mass'
        
%         % MW debug:
%         disp('pPoints')
%         disp(pPoints)
%         disp('L')
%         disp(L)
%         % end MW debug
%         
        % Finding the nearest neighbor of the predictees in the range of sizeNei
        [index, distance] = nn_search(pointset,atria,1:pPoints,sizeNei,Theiler);
        
        % Finding the index of the forward time projections of the neighbors
        % of the set of predictees
        for pp = 1:pPoints
            
            % Local predictor - see note  above for prediction horizon
            lp(pp) = sum(timeSeries(index(pp,:)+u+(dimMax-1)*tauMax))./sizeNei;
            
            % Actual value - see note  above for prediction horizon
            actual(pp) = timeSeries(pp+(dimMax-1)*tauMax+u);
            
        end
        
    otherwise
        
        disp('Please specify in the arguments either Range or Mass search')
        
end

mre = (sum((lp-actual).^2)/pPoints)/std(timeSeries);

%cd(working_dir)