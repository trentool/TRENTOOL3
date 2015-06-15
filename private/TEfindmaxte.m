function opt_u  = TEfindmaxte(data)

% FUNCTION TEFINDMAXTE
% Finds the optimal delay by searching for the maximum TE value over various 
% assumed delays u.
% The function is called from within TEfindDelay. The output is a vector
% with optimal values for u for each channel combination.
%
% The optimal delay is found by looking for the smalles u, which maximes 
% the TE raw value (mean TE value over trials). 
%
%
% You can call this function directly as follows:
%         opt_u = TEfindmaxte(data)
%
% * DEPENDENCIES
%
% * INPUT PARAMETERS
%
%   data = Cell array with as many  cells as assumed delays u's that were 
%          scanned within TEfindDelay. Each cell contains the output of a 
%          call to TEsurrogatestats (TEpermtest and TEresult)
%          (see below for details)
%
%             TEpermtest
%            .TEpermvalues  = matrix with size:
%                             (channelpair,value)
%                           The last dimension "value" includes:
%                           1 - p_values of the statistic within the
%                               distribution given by the permutations
%                           2 - 1 (0), if the statistics is significant at
%                               the prescribed alpha level (or not)
%                           3 - 1 (0), if the statistics is significant
%                               after correction for multiple comparisons
%                               (or not)
%                           4 - 1 (0), mean difference or tvalue of mean
%                               difference depending on cfg.permstatstype
%                           5 - 1 (0), if instantaneous mixing (volume
%                               conduction) exists (or not)
%            .TEmat         = matrix containing raw TE values for each
%                             channel combination (channelcombi x trial for
%                             TE estimation on CPU or channelcombi x 1 for
%                             TE estimation using the ensemble method/GPU
%                             estimation)
%            .dimord        = dimensions of TEpermvalues
%            .cfg           = configuration file used to calculate TE and
%                             permtest
%            .sgncmb        = labels of channel combinations (source ->
%                             target)
%            .numpermutation = number of permutations
%            .ACT           = structure including
%                .act       = ACT matrix (channelcombi x 2 x trial)
%            .nr2cmc        = number of tests to correct for multiple
%                             comparisons
%            .TEprepare     = results of the function TEprepare from the
%                             data
%
%            AND
%
%           TEresult             (= Output structure of the function tranferentropy)
%          .TEmat       = resultmatrix including transfer entropy(TE)
%                         values (channelpairs x u x trial)
%          .MImat       = resultmatrix including mutual information (MI)
%                         values (channelpairs x u x trial)
%          .dimord      = 'channelpair_u_trial'; the dimensions of TEmat
%                         and MImat
%          .cfg         = configuration file used to calculate TE
%          .trials      = trial numbers selected from raw dataset
%          .act         = ACT matrix (channelcombi x 2 x trial)
%          .sgncmb      = labels of channel combinations (source -> target)
%          .TEprepare   = results of the function TEprepare from the
%                         data
%   if instantaneous mixing is found in the data, then another field will
%   be added:
%          .instantaneousmixing = matrix (channel x u) which indicates were
%          the instantaneous mixings were found (1) or not (0).%
%
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
% version2.0
% (C) Michael Wibral, Raul Vicente and Michael Lindner 2012

% CHANGELOG
% Nicu Pampu:bugfix for the nan index (initial was without value when using function max) now set to 1 
% 13-11-03 PW: add u-reconstruction for ensemble method
% 13-27-06 PW: field .groupprepare is added to the ouptput if
% TEprepare_group was run on the data previous to TE analysis
% 13-12-11 PW: field .TEmat containing the raw TE values for each channel
% combination is added to the ouptput, raw TE values are needed for the
% group statistics
% 27-11-14 PW: added option 'max_TE' for opt u selection
% 12-06-15 PW: changed the function to be called from within TEfindDelay

n_u            = length(data);           % no. scanned delays 
n_channelcombi = size(data{1}.sgncmb,1); % no. signal combinations

TGA = data{1};                                   % generate output structure
TGA.TEpermvaluesTmp = nan([n_channelcombi n_u]); % container to collect individual TEpermvalues
TGA.TEpermvalues = zeros([n_channelcombi 5]);    % preallocate the output array
TGA.sgncmb=data{1}.sgncmb;

% collect TEpermvalues tables for all u's into one temporary data structure
TERawmat  = nan(n_channelcombi,n_u);
uvec      = nan(n_u,1);
for uu=1:n_u
    uvec(uu)       = data{uu}.cfg.predicttime_u(1);
    TERawmat(:,uu) = mean(data{uu}.TEmat,2);
end

% loop over channel combinations
TEmat = nan(size(data{1}.TEmat)); % remember the raw TE values, needed for group statistics
opt_u = nan(n_channelcombi,1);
for cc=1:n_channelcombi
         
         % find the index of the optimal u, 'max' returns the first maximum 
         [maxTE, ind_maxTE] = max(TERawmat(cc,:));
                  
         opt_u(cc)   = data{ind_maxTE}.cfg.predicttime_u(1);  % get the u
         TEmat(cc,:) = data{ind_maxTE}.TEmat(cc,:);           % get raw TE values

end

% add u vector, remove volume conduction (volume conduction will result in  u = 0 ms).
if ~strcmp(data{1}.cfg.ensemblemethod,'yes')    
    VolCondIndicator = (ones(size(squeeze(TGA.TEpermvalues(:,5))))-squeeze(TGA.TEpermvalues(:,5))); % 0 for volume conduction, 1 for OK
    opt_u            = (opt_u).*VolCondIndicator;
end;

