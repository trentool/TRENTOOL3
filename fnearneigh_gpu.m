% distance = fnearneigh_gpu(single(pointset),single(pointset),
%	k,TheilerT,nchunks);
%
% FNEARNEIGH_GPU performs a k nearest neighbour search on the GPU. All
% pointsets have to be entered as SINGLE PRECISION values.
%
% *INPUT
%
%	pointset = set of individually embedded points, both pointsets
%		   have to be the same and have to be casted to single
%		   precision
%       k        = number of neighbors for fixed mass search (controls
%                  balance of bias/statistical errors) (default = 4)
%       TheilerT = number of temporal neighbors excluded to avoid serial
%                  correlations (Theiler correction) (default = ACT)
%	nchunks  = number of chunks that are entered into the function
%		   (one chunk is equal to all points for a channel
%		   combinations, trials are ignored) 
%
%
% *OUTPUT
%
%	distance = matrix with distances to the first k neighbours 
%		   (distance(:,k) is later used for range searches)
%
%
% M. Zarzuela, 2013