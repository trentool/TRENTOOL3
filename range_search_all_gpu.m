% ncount = range_search_all_gpu(single(pointset),single(pointset),
%   single(radius),TheilerT,nchunks)
%
% RANGE_SEARCH_ALL_GPU performs a range search on the pointsets provided,
% returns the number of neighbours for each point in the pointset within
% the provided radius. All pointsets and the radius vector have to be 
% entered as SINGLE PRECISION values.
%
%
% *INPUT
%
%	pointset = set of individually embedded points, both pointsets
%		   have to be the same and have to be casted to single
%		   precision
%	radius   = vector of individual distances for every point in the
%		   pointset (has to have as many entries as there are
%		   points in the pointset)
%       TheilerT = number of temporal neighbors excluded to avoid serial
%                  correlations (Theiler correction) (default = ACT)
%	nchunks  = number of chunks that are entered into the function
%		   (one chunk is equal to all points for a channel
%		   combinations, trials are ignored)
%
%
% *OUTPUT
%
%	ncount   = vector containing the number of neighbours for each
%		   point
%
%
% M. Zarzuela, 2013