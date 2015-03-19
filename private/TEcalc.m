function [te mi tel] = TEcalc(cfg,ncount);

% TECALC: Function splits up concatenated neighbour counts 
% returned by TEcallGPUsearch and calculates transfer entropy
% and mutual information for each chunk (i.e. for individual
% data and each surrogate data set for the current 
% channelpair). 'VW_ds' is used as estimator.
%
%
% * REFERENCE INFORMATION
%
%   - transfer entropy
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
% ncount = structure returned by TEcallGPUsearch, contains
%	   neighbour counts for all relevant embedding spaces
%	.ncount_p21_p2
%	.ncount_p21_p2
%	.ncount_p21_21 
%	.ncount_p21_2 
%	.ncount_12_1 
%	.ncount_12_2 
%	.nchunks       = total nr of chunks used in the neighbour searches
%			 (equal to numpermutation+1)
%
%   cfg            = configuration structure - it MUST contain:
%       .k_th      = number of neighbors for fixed mass search (controls
%                    balance of bias/statistical errors) (default = 4)
%       .MIcalc    = determines whether MI is calculated from neighbor
%		     counts
%
%
% * OUTPUT PARAMETERS
%
% te = vector of transfer entropy values, te(1) is the TE value for 
%      original data, te(2:end) are TE values for surrogate data
% mi = vector of mutual information values, mi(1) is the MI for 
%      original data, mi(2:end) are MI values for surrogate data
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
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral, Patricia
% Wollstadt
% Frankfurt 2013

% CHANGELOG
%
% 2014-24-13 PW: I made MI calculation optional (can be switched on/off
% using cfg.MIcalc
% 2014-10-11 PW: I fixed a bug in the assignment of NaNs to the MI results

% some auxiliary values to cut results from GPU search into chunks
nchunks   = ncount.nchunks;		% this is the total no. chunks over all calls to the GPU
chunksize = length(ncount.ncount_p21_p2)/nchunks;
cutpoint  = chunksize;


% data structures for results
te  = zeros(1,ncount.nchunks);
mi  = zeros(1,ncount.nchunks);
if cfg.TELcalc
    tel = zeros(ncount.nchunks,chunksize);
else
    tel = [];
end

for ii=1:ncount.nchunks
	
	%% Transfer entropy
	nc_p21_2  = single(ncount.ncount_p21_2(cutpoint-chunksize+1:cutpoint));
	nc_p21_p2 = single(ncount.ncount_p21_p2(cutpoint-chunksize+1:cutpoint));
	nc_p21_21 = single(ncount.ncount_p21_21(cutpoint-chunksize+1:cutpoint));
	
	te(ii)    = psi(cfg.kth_neighbors)+mean(psi(nc_p21_2+1)-psi(nc_p21_p2+1)-psi(nc_p21_21+1)); 
	if cfg.TELcalc
		tel(ii,:) = psi(cfg.kth_neighbors)+psi(nc_p21_2+1)-psi(nc_p21_p2+1)-psi(nc_p21_21+1); 
	end
    
	%% Mutual information
	if strcmp(cfg.extracond, 'Faes_Method')
		mi(ii) = NaN;			% MW: Quick hack until I have some more time 
    elseif cfg.MIcalc == 0
		mi(ii) = NaN;
	else
		nc_12_1 = single(ncount.ncount_12_1(cutpoint-chunksize+1:cutpoint));
		nc_12_2 = single(ncount.ncount_12_2(cutpoint-chunksize+1:cutpoint));
		L = size(nc_12_1,1);
		mi(ii) = psi(cfg.kth_neighbors)+psi(L)-mean(psi(nc_12_1+1)+psi(nc_12_2+1)); 
	end
	
	cutpoint = cutpoint + chunksize;
end
