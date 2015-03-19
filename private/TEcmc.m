function [significance correctm] = TEcmc(data, correctm, alpha, nrinstmix)

% TEcmc: This function implements the correction for multiple comparisons
% with false discovery rate or the more conservative Bonferroni correction.
%
% This function is called by the functions TEperm
%
% * REFERENCE INFORMATION
%   Genovese, C.R., Lazar, N.A., & Nichols, T. (2002). Thresholding of
%   statistical maps in functional neuroimaging using the false discovery
%   rate. Neuroimage, 15(4), 870-878.
%
%
% * OUTPUT PARAMETERS
%   significance = matrix including significances (1) after correction for
%   multiple comparison
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2010
%

% CHANGELOG
% 03/13/2014: PW changed the check of the no. elements in the data, this
% was not effective before (changed cfg.correctm instead of correctm)
% 01/19/2015: PW fixed a bug in the calculation of the individual
% thresholds (missing brackets)

% if number of multiple comparisons are smaller than 10 Bonferroni correction
% will be used automatically
if numel(data) <=10 && strcmp(correctm, 'FDR');
    %cfg.correctm = 'BONF';             % PW this doesn't lead to the use of BONF below, because it changes cfg.correctm not correctm
    correctm = 'BONF';
    fprintf('\nTRENTOOL WARNING: Number of Data to small for FDR -> Bonf was used instead');
end



if strcmp(correctm, 'FDR')
    
    % reshape and sort data
    dim = size(data);
    nrdata = numel(data);
    data = reshape(data, 1, nrdata);
    [sorteddata, index] = sort(data);
    
    % remove instanteneous mixing
    nrdatacor = nrdata-nrinstmix;
    sorteddatacor = sorteddata(1:end-nrinstmix);
    
    % calculate threshold (exact or by approximating the harmonic sum)
    if nrdata < 1000
        thresh = ((1:nrdatacor)/nrdatacor)  * alpha / sum(1.0./[1:nrdatacor]);
    else
        thresh = ((1:nrdatacor)/nrdatacor)  * alpha / (log(nrdatacor) + -psi(1)); % -psi(1) returns the Euler constant (gamma)       
    end
    
    % compare data to threshold and prepare output 
    significancecor = sorteddatacor <= thresh;
    significance = zeros(1,nrdata);
    significance(1:end-nrinstmix) = significancecor;   
    [dummy, unsorteddata] = sort(index);
    significance = significance(unsorteddata);        
    significance = reshape(significance, dim);
    
elseif strcmp(correctm, 'BONF')
    
    nrofcomp = size(data,1) * size(data,2) -nrinstmix; %* size(TEresult.TEmat,3);    
    significance = (data<=(alpha/nrofcomp));
    
end