function [significance, correctm, nr2cmc] = TEcmc(data, correctm, alpha, nrinstmix)

% TEcmc: This function implements the correction for multiple comparisons
% with false discovery rate or the more conservative Bonferroni correction.
%
% This function is called by the functions TEperm
%
% * REFERENCE INFORMATION
%   Genovese, C.R., Lazar, N.A., & T. Nichols (2002). Thresholding of
%   statistical maps in functional neuroimaging using the false discovery
%   rate. Neuroimage, 15(4), 870-878.
%   Y. Benjamini & Y Hochberg (1995). Controlling the False Discovery Rate:
%   A Practical and Powerful Approach to Multiple Testing. J R Stat Soc B,
%   57(1), 289-300.
%
%
% * OUTPUT PARAMETERS
%   significance = matrix including significances (1) after correction for
%                  multiple comparison
%   correctm     = method used for correction for multiple comparisons
%   nr2cmc       = number multiple comparisons (e.g. used to correct the
%                  alpha level in Bonferroni correction)
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
% 01/19/2015: PW fixed a bug in the FDR algorithm, thresholds and P-Values
% have to be compared iteratively, starting with the smallest p
% 10/01/2015: PW nr2cmc is determined in this function and returned for use
% in calling functions

% if number of multiple comparisons are smaller than 10 Bonferroni correction
% will be used automatically
if numel(data) <=10 && strcmp(correctm, 'FDR');
    %cfg.correctm = 'BONF';             % PW this doesn't lead to the use of BONF below, because it changes cfg.correctm not correctm
    correctm = 'BONF';    
    warning('\nTRENTOOL WARNING: Number of channel combinations (%d) to small for FDR -> Bonf was used instead', ...
        numel(data));
end

nr2cmc = size(data,1) * size(data,2) -nrinstmix; %* size(TEresult.TEmat,3);    

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
    ind = 1;
    while ind <= length(significancecor) && significancecor(ind)        
        significance(ind) = 1;
        ind = ind + 1;
    end
    [dummy, unsorteddata] = sort(index);
    significance = significance(unsorteddata);        
    significance = reshape(significance, dim);
    
elseif strcmp(correctm, 'BONF')
        
    significance = (data<=(alpha/nr2cmc));
    
end