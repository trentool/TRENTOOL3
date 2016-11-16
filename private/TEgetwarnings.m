function msg = TEgetwarnings(cfgTEP, cfgTESS)

% TEGETWARNINGS: This function reads out options set in the cfg structures
% for data preparation and TE estimation. Certain parameter combinations result
% in a warning, which are returned as a string.
%
% * INPUT PARAMETERS
%
%   cfgTEP  - configuration structure passed to TEprepare for data preparation
%   cfgTESS - configuration structure passed to TEsurrogatestats* for TE estimation
%
% * OUTPUT PARAMETERS
%
%   msg - formatted string contatining all raised warnings (concatenated into 
%	  one string
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Patricia Wollstadt, Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt, 2016
%

% CHANGELOG:

msg = '------------------------------------ LOGS & WARNINGS';

%% check if Faes method together with MI calculation is requested
% this analysis returns MI values for channel combinations with volume conduction. When
% a shift-test is performed, all MI values for channels with volume conduction are set to
% NaN; the same is not possible if volume conduction is corrected for by using the 
% Faes method (conditioning on the current sample in the source, x_t).

if strcmp(cfgTEP.extracond, 'Faes_Method') && cfgTESS.MIcalc
	msg = sprintf('%s\n\n# Faes method together with MI calculcation were requested for this analysis. Note that the Faes method corrects for volume conduction in TE estimates only-there is no correction for volume conduction in MI estimates!', msg);
end

%% check if a groupanalysis was performed
% if TE estimation was performed for data prepared for groupanalysis, all estimation
% the embedding dimension was set to be the same for all analyzed subjects (max. over
% all individually optimized embedding)

if cfgTEP.groupanalysis
	msg = sprintf('%s\n\n# TE analysis for group comparison was performed - the embedding dimension was set to the maximum over all channel pairs and subjects to allow for a comparison of estimated TE values between data sets and channels.', msg);
end

msg = sprintf('%s\n\n----------------------------------------------------\n\n', msg);

