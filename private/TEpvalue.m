function [p_value TE_diff] = TEpvalue(statistic,numpermutation)

% TEPVALUE: Checks whether the first value of input statistic (e.g.
% TE/MI value of original data) is an extreme value with respect to the 
% distribution of all other values provided in statistic (values for 
% surrogate data). A p-value is computed for this comparison.
%
%
% * REFERENCE INFORMATION
%
%   - permutation test
%       - Maris & Oostenveld (2007). Nonparametric statistical testing of
%         EEG- and MEG-data. J. of Neuroscience Methods, 164, 177-190.
%
% * INPUT PARAMETERS
%
% statistic = distribution of a arbitrary test statistic, where the first
%	      value is assumed to be the value of interest/the empirically
% 	      derived value
% numpermutation = number of permutations that was used to generate the 
%	      surrogate distribution for this channel combination (if 
%	      p = 0, p is set to 1/numpermutation)
%
% * OUTPUT PARAMETERS
%
% p_value = probability of the empirically derived value given the 
%	        distribution provided in the input
% TE_diff = difference between the empirical TE value and the median of the
%           surrogat distribution
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
% 
% Version 1.0 by Michael Wibral, Patricia Wollstadt
% Frankfurt 2013

% CHANGELOG
%
% 2014 - 04 - 05: PW added the manual setting of 1/numpermutation for
%	count = 0 (no sur values larger than the test statistic)

statistic_emp = statistic(1);
statistic_sur = statistic(2:end);

count   = sum(statistic_sur >statistic_emp);
if count > 0
    p_value = count/length(statistic_sur);
else
    p_value = 1/numpermutation;
end

TE_diff = statistic_emp - median(statistic_sur);
