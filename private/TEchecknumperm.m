function numperm = TEchecknumperm(cfg, nr2cmc, n_set_1, n_set_2)

% TECHECKNUMPERM checks if a value for the number of permutations for 
% statistical comparison of data sets 1 and 2 was provided by the user 
% (in field cfg.numpermuation). If no value was given, the number of 
% permutations is set such that it exceeds the number of required
% permutations by a factor of 20. 
% The number of required permuations is calculated as 1 / (alpha / nr2cmc), 
% where alpha is the alpha level of the permutations test and nr2cmc is the
% factor for a possible Bonferroni correction (usually the number of 
% channel combinations). I.e., to conduct a permutation test on 100 channel
% combinations at an alpha level of 0.05, with Bonferroni correction, at 
% least 2000 permutations are required (so the smallest possible p-value is
% 1/2000 = 0.05/100).
% 
% If a number of permutations is requested by the user, the function tests
% if it is sufficient to conduct a permutation test with the requested
% alpha level and correction for multiple comparisons. The function further
% tests if the number of trials, which are later permuted, is sufficient to
% allow for the required number of permutations (usually this is the case, 
% 10 trials are sufficient to allow for a permutation test on 600 channel 
% combinations at an alpha level of 0.05).
%
% * INPUT
%   cfg     - configuration structure with fields .alpha, .verbosity and 
%             (optionally) .numpermutation and/or .permstatstype
%   nr2cmc  - factor for a possible Bonferroni correction (usually the 
%             number of channel combinations)
%   n_set_1 - number of units of observation in data set 1 (e.g. no. trials
%             or subjects in a group)
%   n_set_2 - number of units of observation in data set 1 (e.g. no. trials
%             or subjects in a group)
%
%
% * OUTPUT
%   numperm - numbers of permutation, if a value was provided in 
%             cfg.numpermutation and passed all checks, this value is
%             returned; if no value is set, numperm is set according to the
%             number of channels and the requested alpha level
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
% Version 1.0 by Patricia Wollstadt, Michael Wibral
% Frankfurt, 2015

%%
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;

%%

warning('off','all') % nchoosek throws a warning for large n
%max_n = max([n_set_1 n_set_2]);
min_n = min([n_set_1 n_set_2]);

if isfield(cfg, 'permstatstype') && strcmp(cfg.permstatstype, 'depsamplesT')
    possiblepermutations = 2 ^ min_n;
else
    possiblepermutations = nchoosek(n_set_1 + n_set_2, min_n);
end
requiredpermutations = ceil(1/(cfg.alpha/nr2cmc));
warning('on','all')

if ~isfield(cfg, 'numpermutation'),
    %numperm = 190100; % for p<0.01 with a possible bonferroni correcetion of 100, 1/(alpha/nr2cmc)
    numperm = requiredpermutations * 20;
    msg = sprintf(['You didn''t specify a number of permutations. ' ...
        'It was set to %d (for p < %0.2f with a possible bonferroni ', ...
        'correcetion of %d).'], ...
        numperm, cfg.alpha, nr2cmc);
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);    
else
    numperm = cfg.numpermutation;
end
   
if numperm < ceil(1/cfg.alpha)
    error(['\nTRENTOOL ERROR: numperm too small - Nr of permutations' ...
        'must be at least %d!'], num2str(numperm));
end

% PW 2015/10/06: This is too strict for group statistics
% if max_n > 31 && numperm > 2 ^ 31
%     error(['\nTRENTOOL ERROR: numperm too huge - Nr of permutations' ...
%         'should be less than %d!'], 2 ^ 31);    
% elseif numperm > 2 ^ max_n
%     error(['\nTRENTOOL ERROR: numperm too huge - Nr of permutations' ...
%         'should be less than %d!'], 2 ^ max_n);   
% end

if numperm < requiredpermutations
    warning(['\n###### Nr of permutations (%d) not sufficient for correction' ...
        'for multiple comparisons (required: %d)! ######'], numperm, requiredpermutations);
end


if possiblepermutations < requiredpermutations
    fprintf('\n')
    warning(['TRENTOOL ERROR: The number of trials/subjects (%d) is too small to allow ' ...
        'for a sufficient number of permutations given the alpha level (%0.2f) and ' ...
        'necessary corrections for multiple comparisons! (Required permutations: %d, ' ...
        ' possible permutations: %d)'], ...
        min_n, cfg.alpha, requiredpermutations, possiblepermutations);
elseif possiblepermutations < requiredpermutations*20
    fprintf('\n')
    warning(['TRENTOOL WARNING: The number of trials/subjects (%d) probably doesn''t' ...
        'allow for a sufficient number of permutations given the alpha' ...
        'level (%0.2f) and necessary corrections for multiple comparisons! ' ...
        '(Required permutations: %d, possible permutations: %d)'], ...
        min_n, cfg.alpha, requiredpermutations, possiblepermutations);
end
    
    
    