function TEsetRandStream

% TESETRANDSTREAM: This function sets the random number stream using the
% appropriate code for the user's MATLAB version. 
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Patricia Wollstadt
% Frankfurt 2014
%


% This function replaces old TRENTOOL code:
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));


v = getVersionString; % get MATLAB version

seed = sum(100*clock);

switch v
    case '7.6'      % code for versions 7.6 and earlier
        rand('twister',seed);
    case '7.7-7.11' % code for versions 7.7 to 7.11
        RandStream('mt19937ar','Seed',seed);
    case '7.12'     % code for versions 7.12 and later
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
end    


function v = getVersionString

v = version;

ind = strfind(v,'.');

major = str2double(v(1:ind(1)-1));
minor = str2double(v(ind(1)+1:ind(2)-1));

if major > 7
    v = '7.12';
else
    if minor >= 12
        v = '7.12';
    elseif minor >= 7 & minor <= 11
        v = '7.7-7.11';
    elseif minor <= 6
        v = '7.6';
    else
        fprintf('\n\n%s\n',v);
        error('TRENTOOL: Unknown MATLAB version!')
    end
end