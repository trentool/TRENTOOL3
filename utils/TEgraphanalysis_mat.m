function result = TEgraphanalysis_mat(cfg, mat)

% TEGRAPHANALYSIS_MAT: Calls TEgraphanalysis using a weighted adjacency
% matrix as input (instead of the FieldTrip/TRENTOOL data format).
%
%
% * INPUT PARAMETERS
%
%   cfg: The configuration MUST contain:
%
%       .threshold    = scalar value (in ms): tolerance that is used to 
%                       define the reconstruction interval around w_crit - 
%                       all paths that have a summed weight that falls 
%                       within this interval are considered an alternative 
%                       path
%       .cmc          = tells the function whether to use links that are
%                       significant after correction for multiple
%                       comparison (1) or links that are significant at the
%                       original alpha level (0)
%  
%   mat: adjacency matrix for analysis, where non-zero entries indicate an
%        edge and entries are interpreted as delays associated with edges
%        Note: delays have to have integer values.
%
%
% * OUTPUT
% 
%   result
%       .TEpermvalues = matrix with size edges x 6 (for the exact
%                       specification see INPUT PARAMETERS in the help for
%                       TEgraphanalysis), each entry corresponds to the
%                       channel combination specified in the same row in
%                       the field .TEprepare.channelcombi/
%                       channelcombilabel. If an alternative path was found 
%                       the following changes are made for the respectice 
%                       channelpair:
%                           1 - p-value is set to 1
%                           2 - significance at the prescribed alpha level 
%                               is set to 0
%                           3 - significance after correction for 
%                               multiple comparison is set to 0
%                           4 - mean difference is set to NaN
%                           5 - is set to 2/3/4 according to the type of
%                               spurious interaction:   
%                               2 = cascade effect
%                               3 = cascade effect triangle
%                               4 = common drive link triangle
%                           6 - delay times are set to 0
%
%       .TEprepare = information on channels and channel combinations:
%                    .channellabel  = list of strings, channel names
%                    .channelcombi  = channel combinations analyzed,
%                    (represented by indices for 'channellabel'), indices
%                    also correspond to indices for the input adjacency 
%                    matrix
%                    .channelcombilabel = channel combinations analyzed,
%                    (represented by channel names as strings)
%
%       .graphanalysis   = contains information on the constructed graph as
%                          n_vertices = number of vertices, 
%                          n_edges    = number of edges
%                          density    = graph density, defined as
%                                       dens = E/(V*(V-1))
%                                       V = n_vertices and E = n_edges
%                          cfg        = user provided analysis parameters
%                                       (see INPUT PARAMETERS)
%
%       .n_spurious_edges = number of edges flagged by the algorithm
%       .flagged_mat      = adjacency matrix with flagged links, numbers
%                           indicate type of spurious interaction:
%                               2 = cascade effect
%                               3 = cascade effect triangle
%                               4 = common drive link triangle
%       .TEsteps          = string, list of analysis steps
%
% PW - 04/21/2015
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
% Frankfurt 2015

% CHANGELOG


%% find significant edges

[rows, cols] =  find(mat);

if any(rows == cols)
    error('There are edges connecting nodes to itself, this is not allowed!')
end

n_edges  = size(rows,1);
edges    = zeros(n_edges,6);
sgncmb   = cell(n_edges,2);

for i=1:length(rows)
    edges(i,:)  = [0 ones(1,3) 0 mat(rows(i), cols(i))];
    sgncmb{i,1} = num2str(rows(i));
    sgncmb{i,2} = num2str(cols(i));    
end


%% prepare input in FieldTrip format

data = [];
data.TEpermvalues           = edges;
data.TEprepare.channellabel = unique(sgncmb);
data.TEprepare.channelcombi = [rows cols];
data.TEprepare.channelcombilabel = sgncmb;

%% call graph analysis

result = TEgraphanalysis(cfg, data);

%% change output to adjacency matrix

flagged_edges = find(result.TEpermvalues(:,5));
flagged_mat   = zeros(size(mat));

for i=flagged_edges'
    flagged_mat(rows(i), cols(i)) = result.TEpermvalues(i,5);
end

result.flagged_mat = flagged_mat;