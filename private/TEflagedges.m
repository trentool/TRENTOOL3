function [data, triangle_edges, triangle_nodes] = TEflagedges(data,altpaths,edges4analysis,edges_original)

% TEFLAGEDGES: Flags spurious edges detected by TEdyn and TEbacktracking. 
% Flags used by TEflagedges:
%      2 = cascade effect
%      3 = cascade effect triangle
%      4 = common drive link triangle
%
% 1 denotes instantaneous mixing/volume conduction found by shift/surrogate 
% testing.
%
%
% * INPUT PARAMETERS
%      data        = original data, spurious edges are flagged inside this
%                    structure
%      altpaths    = alternative paths detected by TEdyn and TEbacktracking
%                    [edge_index s t path_count path_tree]
%      edges4analysis  = analyzed edges (original enumeration - enumeration 
%                        is changed within TEgraphanalysis, a requirement 
%                        of the dynamic programming used in TEdyn)
%      edges_original  = all edges in the original data set with original
%                        enumeration
%
%
% * OUTPUT PARAMETERS
%      data       = original data structure with flagged spurious edges.
%                   Flagging is done by setting the following parameters in
%                   the respective columns of .TEpermvalues:
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
%   triangle_edges  = Nx2 list of edges, that form the two potentially 
%                     spurious edges in a triangle: 
%                     [CE_link_ind CD_link_ind].
%   triangle_nodes  = Nx3 list of nodes that form a triangle in the graph,
%                     nodes 1 and 2 form the non-spurious edge, nodes 2 and
%                     3 form the potential CD link, nodes 1 and 3 form the
%                     potential CE link.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.2 by Patricia Wollstadt
% Frankfurt 2012
%
% PW - 08/06/2012

% CHANGELOG
%
% 2012-06-13: PW corrected indexing of flagged edges
% 2012-07-04: PW added documentation and changes command-line feedback
% 2013-11-03: PW corrected triangle count and commented out warning for
%             already flagged paths (this is n/a if ensemble method results
%             are analyzed)
% 2014-03-21: PW added the list of triangles to the output
% 2014-06-23: PW simplified the code/ added comments

n_altpaths = size(altpaths,1);

% make a list from all cascade effects (flagged with a 2)
% list contains: [edge index] [starting node] [target node] [flag type]
flag_list = [cell2mat(altpaths(:,1:3)) 2*ones(n_altpaths,1)];

%% check for triangles and convert edge indices back to original enumeration

no_triangles   = 0;
triangle_edges = [];
triangle_nodes = [];

for i=1:n_altpaths
    
    % find edge index in original data.sgncmb-structure
    current_edge = edges4analysis(altpaths{i,1},:);
    flag_list(i,1) = get_original_ind(current_edge, edges_original);
    
    current_path_list = altpaths{i,5};
       
    for j=1:length(current_path_list)
        
        % check for triangles (alternative paths of length 3: s -> u -> t)
        if (length(current_path_list{j}) == 3);
            
            no_triangles = no_triangles + 1;
            
            % add additional list to flag_list -> simple common drive effect (4)
            CD_link = [current_path_list{j}(2) current_path_list{j}(3)];
            CD_link_ind = get_original_ind(CD_link, edges_original);
            flag_list = cat(1, flag_list, [CD_link_ind CD_link 4]);     % add CD link
            flag_list(i,4) = 3;                                         % change CE link to simple CE/triangle
            
            % remember relevant edges in the triangle for trivariate analysis
            CE_link = [current_path_list{j}(1) current_path_list{j}(3)];
            CE_link_ind = get_original_ind(CE_link, edges_original);
            triangle_edges = cat(1,triangle_edges,[CE_link_ind CD_link_ind]);
            triangle_nodes = cat(1,triangle_nodes,current_path_list{j}');
        end;
    end;
end;


%% flag all spurious edges (cascade and common drive effects)

% check for duplicates in deletion list
duplicates = [];
for i=1:size(flag_list,1)
    ind = flag_list(i,1) == flag_list(i+1:end,1);
    if (sum(ind)>0)
        ind = find(ind)+i;
        duplicates = [duplicates; ind];
    end;
end;
flag_list(duplicates,:) = [];

% delete list by resetting values in TEpermvalues
data.n_spuriousedges = size(flag_list,1);
for i=1:size(flag_list,1)
    ind = flag_list(i,1);
    data.TEpermvalues(ind,:) = [1 0 0 NaN flag_list(i,4) 0];
end;



function new_ind = get_original_ind(edge, original_list)

%   FUNCTION GET_ORIGINAL_IND: finds the index of the edge given as first
%   argument to the function.
%       edge            = 1x2 array, specifying an edge in
%       original list   = Nx2 array of edges
%
% PW 2014-06-23

new_ind = find(original_list(:,1) == edge(1) & original_list(:,2) == edge(2));
